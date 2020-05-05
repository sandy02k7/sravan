import matplotlib.pyplot as plt
import numpy as np

print("Hello World")

arr = []
time_taken = 0.0;
dr = 0.0004
R = 0.0080
n = int(R/dr)
T = np.ones(n)*1000
Tf = np.ones(n)

ro = 7850.0
h = 24417.98863
bR = 2*h/(dr)
x = (dr*dr*ro)/2.0
y = (dr*dr)
z = (R*dr)
T_water_air = 25.0
for i in range(100000000):
	for j in range(1, n-1):
		cp = 493.71 + (2.3*T[j])
		k =  9.84 + (0.013*T[j])
		dt0 =x*cp/(k)
		dtR = dt0/2.0
		dtr = 2*x*cp/((2*k/y) - (k/z) + (bR))
		dt = min(dt0, min(dtr, dtR))
		a = k/y
		b = k/(2*z)
		cr = ro*cp/(dt) - 2*k/(y)
		d = (dt/(ro*cp))
		Tf[j] = (a*(T[j-1] + T[j+1]) + b*(T[j-1] - T[j+1]) + cr*(T[j]))*d

	cp = 493.71 + (2.3*T[0])
	k =  9.84 + (0.013*T[0])
	dt0 = x*cp/(k)
	dtR = dt0/2.0
	dtr = 2*x*cp/((2*k/y) - (k/z) + (bR))
	dt = min(dt0, min(dtr, dtR))
	a = k/y
	c0 = (ro*cp/(dt)) - (4*k/y)
	d = dt/(ro*cp)
	Tf[0] = (4*a*T[1] + c0*T[0])*d

	cp = 493.71 + (2.3*T[n-1])
	k =  9.84 + (0.013*T[n-1])
	dt0 = x*cp/(k)
	dtR = dt0/2.0
	dtr = 2*x*cp/((2*k/y) - (k/z) + (bR))
	dt = min(dt0, min(dtr, dtR))
	time_taken = time_taken + dt
	aR = 2*k/(y) - k/(z)
	cR = ro*cp/(dt) - (aR) - (bR)
	d = dt/(ro*cp)
	Tf[n-1] = (aR*T[n-2] + bR*T_water_air + cR*T[n-1])*d

	T = Tf
	if(i%1000000 == 0):
		arr.append(Tf[n-1])
		print(Tf[n-1], i/1000000, time_taken)

print(arr)
plt.plot(range(100), arr)
plt.ylabel('Temperature')
plt.show(block=True)
