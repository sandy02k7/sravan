import datetime
import matplotlib.pyplot as plt
import numpy as np

def Cp(T):
        return 470 + (0.02*T) + (0.00038*T*T)
    
def K(T):
        return 52.27 - (0.01541*T) - (0.00002155*T*T)

def transpose(l1, l2): 
  
    # we have nested loops in comprehensions 
    # value of i is assigned using inner loop 
    # then value of item is directed by row[i] 
    # and appended to l2 
    l2 =[[row[i] for row in l1] for i in range(len(l1[0]))] 
    return l2 

arr = []
arr_t = []
dr = 0.0004
R = 0.0080
n = int(R/dr)
time_taken = np.ones(n)*0;
T = np.ones(n)*1000
Tf = np.ones(n)

ro = 7850.0
h = 24417.98863
bR = 2*h/(dr)
x = (dr*dr*ro)/2.0
y = (dr*dr)
z = (R*dr)
T_water_air = 25.0
dt = 1e9

for t in range(20, 1200):
	dt0 = x*Cp(t)/(K(t))
	dtR = dt0/2.0
	dtr = 2*x*Cp(t)/((2*K(t)/y) - (K(t)/z) + (bR))
	dt = min(dt, min(dt0, min(dtr, dtR)))

itr = int(1e8)
for i in range(itr):

	vec = []

	cp0 = Cp(T[0])
	k0 = K(T[0])
	time_taken[0] = time_taken[0] + dt
	a = k0/y
	c0 = (ro*cp0/(dt)) - (4*k0/y)
	d = dt/(ro*cp0)
	Tf[0] = (4*a*T[1] + c0*T[0])*d
	if(i%(itr/100) == 0):
		vec.append(tuple((Tf[0], time_taken[0])))

	for j in range(1, n-1):
		cpr = Cp(T[j])
		kr = K(T[j])
		time_taken[j] = time_taken[j] + dt
		a = kr/y
		b = kr/(2*z)
		cr = ro*cpr/(dt) - 2*kr/(y)
		d = (dt/(ro*cpr))
		Tf[j] = (a*(T[j-1] + T[j+1]) + b*(T[j-1] - T[j+1]) + cr*(T[j]))*d
		if(i%(itr/100) == 0):
			vec.append(tuple((Tf[j], time_taken[j])))

	cpR = Cp(T[n-1])
	kR = K(T[n-1])
	time_taken[n-1] = time_taken[n-1] + dt
	aR = 2*kR/(y) - kR/(z)
	cR = ro*cpR/(dt) - (aR) - (bR)
	d = dt/(ro*cpR)
	Tf[n-1] = (aR*T[n-2] + bR*T_water_air + cR*T[n-1])*d
	if(i%(itr/100) == 0):
		vec.append(tuple((Tf[n-1], time_taken[n-1])))

	T = Tf

	if(i%(itr/100) == 0):
		arr.append(vec);
		now = datetime.datetime.now()
		print(Tf[n-1], i/(itr/100), time_taken[n-1], now.strftime("%Y-%m-%d %H:%M:%S"))

filename = 'data.txt'
with open(filename, 'w') as f:
	arr_t = transpose(arr, arr_t)
	for ind in range(len(arr_t)):
		print(arr_t[ind], file=f)  # Python 3.x

#plt.plot(range(len(arr)), arr)
#plt.ylabel('Temperature')
#plt.show(block=True)
