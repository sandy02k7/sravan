import matplotlib.pyplot as plt
import numpy as np

def transpose(l1, l2): 
  
    # we have nested loops in comprehensions 
    # value of i is assigned using inner loop 
    # then value of item is directed by row[i] 
    # and appended to l2 
    l2 =[[row[i] for row in l1] for i in range(len(l1[0]))] 
    return l2 


arr = []
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
mcp = 0
ccp = 0
mk  = 0
ck  = 0
itr = int(1e9)
for i in range(itr):

	vec = []

	if(T[0] > 802):
		mcp = 2.3
		ccp = 493.71
		mk  = 0.013
		ck  = 9.84
	elif(T[0] <= 802 and T[0] > 677):
		mcp = 0
		ccp = 1400
		mk  = -0.047
		ck  = 75.42
	elif(T[0] <= 677):
		mcp = 0.5066
		ccp = 281.4
		mk  = -0.047
		ck  = 75.42
	cp = ccp + (mcp*T[0])
	k =  ck + (mk*T[0])
	dt0 = x*cp/(k)
	dtR = dt0/2.0
	dtr = 2*x*cp/((2*k/y) - (k/z) + (bR))
	dt = min(dt0, min(dtr, dtR))
	time_taken[0] = time_taken[0] + dt
	a = k/y
	c0 = (ro*cp/(dt)) - (4*k/y)
	d = dt/(ro*cp)
	Tf[0] = (4*a*T[1] + c0*T[0])*d
	if(i%(itr/100) == 0):
		vec.append(tuple((Tf[0], time_taken[0])))

	for j in range(1, n-1):
		if(T[j] > 802):
			mcp = 2.3
			ccp = 493.71
			mk  = 0.013
			ck  = 9.84
		elif(T[j] <= 802 and T[j] > 677):
			mcp = 0
			ccp = 1400
			mk  = -0.047
			ck  = 75.42
		elif(T[j] <= 677):
			mcp = 0.5066
			ccp = 281.4
			mk  = -0.047
			ck  = 75.42
		cp = ccp + (mcp*T[j])
		k =  ck + (mk*T[j])
		dt0 =x*cp/(k)
		dtR = dt0/2.0
		dtr = 2*x*cp/((2*k/y) - (k/z) + (bR))
		dt = min(dt0, min(dtr, dtR))
		time_taken[j] = time_taken[j] + dt
		a = k/y
		b = k/(2*z)
		cr = ro*cp/(dt) - 2*k/(y)
		d = (dt/(ro*cp))
		Tf[j] = (a*(T[j-1] + T[j+1]) + b*(T[j-1] - T[j+1]) + cr*(T[j]))*d
		if(i%(itr/100) == 0):
			vec.append(tuple((Tf[j], time_taken[j])))

	if(T[n-1] > 802):
		mcp = 2.3
		ccp = 493.71
		mk  = 0.013
		ck  = 9.84
	elif(T[n-1] <= 802 and T[n-1] > 677):
		mcp = 0
		ccp = 1400
		mk  = -0.047
		ck  = 75.42
	elif(T[n-1] <= 677):
		mcp = 0.5066
		ccp = 281.4
		mk  = -0.047
		ck  = 75.42
	cp = ccp + (mcp*T[n-1])
	k =  ck + (mk*T[n-1])
	dt0 = x*cp/(k)
	dtR = dt0/2.0
	dtr = 2*x*cp/((2*k/y) - (k/z) + (bR))
	dt = min(dt0, min(dtr, dtR))
	time_taken[n-1] = time_taken[n-1] + dt
	aR = 2*k/(y) - k/(z)
	cR = ro*cp/(dt) - (aR) - (bR)
	d = dt/(ro*cp)
	Tf[n-1] = (aR*T[n-2] + bR*T_water_air + cR*T[n-1])*d
	if(i%(itr/100) == 0):
		vec.append(tuple((Tf[n-1], time_taken[n-1])))

	T = Tf

	if(i%(itr/100) == 0):
		arr.append(vec);
		print(Tf[n-1], i/(itr/100), time_taken[n-1])

filename = 'data.txt'
with open(filename, 'w') as f:
	arr_t = []
	arr_t = transpose(arr, arr_t)
	for ind in range(len(arr_t)):
		print(arr_t[ind], file=f)  # Python 3.x
#plt.plot(range(100), arr)
#plt.ylabel('Temperature')
#plt.show(block=True)
