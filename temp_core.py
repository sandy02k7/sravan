import datetime
import matplotlib.pyplot as plt
import numpy as np

def Cp(T):
        return 600
#        return 470 + (0.02*T) + (0.00038*T*T)
    
def K(T):
		return 45

#    if(T <= 800):
#           return 54 - 0.0333*T
#  else:
#            return 27.3 
#        return 52.27 - (0.01541*T) - (0.00002155*T*T)

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
Rc = 0.014 
Qtot = 130
v_bar = 6
no = 7
Qno = Qtot/no
p = 2*3.14*((Rc + R)/2)
A = 3.14*(Rc*Rc - R*R)
n = int(R/dr)
time_taken = 0
time_taken_arr = []
T = np.ones(n)*1000
Tf = np.ones(n)
Dh = 4*A/p
v_w = Qno/(3600*A)
v_rel = abs(v_w - v_bar)
ro_w = 996
mu =   0.000993
Re = (v_rel*ro_w*Dh)/mu
kw = 0.597
pr = 7.0 
h = 0
length_cl = 6.87

if(Re < 5000):
	h = (1.86*(Re*pr)**(0.33)*kw)/Dh
else:
	h = (0.023*pr**(0.3)*Re**(0.8))/Dh

ro = 7850
abR = 2*11/(dr)
wbR = 2*h/(dr)
x = (dr*dr*ro)/2.0
y = (dr*dr)
z = (R*dr)
T_water_air = 25.0
time_thr = length_cl/v_bar
dt = 1e9

for t in range(225, 1000):
	dt0 = x*Cp(t)/(K(t))
	dtR = dt0/2.0
	dtr = 2*x*Cp(t)/((2*K(t)/y) - (K(t)/z) + (wbR))
	dt = min(dt, min(dt0, min(dtr, dtR)))

itr = int(1e3)
for i in range(itr):

	vec = []

	cp0 = Cp(T[0])
	k0 = K(T[0])
	a = k0/y
	c0 = (ro*cp0/(dt)) - (4*k0/y)
	d = dt/(ro*cp0)
	Tf[0] = (4*a*T[1] + c0*T[0])*d
	if(i%(itr/100) == 0):
		vec.append(Tf[0])

	for j in range(1, n-1):
		cpr = Cp(T[j])
		kr = K(T[j])
		a = kr/y
		b = kr/(2*z)
		cr = ro*cpr/(dt) - 2*kr/(y)
		d = (dt/(ro*cpr))
		Tf[j] = (a*(T[j-1] + T[j+1]) + b*(T[j-1] - T[j+1]) + cr*(T[j]))*d
		if(i%(itr/100) == 0):
			vec.append(Tf[j])

	if(time_taken < time_thr):
		bR = wbR
	else:
		bR = abR
	cpR = Cp(T[n-1])
	kR = K(T[n-1])
	aR = 2*kR/(y) - kR/(z)
	cR = ro*cpR/(dt) - (aR) - (bR)
	d = dt/(ro*cpR)
	Tf[n-1] = (aR*T[n-2] + bR*T_water_air + cR*T[n-1])*d
	if(i%(itr/100) == 0):
		vec.append(Tf[n-1])

	T = Tf
	time_taken  = time_taken + dt

	if(i%(itr/100) == 0):
		arr.append(vec);
		time_taken_arr.append(time_taken)
		now = datetime.datetime.now()
		print(Tf[n-1], i/(itr/100), time_taken, now.strftime("%Y-%m-%d %H:%M:%S"))

filename = 'data.csv'
with open(filename, 'w') as f:
	for ind in range(len(arr)):
		print(arr[ind], file=f)  # Python 3.x
	print(time_taken_arr, file=f)  # Python 3.x

plt.plot(range(len(arr[n-1])), arr[n-1])
plt.xlabel('time')
plt.ylabel('Temperature')
plt.show(block=True)
