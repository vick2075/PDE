import matplotlib.pyplot as plt
from matplotlib import pyplot, cm
from numpy.linalg import inv
#from thomasalgo import thomas
import numpy as np


'''
Finite Difference Method: Alternating Direction Implicit
'''


# input values
# Ex 1
'''
td = 0.0001          #thermal diffusivity

dt = 20.             #delta t in sec
dx = 0.1             #delta x in cm
dy = 0.1             #delta y in cm

tp = 200             #time period
srx = 1              #spatial range x
sry = 1              #spatial range y 

ic = 200             #initial condition u(x,y,0)
lbc = 100            #left boundary condition u(0,y,t)...T4
rbc = 100            #right boundary condition u(srx=1,y,t)...T2
tbc = 100            #top boundary condition u(x,0,t)...T1
bbc = 100            #bottom boundary condition u(x,sry=1,t)...T3
'''
# Ex 2
'''
td = 0.0001          #thermal diffusivity

dt = 20              #delta t in sec
dx = 0.25            #delta x in cm
dy = 0.1             #delta y in cm

tp = 200             #time period
srx = 10             #spatial range x
sry = 1              #spatial range y

ic = 200             #initial condition u(x,y,0)
lbc = 100            #left boundary condition u(0,y,t)...T4
rbc = 100            #right boundary condition u(srx=1,y,t)...T2
tbc = 100            #top boundary condition u(x,0,t)...T1
bbc = 100            #bottom boundary condition u(x,sry=1,t)...T3
'''
# Ex 3
#'''
td = 0.0001          #thermal diffusivity

dt = 100.             #delta t in sec
dx = 0.25             #delta x in cm
dy = 0.25             #delta y in cm

tp = 1000             #time period
srx = 1              #spatial range x
sry = 1              #spatial range y 

ic = 200             #initial condition u(x,y,0)
lbc = 100            #left boundary condition u(0,y,t)...T4
rbc = 100            #right boundary condition u(srx=1,y,t)...T2
tbc = 100            #top boundary condition u(x,0,t)...T1
bbc = 100            #bottom boundary condition u(x,sry=1,t)...T3
#'''

# process constants
eqx = (td*dt)/(dx**2)# (2*dx**2)#
eqy = (td*dt)/(dy**2)# (2*dy**2)#

pkx = 2*(1+eqx)# (1+2*rx)#
pky = 2*(1+eqy)# (1+2*ry)#
mkx = 2*(1-eqx)# (1-2*rx)#
mky = 2*(1-eqy)# (1-2*ry)#

i = srx/dx+1
j = sry/dy+1
kt = tp/dt

var = (int(j)-2) * (int(i)-2)
mmy = int(j)-2
mmx = int(i)-2


t = np.arange(0, tp+dt, dt)
x = np.arange(0, srx+dx, dx)
y = np.arange(0, sry+dy, dy)




# initial & boundary conditions
T = np.zeros(( int(j), int(i) ))
T[1:-1, 1:-1] = ic
T[0, :] = tbc
T[-1, :] = bbc
T[:, 0] = lbc
T[:, -1] = rbc


# A matrices
A1 = np.zeros((mmx, mmx))
A2 = np.zeros((mmy, mmy))
for jj in range(mmx-1):
    A1[jj, jj+1] = -eqx
    A1[jj, jj] = pkx
    A1[jj+1, jj] = -eqx
    
A1[-1, -1] = pkx

for jj in range(mmy-1):
    A2[jj, jj+1] = -eqy
    A2[jj, jj] = pky
    A2[jj+1, jj] = -eqy

A2[-1, -1] = pky




b2 = np.zeros(mmy)
b1 = np.zeros(mmx)
for k in range(len(t)):
        
    
    for ii in range(1, int(i)+1):
        b2  = eqy*T[2:, ii-1] + mky*T[1:-1, ii-1] + eqy*T[:-2, ii-1]
        b2[0] =  b2[0] + eqy*T[0, ii-1]   # eqy*b2[0] + mky*T[0, ii-1] + eqy*b2[0]
        b2[-1] = b2[-1] + eqy*T[-1, ii-1] # eqy*b2[-1] + mky*T[-1, ii-1] + eqy*b2[-1]

        T[1:-1, ii-1] = np.linalg.solve(A2, b2)
        
    
    for ii in range(1, int(j)+1):
        b1 = eqx*T[ii-1, 2:] + mkx*T[ii-1, 1:-1] + eqx*T[ii-1, :-2]
        b1[0] = b1[0] + eqx*T[ii-1, 0]    # eqx*b1[0] + mkx*T[ii-1, 0] + eqx*b1[0]
        b1[-1] = b1[-1] + eqx*T[ii-1, -1] # eqx*b1[-1] + mkx*T[ii-1, -1] +  eqx*b1[-1]

        T[ii-1, 1:-1] = np.linalg.solve(A1, b1)
        
    










fig = pyplot.figure()
ax = fig.add_subplot(111, projection = '3d')

X, Y = np.meshgrid(x, y)
#surf = ax.plot_surface(X, Y, T, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
ax.plot_surface(X, Y, T, cmap='jet')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()
