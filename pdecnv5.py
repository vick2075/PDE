import numpy as np
from matplotlib import pyplot, cm
import matplotlib.pyplot as plt



'''
Finite Difference Method: Crank Nicolson Method
'''

# Ex 1
#'''
#input values
td = 1.                 #thermal diffusivity

dt = 0.5              #delta t in sec
dx = 0.5                #delta x in cm


i = 100
k = 2000


tp = k*dt               #time period
sr = i*dx                #spatial range

ic = 20.                #initial condition u(x,0)
lbc = 50.                 #left boundary condition u(0,t)
rbc = 80.                 #right boundary condition u(sr,t)
    
#process constants
eq = td*dt/dx**2/2

pk = 2.*(1+eq)-1
mk = 2.*(1-eq)-1


t = np.arange(k)*dt
x = np.arange(i)*dx


var = (i*k)-(k*2)-(i-2)

#'''

# Ex 2

#input values
'''
td = 1.e-4          #thermal diffusivity

dt = 100.            #delta t in sec
dx = 0.2            #delta x in cm

tp = 1000           #time period
sr = 1.            #spatial range

ic = 0            #initial condition u(x,0)
lbc = 100          #left boundary condition u(0,t)
rbc = 100          #right boundary condition u(sr,t)
'''
'''
td = 0.835          #thermal diffusivity

dt = 0.1            #delta t in sec
dx = 2.0            #delta x in cm

tp = 12           #time period
sr = 10            #spatial range

ic = 0            #initial condition u(x,0)
lbc = 100           #left boundary condition u(0,t)
rbc = 50          #right boundary condition u(sr,t)
'''
'''
#process constants
eq = td*dt/dx**2

pk = 2*(1+eq)
mk = 2*(1-eq)

i = sr/dx+1
k = tp/dt+1

var = (i*k)-(k*2)-(i-2)

t = np.arange(0, tp+dt, dt)
x = np.arange(0, sr+dx, dx)

'''





T = np.zeros(( int(k), int(i) ))
T[0,:] = ic
T[:,0] = lbc
T[:,-1] = rbc



A = np.zeros((int(i)-2, int(i)-2))
for jj in range(int(i)-3):
    A[jj, jj+1] = -eq
    A[jj, jj] = pk
    A[jj+1, jj] = -eq

A[-1, -1] = pk



b = np.zeros(int(i)-2)

for ii in range(1, int(k)):
    b = eq*T[ii-1, 2:] + mk*T[ii-1, 1:-1] + eq*T[ii-1, :-2]
    b[0] = b[0] + eq*T[ii, 0]
    b[-1] = b[-1] + eq*T[ii, -1]

    T[ii, 1:-1] = np.linalg.solve(A, b)






fig = pyplot.figure()
ax = fig.add_subplot(111, projection = '3d')


X, Y = np.meshgrid(x,t)
#surf = ax.plot_surface(X, Y, T, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
ax.plot_surface(X, Y, T, cmap='jet')
ax.set_xlabel('X')
ax.set_ylabel('Y-time')
ax.set_zlabel('Z-Temp')
plt.show()









plt.figure()

for jj in range(0, int(k), int(k/4)):
    plt.plot(x, T[jj, :], label= f'{jj*dt} s')
                 




plt.xlabel('distance [m]')
plt.ylabel('Temperature [$\degree$ C]')
plt.grid(True, linestyle='-.')
plt.legend()#[f't = {value} s' for value in t])
plt.show()
