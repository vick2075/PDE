import matplotlib.pyplot as plt
from matplotlib import pyplot, cm
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy.integrate import solve_ivp


'''
Method of lines
'''

# Ex 1
'''
#input values
td =  0.835         #thermal diffusivity

dt = 0.1            #delta t in sec
dx = 2              #delta x in cm

tp = 12             #time period
sr = 10             #spatial range

ic = 0              #initial condition u(x,0)
lbc = 100           #left boundary condition u(0,t)
rbc = 50            #right boundary condition u(sr,t)

#process constants
eq = td*dt/dx**2 

pk = 2*(1+eq)
mk = 2*(1-eq)

i = sr/dx+1
k = tp/dt+1
'''

# Ex 2
'''
#input values
td =  0.02          #thermal diffusivity

dt = 0.05           #delta t in sec
dx = 0.01            #delta x in cm

tp = 5.0            #time period
sr = 1.0            #spatial range

ic = 150            #initial condition u(x,0)
lbc = 100           #left boundary condition u(0,t)
rbc = 200           #right boundary condition u(sr,t)

#process constants
eq = td*dt/dx**2 

pk = 2*(1+eq)
mk = 2*(1-eq)

i = sr/dx+1
k = tp/dt+1
'''

# Ex 3
'''
#input values
td =  1.0         #thermal diffusivity

dt = 0.01            #delta t in sec
dx = 0.02              #delta x in cm

tp = 1             #time period
sr = 1             #spatial range

ic = 0              #initial condition u(x,0)
lbc = 100           #left boundary condition u(0,t)
rbc = 0            #right boundary condition u(sr,t)

#process constants
eq = td*dt/dx**2 

pk = 2*(1+eq)
mk = 2*(1-eq)

i = sr/dx+1
k = tp/dt+1
'''

# Ex 4
# Burgers' Eqn
#'''
#input values
td =  0.1 # 0.01 #         #parameter

dt = 0.01            #delta t in sec
dx = 0.02              #delta x in cm

tp = 3             #time period
sr = 1             #spatial range

def ic(x):
    return np.sin(np.pi*x)      #initial condition u(x,0)
lbc = 0                         #left boundary condition u(0,t)
rbc = 0                         #right boundary condition u(sr,t)

#process constants
eq = td*dt/dx**2
#eq1 = ic*dt/dx

pk = 2*(1+eq)
mk = 2*(1-eq)

i = sr/dx+1
k = tp/dt+1
#'''

xx = np.linspace(0, sr, int(i)) 
tt2 = np.linspace(0, tp, int(k))



        








def fu(t, x):
    deriv=[]
    deriv.insert(0, 0)
    
    for ii in range(1,int(i)-1):
        # heat eqn
        #deriv.insert(ii, eq * (x[ii + 1] - 2*x[ii] + x[ii - 1])/dt )

        # burger eqn
        deriv.insert(ii, eq * (x[ii + 1] - 2*x[ii] + x[ii - 1])/dt - \
                     0.5 * (x[ii]*dt/dx) * (x[ii + 1] - x[ii - 1])/dt)

    deriv.insert(int(i)-1, 0)
                
    return deriv








t = 0


x = np.zeros(int(i))

ss = 0
ii = 0
while ss< sr+dx:
    #x[ii] = ic
    x[ii] = ic(ss)
    ss += dx
    ii += 1

x[0] = lbc
x[-1] = rbc






sol = solve_ivp(fu, [t, tp], x,
                method='Radau', t_eval=tt2, atol=1.e-6, rtol=1.e-6)


u = sol.y.T








fig = pyplot.figure()
ax = fig.add_subplot(111, projection = '3d')

X, T = np.meshgrid(xx,tt2)
#surf = ax.plot_surface(X, T, u, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
ax.plot_surface(X, T, u, cmap='jet')
ax.set_xlabel('X')
ax.set_ylabel('Y-time')
ax.set_zlabel('Z-Temp')
plt.show()



plt.figure()

for jj in range(0, int(k), int(k/5)):
    plt.plot(xx, u[jj,:], label='t={0:1.2f}'.format(tt2[jj]))

plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlabel('X position')
plt.ylabel('Temperature')
plt.grid(True, linestyle='-.')
plt.subplots_adjust(top=0.89, right=0.77)
plt.show()



