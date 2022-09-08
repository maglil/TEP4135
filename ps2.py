# Problem set 2 TEP4135

from math import pi
import numpy as np
import matplotlib.pyplot as plt

def streamUniform(X,Y,U,V):
    return U*Y-V*X

def streamSource(X,Y,m,x0,y0):
    return m * np.arctan2( (Y-y0),(X-x0) )

# Problem 1
U = 1
V = 0
q = -1.2*2 # Times two since 1.2 is only for half-plane
m = q/(2*pi)
rs = q/(2*pi*U)
print(rs)

# geometry
xmin = 0
xmax = 2
ymin = 0
ymax = 1

N = 500

x = np.linspace(xmin,xmax,N)
y = np.linspace(ymin,ymax,N)

X,Y = np.meshgrid(x,y)
theta = np.arctan2(Y,X)
r = (X**2 + Y**2)**(1/2)
r[r==0] = np.NAN

# Pressure

#p = -q*U*np.cos(theta)/(pi*r)- q**2/(2*pi*r) #(p-p_inf) /rho
pc = -(q*np.cos(theta)/(U*pi*r) + (q/(U*2*pi*r))**2)
pc[pc<-10] = 0


#plt.contourf(X,Y,p)

# streamfunction

stream = streamUniform(X,Y,U,V) + streamSource(X,Y,m,0,0)

fig,ax = plt.subplots()
cs = ax.contourf(X,Y,pc,levels=30)
fig.colorbar(cs)

ax.contour(X,Y,stream,levels=15,colors='k')

plt.show()
