#script for exploring the Rankine oval

import numpy as np
import matplotlib.pyplot as plt

xmin = -4
xmax = 4
ymin = xmin
ymax = xmax

N = 100

x = np.linspace(xmin, xmax, N)
y = np.linspace(ymin, ymax, N)

X,Y = np.meshgrid(x,y)

def streamUniform(X,Y,U,V):
    return U*Y-V*X

def streamSource(X,Y,m,x0,y0):
    return m * np.arctan2( (Y-y0),(X-x0) )

U = .5
V = 0
m = 1
a = .5
stream = streamUniform(X,Y,U,V) + \
        streamSource(X,Y,m,-a,0) + \
        streamSource(X,Y,-m,a,0)
plt.contour(X,Y,stream, levels = 20, colors = 'k')
#plt.contour(X,Y,streamSource(X,Y,m,1,1) , levels = 30, colors = 'k')

plt.show()
