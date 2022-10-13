# Script for visualizing flow in problem 4-1 TEP4135

import numpy as np
import matplotlib.pyplot as plt

# Geometry
xmin, xmax = -6, 8
ymin, ymax =  0, 5

N = 2000

x = np.linspace(xmin, xmax, N)
y = np.linspace(ymin, ymax, N)

X,Y = np.meshgrid(x,y)

# Elementary flow solutions
def streamUniform(X,Y,U,V):
    return U*Y-V*X

def streamSource(X,Y,m,x0,y0):
    return m * np.arctan2( (Y-y0),(X-x0) )

# Elementary flow parameters
U,V = .5, 0
m1, x01, y01 = -.5, 0, 1
m2, x02, y02 = -.5, 0, -1

# Total stream function
stream = streamUniform(X,Y,U,V) + \
         streamSource(X,Y,m1,x01,y01) + \
         streamSource(X,Y,m2,x02,y02)
plt.contour(X,Y,stream, levels = 20, colors= 'k' )
plt.show()
