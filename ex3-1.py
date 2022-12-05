import numpy as np
import matplotlib.pyplot as plt

xmin = -2
xmax = 2
ymin = xmin
ymax = xmax

N = 100

x = np.linspace(xmin, xmax, N)
y = np.linspace(ymin, ymax, N)

X,Y = np.meshgrid(x,y)

def stream(X,Y,A,alpha):
    R = np.sqrt(X**2 + Y **2)
    theta = (np.arctan2(Y,X) + 2*np.pi) % (2*np.pi) 
    return A*R**alpha*np.sin(alpha*theta)

A = 1
alpha = 0.7
plt.contour(X,Y,stream(X,Y,A,alpha),levels=20,colors='k')
plt.show()
