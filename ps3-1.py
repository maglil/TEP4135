# Script for visualizing flow in problem 4-1 TEP4135

import numpy as np
import matplotlib.pyplot as plt

# Geometry
xmin, xmax = -4, 4
ymin, ymax = -4, 4

N = 2000

x = np.linspace(xmin, xmax, N)
y = np.linspace(ymin, ymax, N)

X,Y = np.meshgrid(x,y)

r = np.sqrt(X**2 + Y**2)
theta = np.arctan2(Y,X)

def stream(r,theta,A, alpha):
    return A*r**alpha*np.sin(theta)

A = 1
alphas = [3, 2, 3/2, 2/3, 1/2]

for alpha in alphas:
    plt.figure()
    plt.contour(X,Y,stream(r,theta,A,alpha),colors='k')
    plt.title(str(alpha))

plt.show()
