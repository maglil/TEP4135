import matplotlib.pyplot as plt
import numpy as np
from math import pi

d = 1
U = 1
B = 1
y = np.linspace(-1*d, 1*d, 100)
x = np.linspace(-3*d, 3*d, 100)
X,Y = np.meshgrid(x,y)

def uniform_x(U,x,y):
    return U*y

def doublet(B,x0,y0,x,y):
    theta = np.arctan2(y-y0,x-x0)
    r = np.sqrt((x-x0)**2 + (y-y0)**2)
    return -B/r*np.sin(theta)

psi = uniform_x(U,X,Y) + doublet(B,0,d,X,Y) + doublet(B,0,-d,X,Y)

plt.contour(psi,X,Y)
#plt.show()

a = 0.5
k = 1.4
M = 2.8

def a_ratio_squared(M,k):
    return 1/M**2*(2/(k+1)*(1+(k-1)/2*M**2))**((k+1)/(k-1))

ac = a/np.sqrt(a_ratio_squared(M,k))
print(ac)

r = 0.5
A = pi*r**2
Q = 1.5
S = 0.02

n = A/Q*(r/2)**(2/3)*S**(1/2)
print(n)
