# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 21:41:38 2022

@author: Magnus Lilledahl
"""

import matplotlib.pyplot as plt
import numpy as np

a = 0.5
b = 6
A = 1
N = 100

psi_vals = [1,2,3]
phi_vals = [2,4,6,8]

def stream(psi):
    x = np.linspace(a,b,N)
    y = psi/(A*x)
    return x,y

def equipots(phi):
    x = np.linspace(a,b,N)
    y = np.sqrt(phi-x**2)
    return x,y

for psi in psi_vals:
    x,y = stream(psi)
    plt.plot(x,y)
    
for phi in phi_vals:
    x,y = equipots(phi)
    plt.plot(x,y,'--')
    
plt.show()
ax = plt.gca()
ax.set_aspect('equal')


# Plot velocity field using arrows as representation

plt.figure()
x = np.linspace(0.1,6,10)
y = x
X,Y = np.meshgrid(x,y)
U = A*X
V = -A*Y
plt.quiver(X,Y,U,V)
plt.show()