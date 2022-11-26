# Problem set 11, TEP4135

from math import sqrt,tan,pi

g = 9.81

# Problem 1
print("Problem 1")

h = 1.1
dh = 0.12
c_approx = sqrt(g*h)
c = sqrt(g*h*(1+dh/h)*(1+0.5*dh/h))
print(f"The small amplitude wave speed approximation is {c_approx} m/s")
print(f"The wave speed is {c} m/s")
      
dv = c*dh/(h+dh)
print(f"The induced velocity is {dv} m/s")

# Problem 2
print("Problem 2")

h = 0.25
c = sqrt(g*h)
theta = 38/2/180*pi
v = c/tan(theta)
print(f"The velocity is {v} m/s.")

# Problem 3
print("Problem 3")

n = 0.022 #roughness parameter
w = 2
h = 3
theta = 0.85/180*pi

A = w*h
P = 2*h+w
R = A/P
S = tan(theta)

v = 1/n * R**(2/3) * S**(1/2)
Q = v*A
print(f"The volume discharge is {Q:.2} m^3/s")

# Problem 3
print("Problem 5")

S = 0.0045
nu_w = 1e-6 #kinematic visosity water
Re = 500

h = (3*Re*nu_w**2/(g*S))**(1/3)
print(f"The max heigh before turbulence is {h:.2} m")

# Problem 7
print("Problem 7")

n = 0.012
Q = 6
b = 3
h = 1
A = b*h
V = Q/A
R = b*h/(b+2*h)

S = (V*n / R**(2/3))**2
print(f"The design slope needs to be {S:.3}")

n2 = 0.016
V2 = 1/n2*R**(2/3)*S**(1/2)

Q2 = V2*A

dQ = (Q-Q2)/Q

print(f"The percentage change in volume discharge is {dQ:.2%}")

# Problem 8
print("Problem 8")

import numpy as np
import matplotlib.pyplot as plt






areas = [1,5,10]

for A in areas:
    y = np.linspace(0.6,1,300)
    th = np.linspace(0.8,1.3,300)
    r = np.linspace(0.4,1.4,300)
    #Y,TH = np.meshgrid(y,th)
    R,TH = np.meshgrid(r,th)
    #A = 1
    #P = A/Y-Y/np.tan(TH) + 2*Y/np.sin(TH)
    P = (1 + 2*R/np.sin(TH))*np.sqrt(A/(1+R/np.tan(TH))/R)
    
    ir,ith = np.unravel_index(P.argmin(),np.shape(P))
    #iy,ith = np.unravel_index(P.argmin(),np.shape(P))

    #ymin = y[iy]
    rmin = r[ir]
    thmin = th[ith]
    #b = A/ymin-ymin/tan(thmin)
    #rmin = ymin/b
    deg_thmin = np.rad2deg(thmin)
    print(f"ratio is {rmin:.4}, angle is {deg_thmin:.4}, correct is {sqrt(3)/2:.3}")

#fig = plt.figure()
#ax = plt.axes(projection ='3d')
#ax.plot_surface(R, TH, P)

plt.contourf(R,TH,P)
plt.show()





