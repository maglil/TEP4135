# Calculations problem set 12

from math import pi, tan, cos, sin, sqrt, log

# Problem 1
y = 0.8
n = 0.015
S = 1/500
w = 2
theta = 30*pi/180

s = y/sin(theta)
b = y/tan(theta)

P = w + 2*s
A = (w+b)*y

R = A/P
V = 1/n*R**(2/3)*S**(1/2)
Q = V*A
print(Q)

# Problem 4
print("Problem 4")

theta = 25*pi/180
y = 0.35
g = 9.81

Fr = 1/sin(theta)
print(Fr)

c = sqrt(g*y)
v = c*Fr
print(v)
yc = (v**2*y**2/g)**(1/3)
print(yc)

Rh = y
eps = 2.4e-3
f = 2*log(14.8*Rh/eps,10)
S = f/8
print(S)

# Problem 4
print("Problem 4")

y1 = 0.4
y2 = 1.4
r = y2/y1

c1 = sqrt(g*y1)

F1 = sqrt( ((1+2*r)**2-1)/8)
v1 = F1*c1
print(v1)

v2 = v1*y1/y2
print(v2)

yc = (v1*y1/sqrt(g))**(2/3)
print(yc)

hf = (y2-y1)**3/(4*y1*y2)
h1 = y1 + v1**2/(2*g)
ratio = hf/h1
print(ratio)
