from scipy.optimize import root_scalar
from math import tan,pi,cos,sin

# Example 10.6
def f(y):
    theta_d = 50
    theta = pi*theta_d/180
    Q = 16
    g = 9.81
    
    return (2*y*Q**2 / (tan(theta)*g))**0.5 - y**2/tan(theta)

res = root_scalar(f,x0=1,x1=2)
print(res.root)

theta_d = 50
theta = pi*theta_d/180
Q = 16
g = 9.81
n = 0.018

y = (2*Q**2*tan(theta)**2/g)**(1/5)
print(y)
V = Q*tan(theta)/y**2
print(V)
A = y**2/tan(theta)
P = 2*y/sin(theta)
R = A/P
S = n**2*V**2/R**(4/3)
print(S)
print(S*180/pi)
b = 2*y/tan(theta)
S = g*n**2*P/(R**(1/3)*b)
print(S)
