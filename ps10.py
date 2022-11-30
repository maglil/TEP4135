# Problem set 10
from math import sqrt
from scipy.optimize import root_scalar

## Problem 1
k = 1.4
R = 287
T01 = 600
T02 = 1000
t = T01/T02

M1 = 0.3

def A_ratio(M,k):
    return 1/M*( (2+(k-1)*M**2)/(k+1))**(1/2*(k+1)/(k-1))

def T_ratio_rayleigh(M,k):
    return (k+1)*M**2*(2+(k-1)*M**2) / (1+k*M**2)**2

def T_ratio(M,k):
    return 1 + (k-1)/2*M**2

def eq(M2,M1,k):
    return T_ratio_rayleigh(M1,k)/T_ratio_rayleigh(M2,k)-t

r = A_ratio(M1,k)
print(f"Area ratio is {r}")

M2 = root_scalar(eq, args=(M1,k), x0 = 0.1, x1 = 0.9, bracket = (0.1,1)).root
print(f"M2 is {M2:.5}")

r2 = A_ratio(M2,k)
print(f"Area ratio after reheat is {r2}")


#T1 = T01/T_ratio(M1,k)
#print(f"Initial temperature is {T1:.2}")

#c = sqrt(k*R*T1)
#v1 = c*M1
#print(f"Initial velocity is {v1:.2}")

#def eq(T2,T02):
#    return T_ratio(v1/sqrt(k*R*T2),k)-T02/T2

#T2 = root_scalar(eq, args=(T02), x0=800, x1=2000).root
#print(f"Temperature after reaheat is {T2:.3}")

#M2 = v1/sqrt(k*R*T2)
#print(M2)
#r2 = A_ratio(M2,k)
#print(f"Area ratio is {r2}")
