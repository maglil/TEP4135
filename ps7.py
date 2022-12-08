# Problem set 7 - TEP4135

from math import sqrt
from scipy.optimize import root_scalar

def Mfun(M,k):
    return 1+0.5*(k-1)*M**2

def p_ratio(M,k,ratio=0):
    return Mfun(M,k)**(k/(k-1))-ratio

def critical_area_ratio(M,k,ratio=0):
    return 1/M * ( (1+0.5*(k-1)*M**2) / (0.5*(k+1)) )**(0.5*(k+1)/(k-1))-ratio

def a_ratio(M,k):
    return Mfun(M,k)**(1/2)

def rho_ratio(M,k):
    return Mfun(M,k)**(1/(k-1))



# Problem 1
print("Problem 1")

T0 = 273+20
p0 = 100e5
R = 287
A = 1e-3
A1 = 1e-4
A_star = 5e-4
A2 = 5e-4
k = 1.4

# Problem 1a
print("-1a-")

a0 = sqrt(k*R*T0)
print(f"Stagnation speed of sound is {a0:.4}")
a1 = a0/a_ratio(1,k)
print(f"Critical speed of sound is {a1:.4}")
rho0 = p0/(R*T0)
print(f'Stagnation density is {rho0:.4}')
rho1 = rho0/rho_ratio(1,k)
print(f'Density is {rho1:.4}')
mdot = rho1*a1*A1

rho2 = mdot/(a1*A2)

rho02 = rho2*rho_ratio(1,k)

p1 = rho02*R*T0
print(f'Pressure in tank is {p1:.4}')

p01 = p1

# Problem 1b
p2 = p01*0.9395
print(p2)
p_star = p01/p_ratio(1,k)
print(p_star)

print("Problem 1a")

M = root_scalar(critical_area_ratio, args=(k,A/A_star), method='bisect', bracket=[1,4]).root
print(f"Mach number at exit for chocked, non-shock flow (design), is {M:.3}")
p_design = p01/p_ratio(M,k)
print(f"The design pressure is {p_design:.3}")

M = root_scalar(critical_area_ratio, args=(k,A/A_star), method='bisect', bracket=[.1,1]).root
print(f"Mach number at exit for highest exit pressure whil still critical, is {M:.3}")
p_design = p01/p_ratio(M,k)
print(f"The exit pressure at this Mach number is {p_design:.3}")
