from math import sqrt
from scipy.optimize import root_scalar

## Problem 1

# Physical quantities, upstream
M1 = 2.4
p1 = 130e3
T0 = 450
k = 1.4
R = 287
cp = 1005

# upstream temperature
T1 = T0/(1+1/2*M1**2*k*R/cp)
print("T1 = " + str(T1))

T1 = T0/(1 + (k-1)/2*M1**2)
print("T1' = " + str(T1))

# Upstream speed of sound
a1 = sqrt(k*R*T1)
print("a = " + str(a1))

# Upstream velocity
v1 = M1*a1
print("v1 = " + str(v1))

# Upstream rho
rho1 = p1/(R*T1)
print("rho1 = " + str(rho1))

# Upstream stagnation pressure
p01 = p1 + v1**2*rho1
print("p01' = " + str(p01))

p01 = p1 * (1 + 1/2*(k-1)*M1**2)**(k/(k-1))
print("p01 = " + str(p01))

print( 1/(1 + 1/2*(k-1)*M1**2)**(k/(k-1)) )
print(p1/p01)

# Downstream stagnation pressure
p02 = p01*((k+1)*M1**2/(2+(k-1)*M1**2))**(k/(k-1))*((k+1)/(2*k*M1**2-(k-1)))**(1/(k-1))
print("p02 = " + str(p02))

# Stagnation pressure loss
dp = p02-p01
print("dp = " + str(dp))

# Downstream mach number
M2 = (((k-1)*M1**2+2)/(2*k*M1**2-(k-1)))**(1/2)
print("M2 = " + str(M2))

# Downstream temperature
T2 = T0/(1+1/2*M2**2*k*R/cp)
print("T2 = " + str(T2))

# Upstream temperature
T1 = T2/( (2+(k-1)*M1**2)*(2*k*M1**2-(k-1))/((k+1)**2*M1**2) )
print(" T1 = " + str(T1))

# Downstream velocity

v2 = M2*sqrt(k*R*T2)
print("v2 = " + str(v2))

# Velocity ratio
print(v1/v2)
print(((k+1)*M1**2)/((k-1)*M1**2+2))

## Problem 2
p02 = 70
p1 = 15
def p_ratio(M1):
    return (1+0.5*(k-1)*M1**2)**(k/(k-1)) * (1/(k+1)*(2*k*M1**2-(k-1))) - p02/p1

M1 = root_scalar(p_ratio,method= 'brentq', bracket=[0,5]).root
print(M1)

## Problem 3
T = 282
T0 = 300

M2sq = 2*(T0/T-1)/((k-1))
print("M2 = " + str(sqrt(M2sq)))

def mach(M):
    return ((k-1)*M**2 + 2)/(2*k*M**2-(k-1)) - M2sq

M1 = root_scalar(mach, method = 'brentq', bracket = [.5,5]).root
print(M1)

def T_ratio_shock(M1):
    return (2 + (k-1)*M1**2)*(2*k*M1**2-(k-1))/((k+1)**2*M1**2)

T1 = T2/T_ratio_shock(M1)
print("T1 = " + str(T1))

v1 = M1*sqrt(k*R*T1)
print("v1 = " + str(v1))
