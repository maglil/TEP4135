import sympy as sym
from math import sqrt, asin, atan, sin, pi, tan, cos
from scipy.optimize import root_scalar

import numpy as np
import matplotlib.pyplot as plt

# Example 9.1

(p, p1, p2, T, 
T1, T2, gamma,
rho, R, dH, cp, cv) = sym.symbols(
    'p p1 p2 T T1 T2 gamma rho R dH, cp, cv')

#Ideal gas law
idealgas = sym.Eq(p , rho*R*T)

# Change in entropy
entropychange = sym.Eq(dH, cp*(T2-T1))

#
p1 = 1.7e6
p2 = 248e3
rho0 = 18
gamma_ideal = 1.4
R_argon = 208
T2e = 400
gamma_Ar = 1.67
cp_Ar = R_argon*gamma_Ar/(gamma_Ar-1)

# Initial values
initvals = {p : p1,
            rho : rho0,
            R : R_argon
            }
# Final values
endvals = {p : p2,
           T : T2,
           R : R_argon
           }

# a) 
T1e = sym.solve(idealgas,T)[0].subs(initvals)
print(T1e)

# b)
rho2e = sym.solve(idealgas,rho)[0].subs(endvals)
print(rho2e)

#c)
res = sym.solve(entropychange,dH)[0].subs({T1:T1e,T2:T2e, cp:cp_Ar})
print(res)

vals = {T2 : 400,
        p1 : 1.7e6,
        p2 : 248e3,
        gamma : 1.67,
        R : 208,
        rho : 18
    }



 
# Obs! Flow is not isentropic
#isentropic = sym.Eq( p2/p1, (T2/T1)**(gamma/(gamma-1)) )
#Ti = sym.solve(isentropic,T1)[0]

# Example 9.8
print("Example 9.8")

k = 1.4
R = 287
T0 = 400
p0 = 120e3

rho0 = p0/R/T0
print(f"Density is {rho0:.3}")

A0 = 6e-4
a0 = sqrt(k*R*T0)
ratio_p_critical = (2/(k+1))**(k/(k-1))
print(f"ratio is {ratio_p_critical:.3}")



p_star = p0*ratio_p_critical
print(f"P_star is  {p_star:.3}")

# Subsonic=90, supersonic= 63 
p = 63e3
# let e = 1 + 1/2*(k-1)*M**2

e = (p0/p)**((k-1)/k)
a = a0/sqrt(e)
M = sqrt(2*(e-1)/(k-1))
rho = rho0/(e**(1/(k-1)))
print(f"Mach number is {M:.3}")
print(f"Sound velocity is {a}")
print(f"Density is {rho:3}")
V = M*a
mdot = rho*V*A0
print(f"mass flow is {mdot:.3}")

# Example 9.9
print("Example 9.9")

def area_ratio(M,k,ratio):
    r = 1/M * ( (1+0.5*(k-1)*M**2) / (0.5*(k+1)) )**(0.5*(k+1)/(k-1))-ratio
    return r

def a_ratio(M,k):
    return sqrt(1 + 0.5*(k-1)*M**2)

def rho_ratio(M,k):
    return (1+0.5*(k-1)*M**2)**(1/(k-1))

def p_ratio(M,k):
    return (1+0.5*(k-1)*M**2)**(k/(k-1))

T0 = 500
p0 = 1e6
Ae = 0.008
A_star = 0.002
rho0 = p0/(R*T0)
a0 = sqrt(k*R*T0)

M = root_scalar(area_ratio, args=(k,Ae/A_star),method='bisect',bracket=[1,4]).root
print(f"Mach number is {M:.4}")

ae = a0/a_ratio(M,k)
print(f"Sound velocity is {ae:.4}")

rhoe = rho0/rho_ratio(M,k)
print(f"density is {rhoe:.4}")

Ve = M*ae

pe = p0/p_ratio(M,k)
print(f"Pressure is {pe:.4}")

md = rhoe*Ve*Ae
print(f"Mass flow is {md:.4}")

def p_star_ratio(k):
    return (2/(k+1))**(k/(k-1))

M2 = root_scalar(area_ratio, args=(k,Ae/A_star),method='bisect',bracket=[0.1,1]).root
print(f"Mach number is {M2:.4}")

p_star = p_star_ratio(k)*p0

p2 = p0/p_ratio(M2,k)
print(f"pressure at critical is {p2:.4}")

k = 1.4
m = np.linspace(0.1,4,100)
r = area_ratio(m,k,0)
#plt.plot(m,r)
#plt.show()

# Example 9.16
h = 5
l = 9
theta = atan(h/l)
M = 1/sin(theta)
print(f'The mach number is {M:.2}')

# Example 9.17
print('--Example 9.17--')

p1 = 69e3
M1 = 2
theta = 10*pi/180
k = 1.4
eps = 0.1

def beta_eq(beta, M1, k, theta):
    return (k+1)*M1**2*np.sin(beta)**2 / ((k-1)*M1**2*np.sin(beta)**2 + 2) - np.tan(beta)/np.tan(beta-theta)

beta = root_scalar(beta_eq, args=(M1, k, theta), x0=.1, x1=pi/4, bracket=[theta+eps, pi/4]).root
beta_deg = beta*180/pi
print(f'Wave angle beta is {beta_deg:.3}')

#b = np.linspace(theta+eps,pi/4,100)
#plt.plot(b, beta_eq(b,M1,k,theta))
#plt.show()

M1n = M1*sin(beta)
M2n = sqrt( ((k-1)*M1n**2+2) / (2*k*M1n**2-(k-1)) )
M2 = M2n/sin(beta-theta)
print(f'Mach number after obliqe shock is {M2:.4}')

p2 = 1/(k+1)*(2*k*M1**2*sin(beta)**2-(k-1))*p1
print(f'Pressure after shock is {p2:.3}')

p2 = k*M1**2*tan(theta)/sqrt(M1**2-1)*p1 + p1
print(f'Linearized approximated pressure after shock is {p2:.3}')

mu = asin(1/M1)
beta = asin(sin(mu) + (k+1)/(4*cos(mu))*tan(theta))
beta_deg = beta*180/pi
print(f'Approximate beta is {beta_deg:.3}')
