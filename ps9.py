from math import log
from scipy.optimize import root_scalar
import numpy as np

import matplotlib.pyplot as plt

# Problem set 9

# Universal constants
R = 287

# Material constants
k = 7/5 #=1.4

# Problem 1

# System properties
f = 0.0357
D = 0.200

# Physical quantities, inital conditions
v_0 = 180
T_0 = 350
p_0 = 250

# Derived physical properties, initial conditions
a_0 = (k*R*T_0)**(1/2)
M_0 = v_0/a_0
M2_0 = M_0**2 # M2 is M**2
nu = 1-5e-5
Re = v_0*D/nu
print("Re = " + str(Re))
print("M_0 = " + str(M_0))

# Solution positition
L = 5.00

### Solve by integrating differential equations ###
print("\nNumerical solution through integration...")

# Solver parameters
N = 1000
dx = L/(N-1);

# Initialize arrays
m2 = np.zeros(N)
p = np.zeros(N)
T = np.zeros(N)
v = np.zeros(N)
rho = np.zeros(N)
x = np.zeros(N)

# Inital conditions
m2[0] = M2_0
p[0] = p_0
T[0] = T_0
v[0] = v_0
rho[0] = p_0/(R*T_0)

# Differentials
def dfannoM(m2,dx,k,D):
    return k*m2*(1+0.5*(k-1)*m2)/(1-m2)*f*dx/D*m2

def dfannop(p,m2,dx,k,D):
    return -k*m2*(1+(k-1)*m2)/(2*(1-m2))*f*dx/D*p

def dfannoT(T,m2,dx,k,D):
    return -(k*(k-1)*m2**2)/(2*(1-m2))*f*dx/D*T

def dfannov(v,m2,dx,k,D):
    return k*m2/(2*(1-m2))*f*dx/D*v

def dfannorho(rho,m2,dx,k,D):
    return -k*m2/(2*(1-m2))*f*dx/D*rho

for i in range(N-1):
    m2[i+1] = m2[i] + dfannoM(m2[i],dx,k,D)
    p[i+1] = p[i] + dfannop(p[i],m2[i],dx,k,D)
    T[i+1]= T[i] + dfannoT(T[i],m2[i],dx,k,D)
    v[i+1] = v[i] + dfannov(v[i],m2[i],dx,k,D)
    rho[i+1] = rho[i] + dfannorho(rho[i],m2[i],dx,k,D)
    x[i+1] = x[i] + dx

print("m = ", (m2[-1])**0.5)
print("p = ", p[-1])
print("T = ", T[-1])
print("v = ", v[-1])

#plt.plot(x,m2**0.5)
#plt.figure()
#plt.plot(T)
#plt.show()

dx = 5

M2 = M2_0
T = T_0
v = v_0

### Analytical solution with numerical root finding ###
print("\nAnalytical + Numerical solution through root finding")

# Fanno parameter, integrated differential equation for M**2 to sonic point M = 1
# Assumin constant average Darcy friction factor f
fp = (1-M2)/(k*M2) + (k+1)/(2*k)*log((k+1)*M2/(2+(k-1)*M2))
print("fp = " + str(fp))

# Sonic length, L^star
Ls = fp*D/f
print("Ls = " + str(Ls))

# Fanno parameter at length L from sonic point
l = 2
fp_l = f*l/D
print("fp_l = " + str(fp_l))

# For finding root(M2) for some specific fanno parameter fp
def fanno(M2):
    k = 1.4
    fp = fp_l
    return (1-M2)/(k*M2) + (k+1)/(2*k)*log((k+1)*M2/(2+(k-1)*M2)) - fp

M2_l = root_scalar(fanno, method="brentq", bracket=[0.1,2]).root
print("M_l = " + str(M2_l**0.5))

# Integrated expressions from (q,M) to (q*,M=1) for any quantity q

def p_ratio(m2):
    return (1/m2*(k+1)/(2+(k-1)*m2))**0.5

def v_ratio(m2):
    return (1/m2*(2+(k-1)*m2)/(k+1))**(-0.5)

def T_ratio(m2):
    return (k+1)/(2+(k-1)*m2)

p_l = p_ratio(M2_l)/p_ratio(M2_0)*p_0
v_l = v_ratio(M2_l)/v_ratio(M2_0)*v_0
T_l = T_ratio(M2_l)/T_ratio(M2_0)*T_0
print("p_l = "+ str(p_l))
print("v_l = "+ str(v_l))
print("T_l = "+ str(T_l))

