from math import log, pi, sqrt
from scipy.optimize import root_scalar
import numpy as np

import matplotlib.pyplot as plt

# Problem set 9

# Universal constants
R = 287

# Material constants
k = 7/5 #=1.4

# Problem 1
print('--Problem set 9--')

print('--Problem 1--')
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

print('--Problem 2--')
M1 = 2.92
d = 25e-3
T0 = 400
p0 = 29.65e5
R = 287
k = 1.4
f = 0.02

A = pi*d**2/4
a0 = sqrt(k*R*T0)
rho0 = p0/(R*T0)

# same as T_ratio
def isentropic(M,k):
    return 1+1/2*(k-1)*M**2

def a_ratio(M,k):
    return isentropic(M,k)**(1/2)

def rho_ratio(M,k):
    return isentropic(M,k)**(1/(k-1))

def p_ratio(M,k):
    return isentropic(M,k)**(k/(k-1))


print('-2a-')

a1 = a0/a_ratio(M1,k)
rho1 = rho0/rho_ratio(M1,k)
V1 = M1*a1

mdot = rho1*V1*A
print(f'The mass flow rate is {mdot:.5}')


# Alternative procedure
p = p0/p_ratio(M1,k)
print(p)
T = T0/isentropic(M1,k)
print(T)
mdot2 = p*M1*sqrt(k/(R*T))*A

print(f'The mass flow rate is {mdot2:.5} (version 2)')

print('-2b-')
def fanno(M,k):
    return (1-M**2)/(k*M**2) + (k+1)/(2*k)*np.log( (k+1)*M**2/(2 + (k-1)*M**2) )

def T_star_fanno_ratio(M,k):
    return (k+1) / (2 + (k-1)*M**2)

def p_star_fanno_ratio(M,k):
    return 1/M*( (k+1)/(2+(k-1)*M**2) )**(1/2)

fp = fanno(M1,k)
Ls = fp*d/f
print(f'The critcal length is {Ls:.3}')


T1 = T0/isentropic(M1,k)
print(T1)
Ts = T1/T_star_fanno_ratio(M1,k)
print(f'The exit temperature is {Ts:.3}')

p1 = p0/p_ratio(M1,k)
ps = p1/p_star_fanno_ratio(M1,k)
print(f'The exit pressure is {ps:.3}')

cp = R/(1-k)

def dT(M,k,f,D,T,dx):
    return -k*(k-1)*M**4*f*dx/(2*(1-M**2)*D)

def ds(cp,T,dT):
    return cp*np.log( (T+dT)/T )

L = Ls
N = 100

dx = L/(N-1)

T = np.zeros(N)
S = np.zeros(N)

T[0] = T1

for i in range(N-1):
    deltaT = dT(M1,k,f,D,T[i],dx)
    T[i+1] = T[i]
    S[i+1] = S[i] + ds(cp,T[i],deltaT)

#plt.plot(T,S)
#plt.show()

print('--Problem 3--')

def sutherland(mu0,T,T0,S):
    return mu0*(T/T0)**(3/2)*(T0+S)/(T+S)


def darcy(rho,V,d,mu):
    return 64*mu/(rho*V*d)

R = 287
K = 1.4
cp = 1005

mu0 = 1.716e-5
T0 = 273.15
S = 110.4

T = 300
p = 150e3
M1 = 0.4
D = 3e-2
L = 3

rho = p/(R*T)
a = sqrt(k*R*T)
V = M1*a

mu = sutherland(mu0,T,T0,S)
print(f'The viscosity is {mu:.4}')

Re = rho*V*D/mu
print(f'The Reynolds number is {Re:.4}')

#f = darcy(rho,V,d,mu)
f = 0.0136

print(f'The friction factor is {f:.4}')

fp = f * L /D

def eq(M2,M1,fp,k):
    return fanno(M1,k)-fanno(M2,k)-fp

M2 = root_scalar(eq, args=(M1,fp,k), bracket = [M1,1]).root
print(f'The exit mach number is {M2:.3}')
print(fanno(M1,k))
print(fp)
Ms = np.linspace(M1,1,100)
plt.plot(Ms, fanno(M1,k)-fanno(Ms,k)-fp)
plt.show()

print('-ii-')
Lc = D/f*fanno(M1,k)
print(f'The critical length is {Lc:.4}')
