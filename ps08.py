import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar
from math import pi, asin, sin, sqrt, tan, atan, cos

k = 1.4
K = (k+1)/(k-1)
M1 = 1.5
theta = 20*pi/180
p = 80
T = 273-20

def omega(M):
    return np.sqrt(K)*np.arctan( ((M**2-1)/K)**(1/2) ) - np.arctan( (M**2-1)**(1/2))

def eq(M2, M1):
    return omega(M2)-omega(M1)-theta

M2 = root_scalar(eq, args=(M1), x0=1.1, x1=10, bracket=(1.1,10)).root
print(M2)

def p_ratio(M,k):
    return (1+1/2*(k-1)*M**2)**(k/(k-1))

def T_ratio(M,k):
    return 1 + (k-1)/2*M**2

p0 = p*p_ratio(M1,k)
print(p0)

p2 = p0/p_ratio(M2,k)
print(p2)

T0 = T*T_ratio(M1,k)
T2 = T0/T_ratio(M2,k)
print(T2)

mu1 = asin(1/M1)*180/pi
print(mu1)

mu2 = asin(1/M2)*180/pi
print(mu2-20)

# Problem 2
print('--Problem 2--')

# Given quantities
mu = 30*pi/180 # angles in radians
b = 50*pi/180

M1 = 1/sin(mu)
print(f'The incoming Mach number is {M1:.3}')

Mn1 = M1*sin(b)

Mn2 = sqrt( ((k-1)*Mn1**2+2) / (2*k*Mn1**2 - (k-1) ) )
print(f'Normal component of Mach number after shock is {Mn2:.3}')

theta = atan (2/tan(b)*(M1**2*sin(b)**2-1) / (M1**2*(k+cos(2*b)) + 2))
theta_deg = theta*180/pi
print(f'The inclincation is {theta_deg:.3} degrees')

# Alternative way finding roots
def rho_ratio_theta(beta, theta):
    return np.tan(beta)/np.tan(beta-theta)

def rho_ratio_mach(M,k,beta):
    return (k+1)*M**2*np.sin(beta)**2 / ( (k-1) *M**2*np.sin(beta)**2 + 2 )

def eq(theta,M,k,beta):
    return rho_ratio_mach(M,k,beta)-rho_ratio_theta(beta,theta)

#theta = np.linspace(0,b/2,100)
#plt.plot(theta,eq(theta,M1,k,b))
#plt.show()

theta = root_scalar(eq,args=(M1,k,b), bracket=[0,b/2]).root
theta_deg = theta*180/pi
print(f'The inclination is {theta_deg:.3} degrees')

# Find alpha 3
# First find mach number after shock
M2 = Mn2/sin(b-theta)
print(f' The Mach number after the shock is {M2:.3}')

mu2 = asin(1/M2)
mu2_deg = mu2*180/pi
print(f'The mach angle after the shock is {mu2_deg:.3}')

a = (mu2 + theta)*180/pi
print(f'Angle relative to horizonal is {a:.3}')

print('--Problem 2b--')

def p_ratio_shock(M,beta,k):
    return 1/(k+1)*( 2*k*M**2*sin(beta)**2-(k-1))

def p_ratio(M,k):
    return (1+1/2*(k-1)*M**2)**(k/(k-1))

print(p_ratio(M1,k))
print(p_ratio_shock(M1,b,k))

p_stag_ratio = p_ratio_shock(M1,b,k)*p_ratio(M2,k)/p_ratio(M1,k)
print(f'Stagnation pressure ratio is {p_stag_ratio:.3}')

print('--Problem 3--')
R = 287
k = 1.4

T = 250
p = 50e3
theta = 10*pi/180
beta = 39*pi/180

a = sqrt(k*R*T)
print(f'The sound velocity is {a:.3}')

def eq(M,k,theta,beta):
    return rho_ratio_mach(M,k,beta)-rho_ratio_theta(beta,theta)

M1 = root_scalar(eq,args=(k,theta,beta), bracket=[0.1,3]).root
print(f'The design Mach number is {M1:.3}')

V1 = M1*a
print(f'The design velocity is {V1:.3}')

print(f'-3b-')
A = 0.5*sin(beta)
rho = p/(R*T)
Q = rho*V1*A
print(f'The mass flow rate is {Q:.5}')

print(f'-3c-')
M1n = M1*sin(beta)
print(f'The normal Mach number before shock is {M1n:.3}')

M2n = sqrt( ((k-1)*M1n**2+2) / (2*k*M1n**2 - (k-1) ) )
print(f'The normal Mach number afte the shock is {M2n:.3}')

M2 = M2n/sin(beta-theta)
print(f'The Mach number after the shock is {M2:.3}')

print('--Problem 4--')

def theta(M1,b,k): 
    return np.arctan (2/np.tan(b)*(M1**2*np.sin(b)**2-1) / (M1**2*(k+np.cos(2*b)) + 2))

b = np.linspace(0,pi/2,100)
Ms = (1.2,1.6, 2, 2.4, 2.8, 3.2, 5, 10, 20)
for M in Ms:
    b_deg = b*180/pi
    t = theta(M,b,k)
    t_deg = t*180/pi
    plt.plot(b_deg,t_deg, label=str(M))
    print(M)
    plt.grid(linestyle='-')
    plt.legend()
plt.ylim(0,60)
plt.xlim(0,100)
plt.show()
