from scipy.optimize import root_scalar

def T_ratio(M,k):
    return 1+(k-1)/2*M**2

def p_ratio(M,k):
    return T_ratio(M,k)**(k/(k-1))

def rho_ratio(M,k):
    return T_ratio(M,k)**(1/(k-1))

def a_ratio(M,k):
    return T_ratio(M,k)**(1/2)

print('--Example 9.8--')

p0 = 120e3
T0 = 400
k = 1.4
R = 287
A = 6e-4 

rho0 = p0/(R*T0)
a0 = (k*R*T0)**(1/2)

p_star = p0/p_ratio(1,k)
print(f'Critical pressure is {p_star:.2}')

print('-a-')
pb = 90e3

def eq(M,k):
    return p_ratio(M,k)-p0/pb

M = root_scalar(eq,args=(k), bracket = [0,1]).root
print(f'Exit mach number is {M:.2}')

rho = rho0/rho_ratio(M,k)
a = a0/a_ratio(M,k)
V = M*a

mdot = rho*V*A
print(f'The mass flow at pb = 90 kPa is {mdot:.2}')


print('-b-')

# Using formula for mdot_max

def mdot_max(rho0, T0, A_star, R, k):
    return rho0*(2/(k+1))**(1/(k-1))*A_star*(2*k*R*T0/(k+1))**(1/2)

rho0 = p0/(R*T0)

mdot_m = mdot_max(rho0,T0,A,R,k)
print(f'Mass flow at 45 kPa (chocked) is {mdot_m:.2}')



