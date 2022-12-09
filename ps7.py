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
Ae = 1e-3
A_star = 5e-4
A2 = 5e-4
Ae = 1e-3
k = 1.4

# Problem 1a
print("-1a-")

a0 = sqrt(k*R*T0)
print(f"Stagnation speed of sound is {a0:.4}")
a1 = a0/a_ratio(1,k)
print(f"Speed of sound first nozzle (at M=1) is {a1:.4}")
rho0 = p0/(R*T0)
print(f'Stagnation density is {rho0:.4}')
rho1 = rho0/rho_ratio(1,k)
print(f'Density in first nozzle is {rho1:.4}')
mdot = rho1*a1*A1
mdot_lf = (2/k+1)**((k+1)/(2*(k-1)))*A1*p0*sqrt(k/(R*T0))
print(f'Mass flow is {mdot:.4}')
print(f'Mass flow from LF is {mdot:.4}')


rho2 = mdot/(a1*A2)

print(f'Density in second throat is {rho2:.4}')
rho02 = rho2*rho_ratio(1,k)
print(f'The stagnation density is then {rho02:.4}')

p1 = rho02*R*T0
print(f'Pressure in tank is {p1:.4}')

p01 = p1

print('-1b-')
p2 = p01*0.9395
print(p2)
p_star = p01/p_ratio(1,k)
print(p_star)

print('-1c-')
A_shock = 6.5e-4

print('Table method')

r = A_shock/A_star
print(f'Critical area ratio at shock is {r:.3}')

# Mach number before shock
M_shock_1 = 1.6 + (r-1.2502)/(1.3376-1.2502)*(1.7-1.6)
print(f'The Mach number before shock is {M_shock_1}:2')

# Pressure p1 before shock
pr =  0.22 # isentropic, eyeball
p1 = pr*p01

# Across shock from table
M2 = 0.65 # eyeball
# Pressure ratio
prs = 3.0 # shock, eyeball
# Pressure after shock
p2 = p1*prs
print(p2)
# Change in A_star acrooss shock
A_shock_ratio = 1.13
# New A*
A2s = A_star*A_shock_ratio
print(A2s)


# Isentropic Area ratio exit
Ae_ratio = Ae/A2s
print(f'Exit A star ratio is {Ae_ratio:.3}')
pe_ratio = 0.91
pe = p2*pe_ratio
print(f'The exit pressur is {pe:.4}')

print('Problem 2')
T0 = 273.15+30
A1_star = 1e-4
A_shock = 1.25e-4
A_exit = 2e-4

print('-2b-')
area_ratio_shock = A_shock/A1_star
print(f'Area ratio at shock is {area_ratio_shock:.4}')
M1 = 1.6 # eyeballing, from table
print(f'Mach number before shock is {M1:.3}')
a_star_ratio = 1.1171
A2_star = A1_star*a_star_ratio
print(A2_star)
area_ratio_exit = A_exit/A2_star
print(f'Area ratio at exit is {area_ratio_exit:.3}')
Me = 0.35 #eyeballing from table
print(f'Exit mach number is {Me:.4}')

        

print("Problem 1a")

M = root_scalar(critical_area_ratio, args=(k,A/A_star), method='bisect', bracket=[1,4]).root
print(f"Mach number at exit for chocked, non-shock flow (design), is {M:.3}")
p_design = p01/p_ratio(M,k)
print(f"The design pressure is {p_design:.3}")

M = root_scalar(critical_area_ratio, args=(k,A/A_star), method='bisect', bracket=[.1,1]).root
print(f"Mach number at exit for highest exit pressure whil still critical, is {M:.3}")
p_design = p01/p_ratio(M,k)
print(f"The exit pressure at this Mach number is {p_design:.3}")



