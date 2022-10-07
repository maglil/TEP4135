# Problem 3a

v = 150
gamma = 1.4 # air
R = 287 # air
T = 277
p = 70

M = v / (gamma*R*T)**0.5
print(M)

# critical temperature
Tc = 2/(gamma+1)*(1+ (gamma-1)/2*M**2)*T
print(Tc)

#stagnation pressure
p0 = ( 1 + (gamma-1)/2*M**2)**(gamma/(gamma-1))*p
print(p0)

# stagnation temperature
T0 = (1+ (gamma-1)/2*M**2)*T
print(T0)

# relativ density rho*/rho
drho = (2/(1+gamma))*(1+(gamma-1)/2*M**2)**(1/(gamma-1))
print("drho: " + str(drho))

# sonic pressure
pc = (2/(gamma+1))**(gamma/(gamma-1))*p0
print(pc)

# critical velocity
c = (gamma*R*Tc)**0.5
print(c)

# relative area decrease Ac/A = rho*v/rhoc*v_c
dA = v/c/drho
print("dA: " + str(dA))

# area decrease
decA = 1-dA
print("dec: " + str(decA))

dp = (2/(gamma+1))**(gamma/(gamma-1))
print("dp: " + str(dp))

# Problem 4
print("Problem 4")

A = 0.05
v = 180
p = 500e3
T = 470

c = (gamma*R*T)**0.5

# Mach number
M = v/c
print("M: " + str(M))

# Stagnation temperature
T0 = T*(1+(gamma-1)/2*M**2)
print("T0: " + str(T0))

# Critical temperature
Tc = T0*(2/(gamma+1))

# Stagnation pressure
p0 = p*(T0/T)**(gamma/(gamma-1))
print("p0: " + str(p0))

# density
rho = p/(R*T)

# stagnation density
rho0 = (T0/T)**(1/(gamma-1))*rho

# crtical density
rhoc = rho0*(2/(gamma+1))**(1/(gamma-1))

# Speed of sound
c = (gamma*R*Tc)**0.5

Ac = rho*v*A/(rhoc*c)
print("Ac: " + str(Ac))

# Mass flux
m = rho*v*A
print("mdot: " + str(m))


