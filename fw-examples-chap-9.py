from scipy.optimize import root_scalar
import numpy as np
import matplotlib.pyplot as plt

gamma = 1.4
R = 287

A1 = 0.05
v1 = 180
T1 = 470
p1 = 500e3

# Example 9.4

# b)
M = v1/(gamma*R*T1)**(1/2)
print("M: " + str(M))

#a)
T0 = (1 + (gamma-1)/2*M**2)*T1
print("T0: " + str(T0))

# c)
p0 = (T0/T1)**(gamma/(gamma-1))*p1
print("p0: " + str(p0))

# d)
m = p1/R/T1*v1*A1
print("m: " + str(m))

f = 1/M*( (1+1/2*(gamma-1)*M**2)/(gamma+1)*2 )**((gamma+1)/(2*(gamma-1)))
Ac = A1/f
print("Ac: " + str(Ac))

A2 = 0.036

def r(M):
    k=1.4
    return 1/M*( (1 + 1/2*(k-1)*M**2)/(1/2*(k+1)) )**((k+1)/(2*(k-1))) - A2/Ac

Ma = np.linspace(0.1,3,100)
plt.plot(Ma,r(Ma))
plt.grid(True)


plt.show()
root = root_scalar(r, bracket=[1.1,2.9],method='bisect')
print("root: " + str(root))

