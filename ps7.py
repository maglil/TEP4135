# Problem set 7 - TEP4135

from scipy.optimize import root_scalar

def Mfun(M,k):
    return 1+0.5*(k-1)*M**2

def p_ratio(M,k,ratio=0):
    return Mfun(M,k)**(k/(k-1))-ratio

def critical_area_ratio(M,k,ratio=0):
    return 1/M * ( (1+0.5*(k-1)*M**2) / (0.5*(k+1)) )**(0.5*(k+1)/(k-1))-ratio

# Problem 1
print("Problem 1")

A = 1e-3
A_star = 5e-4
k = 1.5

# Problem 1a
p01 = 20

# Problem 1b
print("Problem 1a")

M = root_scalar(critical_area_ratio, args=(k,A/A_star), method='bisect', bracket=[1,4]).root
print(f"Mach number at exit for chocked, non-shock flow (design), is {M:.3}")
p_design = p01/p_ratio(M,k)
print(f"The design pressure is {p_design:.3}")

M = root_scalar(critical_area_ratio, args=(k,A/A_star), method='bisect', bracket=[.1,1]).root
print(f"Mach number at exit for highest exit pressure whil still critical, is {M:.3}")
p_design = p01/p_ratio(M,k)
print(f"The exit pressure at this Mach number is {p_design:.3}")
