from sympy import *
p1, p2, rho1, rho2, V1, k, y= symbols('p1 p2 rho1 rho2 V1 k y')

eq1 = Eq(p1-p2, (rho1/rho2-1)*rho1*V1**2)
sol = solveset(eq1,rho2)
print(sol)
exp1 = p1/rho1 - p2/rho2
exp2 = V1**2/(2*k/(k-1))*(rho1**2/rho2**2-1)
eq2 = Eq(exp1,exp2)
#eq2 = Eq(p1/rho1 - p2/rho2, V1**2/(2*k/(k-1))*(rho1**2/rho2**2-1))
eq3 = Eq(y,p1/p2)


#eq4 = eq2.subs(p2,sol)
#print(eq4)
