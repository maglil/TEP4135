import sympy as sym

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
