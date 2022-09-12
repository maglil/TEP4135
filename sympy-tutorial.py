from sympy import *

t = symbols('t')

y = Function('y')

f = dsolve(Eq(y(t).diff(t,t)-y(t),exp(t)),y(t))

ev = Matrix([[1, 2],[3, 4]]).eigenvals()

ev
