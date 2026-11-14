"""Tiny exact arithmetic receipt (requires SymPy)."""
from collections import Counter
import sympy as s

t, X = s.symbols('t X')
f = 30*t**3-5*t**2-10*t-1
assert all(f.subs(t, j) % 11 for j in range(11))
assert s.discriminant(f, t) == 70700
D = (10+s.I)/(10-s.I)
c = X**3-2*X**2+2*D*X-D
assert s.cancel((10-s.I)*(1-s.I*t)**3*c.subs(X,(1+s.I*t)/(1-s.I*t))+2*s.I*f) == 0
R = {1,3,4,5,9}
N = set(range(1,11))-R
assert Counter((x*y)%11 for x in R for y in R) == Counter({j:5 for j in R})
assert Counter((x+y)%11 for x in R for y in N) == Counter({0:5,**{j:2 for j in range(1,11)}})
phi = s.cyclotomic_poly(11, X)
assert s.rem(sum(X**j for j in R)*sum(X**j for j in N)-3,phi,X)==0
print('irreducible real cubic, exact trace 2, and cyclotomic gain of norm 3: PASS')
