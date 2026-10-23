from itertools import combinations, product
from fractions import Fraction
from math import gcd
import json
import sympy as s

def need(x,m):
 if not x:raise RuntimeError(m)
def fac(b):
 return tuple(next(k for k in range(12) if b%p**(k+1)) for p in (2,3,5))
def h(e,f):return min(a//b for a,b in zip(e,f) if b)
bases=(2,3,4,5,6,10,12,15,18,20)
profiles={b:fac(b) for b in bases}
box=list(product(range(7),repeat=3));count=0
for b,c in combinations(bases,2):
 f,g=profiles[b],profiles[c]
 attained={(h(e,f),h(e,g)) for e in box}
 for k,l in product(range(4),repeat=2):
  lower=tuple(max(k*a,l*b) for a,b in zip(f,g))
  need(((h(lower,f),h(lower,g))==(k,l))==((k,l) in attained),'least criterion')
  count+=1
z=s.symbols('z');polychecks=0
for a,b in product(range(1,8),repeat=2):
 n=a+b;d=gcd(a,b);P=s.Rational(n**n,a**a*b**b)*z**a*(1-z)**b
 rhs=s.Rational(n**n,a**a*b**b)*z**(a-1)*(1-z)**(b-1)*(a-n*z)
 need(s.expand(s.diff(P,z)-rhs)==0,'derivative')
 aa=a//d;bb=b//d;nn=aa+bb
 Q=s.Rational(nn**nn,aa**aa*bb**bb)*z**aa*(1-z)**bb
 need(s.expand(P-Q**d)==0,'decomposition')
 need(P.subs(z,s.Rational(a,n))==1,'critical value')
 polychecks+=1
print(json.dumps({'status':'PASS','independent_three_prime_constraint_pairs':count,'polynomial_derivative_decomposition_checks':polychecks,'sphere_chi':2,'orbifold_10_2_7':str(Fraction(1,10)+Fraction(1,2)+Fraction(1,7)-1)},indent=2,sort_keys=True))
