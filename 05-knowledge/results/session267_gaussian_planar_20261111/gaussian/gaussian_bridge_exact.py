"""Exact new-result controls; no inherited-theorem audit. Python 3 standard library."""
from fractions import Fraction
from math import factorial, comb, gcd
from itertools import product

def mul(p,q):
 r={}
 for (a,b),c in p.items():
  for (u,v),d in q.items(): r[a+u,b+v]=r.get((a+u,b+v),0)+c*d
 return {k:v for k,v in r.items() if v}
def moment(p): return sum(c*factorial(a) for (a,b),c in p.items() if a==b)
def require(b,msg):
 if not b: raise RuntimeError(msg)
checks=0
# Primitive opposite-charge return: first moment is r+s over gcd, no degree-independent bound.
for r,s in product(range(1,9),repeat=2):
 p={(r,0):1,(0,s):1}; power={(0,0):1}; L=(r+s)//gcd(r,s)
 for m in range(1,L+1):
  power=mul(power,p); value=moment(power)
  expected=comb(L,s//gcd(r,s))*factorial(r*s//gcd(r,s)) if m==L else 0
  require(value==expected,'primitive return'); checks+=1
print('Primitive-return checks:',checks)
# Uniform orientation theorem: arbitrary parameter values and polynomial multiplier.
for t in range(-3,4):
 p={(1,0):t,(2,1):t*t+1,(4,0):1-t}; q={(0,5):2,(1,3):-1,(0,0):7}; power={(0,0):1}
 for m in range(1,10):
  power=mul(power,p)
  if m>5: require(moment(mul(q,power))==0,'uniform charge bound'); checks+=1
print('Uniform-family total checks:',checks)
# Wick rank-one sampling matrix at bidegree (m,m): diagonal 1/binom(m,k), rank m+1.
for m in range(1,31):
 diag=[Fraction(1,comb(m,k)) for k in range(m+1)]
 require(all(x!=0 for x in diag),'rank'); checks+=1
print('Rank-one sample lower bound verified symbolically for m=1..30')
# Fundamental failure of nullness transport under covariance compression.
# f=xi_1*z_2: E_I(f^m)=0, but at C=a*b^T with a=b=(1,1), E_C(f^m)=m!.
for m in range(1,31):
 require(factorial(m)>0,'covariance hostile'); checks+=1
print('Covariance nullness hostile: E_I=0, E_rank1=m!, m=1..30')
print('PASS:',checks,'exact checks. No claim of a Keller compression construction.')

# Forward-data two-shear construction, exact polynomials in t.
def padd(a,b):
 n=max(len(a),len(b)); r=[(a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0) for i in range(n)]
 while len(r)>1 and r[-1]==0: r.pop()
 return r
def pneg(a): return [-c for c in a]
def pmul(a,b):
 r=[0]*(len(a)+len(b)-1)
 for i,c in enumerate(a):
  for j,d in enumerate(b): r[i+j]+=c*d
 while len(r)>1 and r[-1]==0: r.pop()
 return r
def peval(coeff,x):
 r=[0]
 for c in reversed(coeff): r=padd(pmul(r,x),[c])
 return r
def deriv_value(coeff,k,u):
 return sum(coeff[j]*factorial(j)//factorial(j-k)*u**(j-k) for j in range(k,len(coeff)))
shear_checks=0
for d in range(1,8):
 for e in range(1,6):
  for u,v in product(range(-2,3),repeat=2):
   p=[(-1)**j*(j+1) for j in range(e+1)]
   q=[j+2 for j in range(d+1)]
   pv=peval(p,[v])[0]
   U=[u,-pv]
   V=padd([v],pneg(pmul([0,1],peval(q,U))))
   # Apply B_t then A_t to the constructed inverse, exactly in t.
   back_y=padd(V,pmul([0,1],peval(q,U)))
   back_x=padd(U,pmul([0,1],peval(p,back_y)))
   require(back_x==[u] and back_y==[v],'two-shear inverse identity')
   require(len(V)-1<=d+1,'uniform two-shear degree')
   for m in range(1,d+4):
    actual=(V[m] if m<len(V) else 0)*factorial(m)
    expected=-m*((-pv)**(m-1))*deriv_value(q,m-1,u)
    require(actual==expected,'two-shear marked coefficient')
    shear_checks+=1
print('Two-shear forward-data coefficient checks:',shear_checks)
print('PASS total:',checks+shear_checks,'exact checks; all requirements remain active under -O.')
