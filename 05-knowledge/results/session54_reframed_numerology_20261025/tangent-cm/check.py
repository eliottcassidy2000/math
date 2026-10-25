#!/usr/bin/env python3
"""Exact certificate for tangent-Legendre CM classification; no floating arithmetic."""
import json, math
from itertools import product
from fractions import Fraction as Q
from pathlib import Path
import sympy as s
z,t,J=s.symbols('z t J')
def require(ok,msg):
 if not ok: raise RuntimeError(msg)
def squares(n):
 u={a*a%n for a in range(n) if math.gcd(a,n)==1}
 return len(u),len({min(a,(-a)%n) for a in u})
def polys(n):
 # Phi_n(z)=z^(phi(n)/2)*Psi_n(z+z^-1), independently reconstructed below.
 phi=s.Poly(s.cyclotomic_poly(n,z),z)
 psi=s.Poly(s.resultant(phi.as_expr(),z*z-t*z+1,z),t).sqf_part().monic()
 degree=psi.degree()
 reconstructed=s.Poly(s.cancel(z**degree*psi.as_expr().subs(t,z+1/z)),z)
 require(reconstructed==phi, f'cyclotomic reconstruction n={n}')
 # The resultant is a power of the minimal polynomial, as psi is irreducible.
 rel=J*(2-t)**2*(2+t)-64*(10+3*t)**3
 raw=s.Poly(s.resultant(psi.as_expr(),rel,t),J)
 minimal=raw.sqf_part().monic()
 require(minimal.is_irreducible,f'irreducibility n={n}')
 # Independent field identity: substitute j rational function into minimal polynomial.
 ring=s.QQ.poly_ring(t)
 modulus=s.Poly(psi,t,domain=s.QQ)
 denominator=s.Poly((2-t)**2*(2+t),t,domain=s.QQ)
 numerator=s.Poly(64*(10+3*t)**3,t,domain=s.QQ)
 jf=(numerator*s.invert(denominator,modulus)).rem(modulus)
 value=s.Poly(0,t,domain=s.QQ)
 for coefficient in minimal.all_coeffs():
  value=(value*jf+s.Poly(coefficient,t,domain=s.QQ)).rem(modulus)
 require(value.is_zero,f'field substitution n={n}')
 return psi,minimal
expected=[3,4,5,6,7,8,9,10,12,13,14,15,16,18,20,21,24,26,28,30,36,40,42,48,56,60,72,84,120,168]
ns=[]
for a,b,c,d,e in product(range(6),range(3),range(2),range(2),range(2)):
 n=2**a*3**b*5**c*7**d*13**e
 if n>=3 and squares(n)[1]<=3:ns.append(n)
ns.sort();require(ns==expected,'complete group screen')
rows=[];cache={}
for n in ns:
 psi,p=polys(n);cache[n]=(psi,p)
 require(psi.degree()==p.degree(),f'no degree loss n={n}')
 rows.append(dict(n=n,real_degree=psi.degree(),j_degree=p.degree(),square_subgroup=squares(n)[1],integral=all(c.q==1 for c in p.all_coeffs()),polynomial=str(p.as_expr())))
elementary=[r['n'] for r in rows if r['square_subgroup']==1]
require(elementary==[3,4,5,6,8,10,12,24],'elementary screen')
survivors=[r['n'] for r in rows if r['square_subgroup']==1 and r['integral']]
require(survivors==[4,8,12,24],'integrality screen')
# All conjugates are >1728: exact numerator for j(t)-1728 is a square;
# equality t=-14/9 is not a cyclotomic integer.
require(s.expand(64*(10+3*t)**3-1728*(2-t)**2*(2+t))==s.expand(256*(9*t+14)**2),'positive conjugates')
# Thus each root is smaller than its trace. Bound the dominant CM conjugate
# by trace and |j|-exp(pi sqrt(|D|)) <=2079 (cited proposition).
# Certify exp(pi sqrt(B)) > trace+2079 using rational Taylor LOWER bounds,
# pi>157/50, and lower approximations r to sqrt(B).
bounds=[]
for n,B,r in [(8,15,Q(387,100)),(12,19,Q(435,100)),(24,32,Q(565,100))]:
 p=cache[n][1];trace=int(-p.all_coeffs()[1]);v=Q(157,50)*r
 require(r*r<B,'sqrt lower bound')
 exp_lower=sum((v**k/Q(math.factorial(k))) for k in range(100))
 require(exp_lower>trace+2079,'exponential lower bound')
 bounds.append(dict(n=n,trace=trace,discriminant_abs_strict_upper=B,exp_argument_lower=str(v),taylor_terms=100))
# Reduced primitive positive binary forms, canonical boundary convention.
def forms(D):
 out=[]
 for a in range(1,math.isqrt((-D)//3)+1):
  for b in range(-a,a+1):
   if (b*b-D)%(4*a):continue
   c=(b*b-D)//(4*a)
   if a>c or math.gcd(math.gcd(a,abs(b)),c)!=1:continue
   if (abs(b)==a or a==c) and b<0:continue
   out.append([a,b,c])
 return out
classes={str(-d):forms(-d) for d in range(3,32) if (-d)%4 in (0,1)}
require(all(len(v)<2 for k,v in classes.items() if -int(k)<15),'no degree2 below15')
require([int(k) for k,v in classes.items() if -int(k)<19 and len(v)==2]==[-15],'degree2 below19')
require(all(len(v)<4 for v in classes.values()),'no degree4 below32')
# sole possible n12 discriminant -15 excluded by exact published H_-15.
H15=J**2+191025*J-121287375
require(s.gcd(cache[12][1].as_expr(),H15)==1,'H-15 mismatch')
box=[]
for m in range(5,25):
 n=m//math.gcd(m,2)
 if n not in cache:cache[n]=polys(n)
 p=cache[n][1]
 box.append(dict(m=m,n=n,degree=p.degree(),denominator=int(s.ilcm(*[c.q for c in p.all_coeffs()])),polynomial=str(p.as_expr())))
require(cache[4][1].as_expr()==J-8000,'positive CM control')
require(cache[3][1].as_expr()==J-Q(21952,9),'m6 hostile')
require(cache[6][1].as_expr()==J-Q(140608,3),'m12 hostile')
result=dict(status='EXACT CERTIFICATE PASSED; classification relies on explicitly cited CM theorems',candidate_conductors=rows,elementary_conductors=elementary,integral_elementary_conductors=survivors,discriminant_bounds=bounds,reduced_forms_abs_discriminant_below32=classes,requested_box=box,CM_integer_m_greater_than4=[8])
print(json.dumps(result,indent=2))
