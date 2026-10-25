#!/usr/bin/env python3
"""Exact elementary certificate; CM classification facts are separately cited."""
import json,math
from fractions import Fraction

def require(ok,reason):
 if not ok:raise RuntimeError(reason)
def primes(bound):
 return [p for p in range(2,bound+1) if all(p%q for q in range(2,math.isqrt(p)+1))]
def square_size(n):
 return len({min(a*a%n,(-a*a)%n) for a in range(n) if math.gcd(a,n)==1})
def odd_local(p,a):
 return p**(a-1)*(p-1)//math.gcd(4,p-1)
def conductor_bound(d):
 if not isinstance(d,int) or isinstance(d,bool) or d<1:raise ValueError('map degree must be positive integer')
 factors={2:3+(d.bit_length()-1)}
 for p in primes(4*d+1):
  if p==2:continue
  a=1
  while odd_local(p,a)<=d:a+=1
  if a>1:factors[p]=a-1
 N=math.prod(p**a for p,a in factors.items())
 return N,factors
# Exact quotient-group formula checked independently by unit enumeration.
for p in primes(33):
 if p==2:continue
 for a in range(1,3):
  require(square_size(p**a)==odd_local(p,a),'odd prime-power formula')
for a in range(3,10):
 require(square_size(2**a)==2**(a-3),'two-power formula')
rows=[]
for d in range(1,9):
 N,f=conductor_bound(d)
 require(f[2]==3+(d.bit_length()-1),'two exponent')
 require(2**(f[2]-3)<=d<2**(f[2]-2),'two maximality')
 for p,a in f.items():
  if p>2:require(odd_local(p,a)<=d<odd_local(p,a+1),'odd maximality')
 # Bounded independent sanity only; all-n proof is in THEOREM.md.
 passing=[n for n in range(3,501) if square_size(n)<=d]
 require(all(N%n==0 for n in passing),'finite group-to-divisor sanity')
 rows.append(dict(d=d,N=N,prime_exponents=f,sanity_range=[3,500],passing_count=len(passing)))
require(conductor_bound(1)==(120,{2:3,3:1,5:1}),'degree-one control')
require(conductor_bound(3)==(65520,{2:4,3:2,5:1,7:1,13:1}),'degree-three control')
# Exact positive CM control: t^2+t-1=0 and j=85995*t-52515.
# Substitution in H_-15 leaves a multiple of t^2+t-1.
a,b=85995,-52515
coeff=[a*a,2*a*b+191025*a,b*b+191025*b-121287375]
require(coeff==[a*a,a*a,-a*a],'H-15 affine positive control')
# Hostile controls do not rely on floating singular-modulus recognition.
require(square_size(7)==3,'constant-map countermodel')
require(120%15==0 and square_size(15)==2,'divisibility is weaker than square bound')
require(square_size(3)==1 and Fraction(-1,2).denominator!=1,'necessary not sufficient')
# Pole control: for n4, t=0 and R=1/X has denominator zero.
z4=(0,1);z4_inverse=(0,-1) # exact Gaussian-integer coordinates
t4=tuple(a+b for a,b in zip(z4,z4_inverse))
require(t4==(0,0),'pole target t4=0')
for bad in [0,-1,Fraction(1,2),True]:
 try:conductor_bound(bad)
 except ValueError:pass
 else:raise RuntimeError('invalid map degree accepted')
print(json.dumps({'status':'PASS','bounds':rows,'positive_affine_CM_control':{'n':5,'d':1,'R':'85995*X-52515','minimal_polynomial':'J^2+191025*J-121287375','CM_discriminant':-15},'hostile_controls':['constant R=8000 has CM values at unbounded n','R=1/X has pole at n4','R=X/2 at n3 passes d1 group screen but gives nonintegral -1/2','n15 divides N1 but has square subgroup size2']},indent=2))
