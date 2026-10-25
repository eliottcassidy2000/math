#!/usr/bin/env python3
"""Separate addendum checks. No writes to frozen base artifacts."""
from pathlib import Path
from itertools import product
from math import prod, isqrt, gcd
import json

def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def mul(z,w):return z[0]*w[0]-z[1]*w[1],z[0]*w[1]+z[1]*w[0]
def power(z,n):
 r=(1,0)
 for _ in range(n):r=mul(r,z)
 return r
def gaussian_prime(p):
 for b in range(1,isqrt(p)+1):
  a=isqrt(p-b*b)
  if a*a+b*b==p:return a,b
 raise RuntimeError('split prime missing')
def val(n,p):
 a=0
 while n%p==0:n//=p;a+=1
 return a
rows=[]
for ell,ps in ((3,(5,)),(5,(5,)),(5,(5,13)),(7,(5,)),(7,(5,13))):
 a=(ell-1)//2;k=len(ps);e=prod(p**a for p in ps);gapset=set();fibers={};indexed={}
 for j in product(range(-a,a+1),repeat=k):
  if not any(j) or next(v for v in j if v)<0:continue
  z=(1,0)
  for p,v in zip(ps,j):
   w=gaussian_prime(p);w=w if v>=0 else (w[0],-w[1]);z=mul(z,power(w,2*abs(v)));z=(z[0]*p**(a-abs(v)),z[1]*p**(a-abs(v)))
  d=abs(2*z[0]*z[1]);need(d>0 and d not in gapset,'gap uniqueness');gapset.add(d);indexed[j]=d
  need(isqrt(e*e-d)**2==e*e-d and isqrt(e*e+d)**2==e*e+d,'endpoint squares')
  need(gcd(e*e,d)==prod(p**(2*(a-abs(v))) for p,v in zip(ps,j)),'full exponent gcd')
  v=tuple(x%ell for x in j);iv=pow(next(x for x in v if x),-1,ell);pt=tuple(x*iv%ell for x in v)
  fibers.setdefault(pt,[]).append((j,d))
 need(len(gapset)==(ell**k-1)//2,'signed quotient count')
 need(len(fibers)==(ell**k-1)//(ell-1),'projective quotient count')
 need(all(len(v)==a for v in fibers.values()),'uniform projective fiber')
 if e<=5000:
  brute={e*e-u*u for u in range(1,e) if isqrt(2*e*e-u*u)**2==2*e*e-u*u}
  need(brute==gapset,'independent endpoint census')
 if ell==5 and ps==(5,):
  need(indexed=={(1,):600,(2,):336},'minimal valuation hostile')
  need(val(600,5)==2 and val(336,5)==0,'distinct valuation in same projective fiber')
 rows.append({'ell':ell,'a':a,'actual_split_primes':ps,'center':e,'gap_count':len(gapset),'projective_points':len(fibers),'fiber_size':a})
result={'status':'PASS','checks':rows,'minimal_hostile':{'ell':5,'center':25,'indices':[1,2],'gaps':[600,336],'v5':[2,0],'projective_point':'the unique point of P^0(F_5)'}}
Path(__file__).with_name('uniform-verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
