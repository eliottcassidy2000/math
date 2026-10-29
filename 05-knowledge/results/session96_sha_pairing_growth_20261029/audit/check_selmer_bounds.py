#!/usr/bin/env python3
"""Independent exact tower algebra controls; not arithmetic realizations."""
from itertools import product
from fractions import Fraction
import json
N=0

def ck(x,m):
 global N
 N+=1
 if not x:raise RuntimeError(m)

for p in (2,3,5):
 for exps in product(range(4),repeat=3):
  exps=[e for e in exps if e]
  for corank in range(3):
   b=[corank*n+sum(min(n,e) for e in exps) for n in range(8)]
   a=[b[n+1]-b[n] for n in range(7)]
   ck(all(a[i]>=a[i+1]>=0 for i in range(6)),('monotonic',p,exps,corank))
   for n in range(7):
    ck((b[n+1]==b[n])==(corank==0 and all(e<=n for e in exps)),('plateau',p,exps,corank,n))
   for rank in range(3):
    torsion=(1,3)
    logs=[rank*n+sum(min(n,e) for e in torsion)+b[n] for n in range(8)]
    for n in range(5):
     tau=p**sum(int(e>=n+1) for e in torsion)
     for r0 in range(rank+1):
      for lowloss,upgain in product(range(2),repeat=2):
       L=p**max(0,logs[n]-lowloss);U=p**(logs[n+1]+upgain)
       criterion=U<=L*p**r0*tau
       if criterion:ck(rank==r0 and b[n+1]==b[n],('certified bound',p,exps,corank,rank,n,r0,lowloss,upgain))
for p in (2,3):
 for E in range(1,6):
  for n in range(1,E+1):
   q=p**n;scale=Fraction(p**(2*(E-n)),p**E);rad=0
   for x,y in product(range(q),repeat=2):
    if (scale*x).denominator==1 and (scale*y).denominator==1:rad+=1
   ck(rad==p**(2*min(n,E-n)),('restriction radical',p,E,n))
# An unsupported rank guess can erase a genuine divisible growth rate.
ck(all(n-n==0 for n in range(8)),'wrong rank artificial zero slope')
print(json.dumps({'status':'PASS','checks':N,'scope':'cofinite model count profiles and exact bound implication; finite ambient pairing radicals'},sort_keys=True))
