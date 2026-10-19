#!/usr/bin/env python3
"""Exact finite controls for arbitrary-order simple resonance and iteration LTE."""
import json

def vp(x,p):
 k=0
 if x==0:raise ValueError('valuation of zero requested')
 while x%p==0:x//=p;k+=1
 return k
checks=[]
for p in (3,5,7,11,13):
 lam=2;r=next(k for k in range(1,p) if pow(lam,k,p)==1)
 if (lam-1)%p==0 or vp(lam**r-1,p)!=1:raise ValueError('bad control multiplier')
 mod=p**9;tested=0
 for A in range(-2,3):
  for B in range(-2,3):
   def g(x):return (lam*x+A*x*x+B*x*x*x)%mod
   for s in (1,2,3):
    for unit in (1,2,p-1):
     if unit%p==0:continue
     t=p**s*unit
     for n in (1,2,3,p,p+1,p*p,2*p*p):
      x=t
      for _ in range(r*n):x=g(x)
      delta=(x-t)%mod
      expected=s+1+vp(n,p)
      if delta==0 or vp(delta,p)!=expected:raise ValueError((p,A,B,s,unit,n,delta,expected))
      tested+=1
 checks.append(dict(p=p,lambda_value=lam,resonance_order=r,checks=tested))
# A boundary hostile: omit lambda-1-unit, and a nonzero fixed point appears.
# g(t)=(1+p)t-t^2 fixes both 0 and p, even though v_p(lambda-1)=1.
if (4*3-3*3)!=3:raise ValueError('boundary hostile failed')
print(json.dumps(dict(status='PASS',checks=checks,total_checks=sum(r['checks']for r in checks),boundary_hostile='p=3, g(t)=4t-t^2, alpha=0, lambda=4; t=3 is another fixed point, so the lambda-1-unit hypothesis cannot simply be removed.'),indent=2))
