#!/usr/bin/env python3
"""Finite exact checks; no assertions; stdlib only."""
from math import prod, isqrt, comb
import json
from pathlib import Path
B=100000
spf=list(range(B+1))
for p in range(2,isqrt(B)+1):
 if spf[p]==p:
  for n in range(p*p,B+1,p):
   if spf[n]==n: spf[n]=p

def factors(n):
 out=[]
 while n>1:
  p=spf[n]; a=0
  while n%p==0:n//=p;a+=1
  out.append((p,a))
 return out
polys=[lambda a:12*a*a+8*a+1, lambda a:2*a*a+4*a+1, lambda a:4*a+1,lambda a:2*a+1]
weights=[1,-6,-3,14]
N=[0]*(B+1); P=[0]*(B+1)
for n in range(1,B+1):
 f=factors(n); ex=[a for p,a in f if p%4==1]
 num=sum(w*prod(q(a) for a in ex) for w,q in zip(weights,polys))-6
 if num%48: raise RuntimeError(('N integrality',n))
 N[n]=num//48
 if all(p%4==1 for p,a in f) and f:
  num=prod(24*a-4 for a in ex)-6*prod(4*a+2 for a in ex)-3*4**len(ex)+14*2**len(ex)
  if num%48: raise RuntimeError(('P integrality',n))
  P[n]=num//48
 if min(N[n],P[n])<0:raise RuntimeError(('nonnegative',n))
ancestry=[0]*(B+1)
for d in range(1,B+1):
 for n in range(d,B+1,d): ancestry[n]+=P[d]
if ancestry!=N:raise RuntimeError('divisor ancestry failed')
# Independent rational-function expansions: Q(t)/(1-t)^r.
Qs=[[1,18,5],[1,4,-1],[1,3],[1,1]]; rs=[3,3,2,2]
for q,Q,r in zip(polys,Qs,rs):
 for a in range(101):
  got=sum(v*comb(a-j+r-1,r-1) for j,v in enumerate(Q) if j<=a)
  if got!=q(a):raise RuntimeError(('root local',a))
  gotP=sum(v*comb(a-j+r-2,r-2) for j,v in enumerate(Q) if j<=a)
  want=q(a)-(q(a-1) if a else 0)
  if gotP!=want:raise RuntimeError(('primitive local',a))
# Product coefficient reconstruction from independently expanded Euler locals.
def euler_coeff(n,Q,r,primitive):
 v=1
 for p,a in factors(n):
  if p%4!=1:
   if primitive:return 0
   continue
  rr=r-int(primitive)
  v*=sum(c*comb(a-j+rr-1,rr-1) for j,c in enumerate(Q) if j<=a)
 return v
for n in range(1,B+1):
 fn=sum(w*euler_coeff(n,Q,r,False) for w,Q,r in zip(weights,Qs,rs))-6
 gn=sum(w*euler_coeff(n,Q,r,True) for w,Q,r in zip(weights,Qs,rs))-6*int(n==1)
 if fn!=48*N[n] or gn!=48*P[n]:raise RuntimeError(('Euler coefficient',n))
# H split factors must have constant 1 and linear 0.
cs=[10,3,2,1]
for Q,r,c in zip(Qs,rs,cs):
 hpower=2*c+1-r
 if Q[0]!=1 or Q[1]-hpower!=0:raise RuntimeError('H linear cancellation')
samples=[1,5,25,65,130,195,325,845,1105,4225,71825]
report={'status':'PASS','bound':B,'local_exponents_checked':100,'checks':['nonnegative integer closed coefficients','divisor ancestry at all n <= bound','independent rational local expansions','Euler product reconstruction at all n <= bound','H split linear cancellation'],'samples':[{'n':n,'factorization':factors(n),'N':N[n],'P':P[n]} for n in samples], 'first_nonzero':next(n for n in range(1,B+1) if N[n]), 'sums':{'N':sum(N),'P':sum(P)}}
print(json.dumps(report,indent=2))
Path(__file__).with_name('certificate.json').write_text(json.dumps(report,indent=2)+'\n')
