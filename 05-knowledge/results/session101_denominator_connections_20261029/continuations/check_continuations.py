#!/usr/bin/env python3
"""Exact stdlib tests; gates survive python -O. Output deterministic JSON."""
from fractions import Fraction
from itertools import product
from collections import Counter
from math import gcd, lcm
from pathlib import Path
import json
checks=0
def gate(x,msg):
 global checks
 checks+=1
 if not x: raise RuntimeError(msg)
def state(w,b=3):
 A=B=1;C=0
 for v in w:A,B,C=b*A,2**v*B,b*C+B
 return A,B,C
def decode(k,H,C,b=3):
 out=[]
 while k>1:
  rest=C-b**(k-1)
  if rest<=0:raise ValueError("not a positive word carry")
  a=vp(rest,2)
  if a<1:raise ValueError("not a positive valuation letter")
  out.append(a);H-=a;C=rest//2**a;k-=1
 if C!=1 or H<1:raise ValueError("invalid terminal word state")
 return tuple(out+[H])
def closed(w,b=3):
 A,B,C=state(w,b);return Fraction(C,B-A)
def independent(w,b=3):
 # compose affine maps over Fraction, unlike integer recurrence
 slope=Fraction(1);intercept=Fraction(0)
 for v in w:slope=b*slope/2**v;intercept=(b*intercept+1)/2**v
 return intercept/(1-slope)
def vp(n,p):
 if not n:return 10**9
 z=0
 while n%p==0:z+=1;n//=p
 return z
def factors(n):
 ans=[];p=2
 while p*p<=n:
  if n%p==0:
   e=0
   while n%p==0:n//=p;e+=1
   ans.append((p,e))
  p+=1
 if n>1:ans.append((n,1))
 return ans
def row(w):
 A,B,C=state(w);x=closed(w);y=closed(w+(1,))
 return dict(word=w,C=C,source=str(x),source_denominator=x.denominator,append_one=str(y),append_denominator=y.denominator)
words=[w for k in range(1,5) for w in product(range(1,5),repeat=k)]
suffixes=[()] + [w for k in range(1,3) for w in product(range(1,4),repeat=k)]
comparisons=0
for b in (3,5,7):
 for u in words:
  A,B,C=state(u,b); x=closed(u,b)
  gate(x==independent(u,b),'independent rational closure')
  gate(decode(len(u),sum(u),C,b)==u,'exact carry reconstruction')
  gate(gcd(x.denominator,2*b)==1,'unit denominator')
  t=x
  for i,v in enumerate(u):
   gate(vp(b*t.numerator+t.denominator,2)==v,'actual valuation')
   t=(b*t+1)/2**v
   gate(t.denominator==x.denominator,'orbit denominator invariance')
   gate(t==closed(u[i+1:]+u[:i+1],b),'rotation')
  for v in suffixes:
   Av,Bv,Cv=state(v,b); N=Av*C+B*Cv; Delta=B*Bv-A*Av
   gate(state(u+v,b)==(A*Av,B*Bv,N),'composition')
   gate(closed(u+v,b)==Fraction(N,Delta),'suffix closure')
   for p,e in factors(abs(Delta)):
    root=(-B*Cv*pow(Av,-1,p**e))%(p**e)
    val=max(0,e-vp(C-root,p))
    gate(vp(closed(u+v,b).denominator,p)==val,'p-adic cancellation ball')
   comparisons+=1
fibre=[(a,b,7-a-b) for a in range(1,6) for b in range(1,7-a)]
rows=[row(w) for w in fibre]
gate(all(r['source_denominator']==101 for r in rows),'101 fibre')
spectrum=dict(sorted(Counter(r['append_denominator'] for r in rows).items()))
gate(spectrum=={7:1,25:1,35:5,175:8},'full spectrum')
anchor=[row((1,1,5)),row((1,3,3))]
gate([r['append_denominator'] for r in anchor]==[35,175],'anchor')
rotations=[row((1,3,3)),row((3,3,1)),row((3,1,3))]
gate([r['append_denominator'] for r in rotations]==[175,175,7],'necklace hostile')
# Exact modulus for full residue responses for a specified finite menu.
u=(1,1,5); A,B,C=state(u)
menu=[(1,),(2,),(1,2)]
mods=[abs(B*state(v)[1]-A*state(v)[0]) for v in menu]
L=lcm(*mods)
for shift in range(-100,101):
 equal=all(state(v)[0]*shift%m==0 for v,m in zip(menu,mods))
 gate(equal==(shift%L==0),'exact full-residue menu modulus')
# Two input numerators agreeing mod p^t have equal cancellation depths below t.
for C1 in range(1,80,2):
 for C2 in range(1,80,2):
  for p,e in [(5,2),(7,1)]:
   N1,N2=3*C1+128,3*C2+128;t=vp(C1-C2,p)
   r1,r2=min(e,vp(N1,p)),min(e,vp(N2,p))
   gate(min(r1,t)==min(r2,t),'cancellation precision stability')
result=dict(status='FINITE-EXACT',universe=dict(multipliers=[3,5,7],prefix_lengths=[1,4],prefix_letters=[1,4],suffix_lengths=[0,2],suffix_letters=[1,3]),prefix_count=len(words),suffix_count=len(suffixes),composition_comparisons=comparisons,checks=checks,anchor=anchor,necklace_hostile=rotations,clock_3_7_fibre=rows,append_one_spectrum=spectrum,finite_menu=dict(suffixes=menu,moduli=mods,lcm=L),limitations=['No claim of minimal denominator-only state','No integer Collatz convergence or cycle exclusion theorem'])
print(json.dumps(result,indent=2,sort_keys=True))
