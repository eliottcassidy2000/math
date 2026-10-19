#!/usr/bin/env python3
"""Independent finite search: unrestricted bounded-exponent smooth PAIRS.
No local triple table is used. Enumerate all half-sums u,v with p-exponents
between 0 and (N-1)v_p(b) (or (N-1)(v_2(b)-1)). Derive w and test it.
"""
from itertools import product
from math import gcd
from fractions import Fraction
import json,time,argparse

def factor(b):
 fs=[]
 for p in range(2,b+1):
  if b%p:continue
  e=0
  while b%p==0:b//=p;e+=1
  fs.append((p,e))
  if b==1:break
 return fs

def search(b):
 start=time.time();fs=factor(b);N=2**len(fs)+2;h=b//2
 magnitudes=[1]
 for p,e in fs:
  base=p**(e-(p==2))
  magnitudes=[a*base**k for a in magnitudes for k in range(N)]
 magset=set(magnitudes);cycles=set();triangles=set();divisible_pairs=0
 for A in magnitudes:
  for B in magnitudes:
   H,R=divmod(h*B,A)
   if R:continue
   divisible_pairs+=1
   for s in (1,-1):
    u=s*A
    for t in (1,-1):
     v=t*B;w=v-h+s*t*H
     if abs(w) not in magset:continue
     a=(u+v-w,u+w-v,v+w-u)
     if len(set(a))!=3 or any(gcd(x,b)!=1 for x in a):continue
     C=b*a[1]-a[0]**2
     assert a[1]**2+C==b*a[2]
     triangles.add((*a,C))
     orbit=[];seen=set();x=a[0]
     while x not in seen and len(orbit)<N:
      seen.add(x);orbit.append(x)
      q,r=divmod(x*x+C,b)
      if r or gcd(q,b)!=1:break
      x=q
     else:
      if x==a[0] and len(orbit)>=3:
       cyc=min(tuple(orbit[k:]+orbit[:k]) for k in range(len(orbit)))
       cycles.add((C,cyc))
 return dict(b=b,N=N,magnitude_count=len(magnitudes),pair_count=len(magnitudes)**2,divisible_pairs=divisible_pairs,primitive_triangle_candidates=len(triangles),cycles=[dict(C=C,c=str(Fraction(C,b*b)),numerators=list(a),period=len(a)) for C,a in sorted(cycles)],elapsed_seconds=round(time.time()-start,4))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--max-b',type=int,default=400);p.add_argument('--b',type=int,nargs='*');p.add_argument('--output',default='session7-cycle-review/independent_results.json');args=p.parse_args()
 bs=args.b or [b for b in range(4,args.max_b+1,4) if len(factor(b))==3]
 results=[]
 for b in bs:
  r=search(b);results.append(r)
  print(json.dumps(r),flush=True)
 out=dict(method='Independent full bounded-exponent smooth-pair enumeration without triple signature table',results=results)
 with open(args.output,'w') as f:json.dump(out,f,indent=2)
