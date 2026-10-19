#!/usr/bin/env python3
"""Independent search using projection of directly enumerated binary words.
The two-edge table is obtained without the author's closed-form triple table.
The third half-sum is calculated from the recurrence, then checked only for
bounded smooth support. No claimed triple-signature formula is imported.
"""
from itertools import product
from math import gcd
from fractions import Fraction
from independent_search import factor
import json,time,argparse

def primitive(w):
 n=len(w)
 return not any(n%d==0 and all(w[i]==w[i%d] for i in range(n)) for d in range(1,n))
def exponent(w,i,j):
 n=len(w)
 if w[i%n]==w[j%n]:return 0
 return next(k for k in range(1,n) if w[(i+k)%n]!=w[(j+k)%n])
def pairs(N):
 result=set()
 for n in range(3,N+1):
  for z in product((0,1),repeat=n-1):
   w=(0,)+z
   if primitive(w):result.add((exponent(w,0,1),exponent(w,0,2)))
 return sorted(result)

def search(b,pairtable):
 start=time.time();fs=factor(b);N=2**len(fs)+2;h=b//2
 magnitudes=[1];local=[]
 for p,e in fs:
  base=p**(e-(p==2))
  powers=[base**k for k in range(N)]
  magnitudes=[a*z for a in magnitudes for z in powers]
  local.append([(powers[i],powers[j]) for i,j in pairtable])
 magset=set(magnitudes);cycles=set();triangles=set();count=0
 for parts in product(*local):
  A=B=1
  for a,bb in parts:A*=a;B*=bb
  H,R=divmod(h*B,A)
  if R:raise ValueError('derived pair not divisible')
  count+=1
  for s in (1,-1):
   u=s*A
   for t in (1,-1):
    v=t*B;w=v-h+s*t*H
    if abs(w) not in magset:continue
    a=(u+v-w,u+w-v,v+w-u)
    if len(set(a))!=3 or any(gcd(x,b)!=1 for x in a):continue
    C=b*a[1]-a[0]**2
    if a[1]**2+C!=b*a[2]:raise ValueError('recurrence mismatch')
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
 return dict(b=b,N=N,local_pair_count=len(pairtable),unsigned_pair_count=count,primitive_triangle_candidates=len(triangles),cycles=[dict(C=C,c=str(Fraction(C,b*b)),numerators=list(a),period=len(a)) for C,a in sorted(cycles)],elapsed_seconds=round(time.time()-start,4))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--max-b',type=int,default=1000);p.add_argument('--b',type=int,nargs='*');p.add_argument('--output',default='session7-cycle-review/projected_results.json');args=p.parse_args()
 bs=args.b or [b for b in range(4,args.max_b+1,4) if len(factor(b))>=3]
 tables={N:pairs(N) for N in sorted({2**len(factor(b))+2 for b in bs})}
 results=[]
 for b in bs:
  r=search(b,tables[2**len(factor(b))+2]);results.append(r)
  print(json.dumps(r),flush=True)
 out=dict(method='Independent exhaustive primitive-word projection and exact two-half-sum enumeration',pair_tables=tables,results=results)
 with open(args.output,'w') as f:json.dump(out,f,indent=2)
