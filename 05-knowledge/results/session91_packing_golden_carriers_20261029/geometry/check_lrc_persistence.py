#!/usr/bin/env python3
"""Exact rational LRC contacts, reflection persistence and quadratic gap."""
from fractions import Fraction as F
from itertools import combinations,combinations_with_replacement
import json
CHECKS=0
def need(v,msg='failed'):
 global CHECKS
 CHECKS+=1
 if not v:raise AssertionError(msg)
def dist(x):
 x=x-x.numerator//x.denominator
 return min(x,1-x)
def f(V,t):return min(dist(v*t) for v in V)
def pair_max(V):
 times={F(k,u+v) for u,v in combinations_with_replacement(V,2) for k in range(1,(u+v)//2+1)}
 return max((f(V,t),t) for t in times)
def cell_max(V):
 cuts=sorted({F(k,2*v) for v in V for k in range(2*v+1)})
 candidates=set(cuts)
 for l,r in zip(cuts,cuts[1:]):
  mid=(l+r)/2;lines=[]
  for v in V:
   a=(v*mid).numerator//(v*mid).denominator
   lines.append((v,-a) if v*mid-a<F(1,2) else (-v,a+1))
  for (a,b),(c,d) in combinations(lines,2):
   if a!=c:
    t=F(d-b,a-c)
    if l<=t<=r:candidates.add(t)
 return max(f(V,t) for t in candidates)
def merge(rows):
 out=[]
 for a,b in sorted(rows):
  if out and a<=out[-1][1]:out[-1]=(out[-1][0],max(b,out[-1][1]))
  else:out.append((a,b))
 return out
def good(V,delta):
 rows=[(F(0),F(1))]
 for v in V:
  allowed=[(F(k+delta,v),F(k+1-delta,v)) for k in range(v)]
  rows=merge((max(a,c),min(b,d)) for a,b in rows for c,d in allowed if max(a,c)<=min(b,d))
 return rows
def image_trace(upper,lower):
 hit=[(a,b) for a,b in lower if any(max(a,c)<=min(b,d) for c,d in upper)]
 return sum(a+b==1 for a,b in hit),len(hit)
def above_alpha(r):
 # r > (3-sqrt5)/2, using only rational arithmetic.
 z=3-2*r
 return z<0 or z*z<5

def main():
 levels=[F(1,12),F(1,6),F(1,4),F(1,3),F(2,5),F(1,2)]
 cases=0
 for r in range(1,5):
  for V in combinations(range(1,10),r):
   lam,t=pair_max(V);M=max(V);cases+=1
   need(lam==cell_max(V),'independent sawtooth arrangement')
   need(t.denominator<=2*M and lam.denominator<=2*M,'denominator')
   c=F(1,12*M*M)
   need(above_alpha(lam-c) or not above_alpha(lam+c),'golden gap')
   sets={d:good(V,d) for d in levels}
   for a in levels:
    for b in levels:
     if a<b:continue
     trace,rank=image_trace(sets[a],sets[b]);odd=all(v%2 for v in V)
     need(trace==odd,'cross-threshold reflection trace')
     if not odd:
      need(rank%2==0,'regular representation even image rank')
      need(all(not (l<=F(1,2)<=r) for l,r in sets[a]),'fixed point excluded')
      left=[(l,r) for l,r in sets[a] if r<F(1,2)]
      right=[(l,r) for l,r in sets[a] if l>F(1,2)]
      need(sorted((1-r,1-l) for l,r in left)==right,'natural doubled half')
 controls={str(d):[(str(a),str(b)) for a,b in good((1,2),d)] for d in (F(2,5),F(1,3),F(1,4))}
 need(good((1,2),F(2,5))==[])
 need(good((1,2),F(1,3))==[(F(1,3),F(1,3)),(F(2,3),F(2,3))])
 need(good((1,2),F(1,4))==[(F(1,4),F(3,8)),(F(5,8),F(3,4))])
 ap=good(tuple(range(1,14)),F(1,14))
 need(ap==[(F(k,14),F(k,14)) for k in (1,3,5,9,11,13)],'AP13 boundary')
 print(json.dumps(dict(status='PASS',checks=CHECKS,speed_sets=cases,controls=controls,AP13=[str(a) for a,b in ap]),indent=2))
if __name__=='__main__':main()
