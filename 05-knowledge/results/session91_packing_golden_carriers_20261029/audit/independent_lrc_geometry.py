#!/usr/bin/env python3
"""Independent rational sawtooth arrangement and reflected component images."""
from fractions import Fraction as F
from itertools import combinations
import json
N=0
def ck(x,m):
 global N
 N+=1
 if not x:raise RuntimeError(m)
def dist(x):
 x=x-x.numerator//x.denominator
 return min(x,1-x)
def value(V,t):return min(dist(v*t) for v in V)
def max_cells(V):
 cuts=sorted({F(k,2*v) for v in V for k in range(2*v+1)});best=F(0);clock=F(0)
 for l,r in zip(cuts,cuts[1:]):
  m=(l+r)/2;lines=[]
  for v in V:
   k=(v*m).numerator//(v*m).denominator
   lines.append((v,-k) if v*m-k<F(1,2) else(-v,k+1))
  cand={l,r}
  for (a,b),(c,d) in combinations(lines,2):
   if a!=c:
    t=F(d-b,a-c)
    if l<=t<=r:cand.add(t)
  for t in cand:
   z=value(V,t)
   if z>best:best,clock=z,t
 return best,clock

def good(V,delta):
 cur=[(F(0),F(1))]
 for v in V:
  allowed=[((k+delta)/v,(k+1-delta)/v) for k in range(v) if delta<=F(1,2)]
  out=[]
  for a,b in cur:
   for c,d in allowed:
    x,y=max(a,c),min(b,d)
    if x<=y:
     if out and x<=out[-1][1]:out[-1]=(out[-1][0],max(y,out[-1][1]))
     else:out.append((x,y))
  cur=out
 return cur

def main():
 rows=0
 for k in range(1,5):
  for V in combinations(range(1,9),k):
   lam,t=max_cells(V);M=max(V)
   ck(t.denominator<=2*M and lam.denominator<=2*M,'rational clock bound')
   ck(any((t*(a+b)).denominator==1 for a in V for b in V),'pair-sum time')
   eps=F(1,12*M*M);poly=lambda x:x*x-3*x+1
   ck(poly(lam)!=0,'irrational target')
   ck(poly(lam+eps)>0 if poly(lam)>0 else poly(lam-eps)<0,'exact golden gap')
   levels=sorted({F(1,10),F(1,4),F(1,3),F(2,5),F(1,2),lam,lam/2,lam+F(1,100)})
   sets={a:good(V,a) for a in levels}
   for a in levels:
    S=sets[a]
    ck([(1-y,1-x) for x,y in S][::-1]==S,'reflected intervals')
    ck((bool(S))==(a<=lam),'arrangement versus optimum')
    for b in levels:
     if b>a:continue
     image=[(x,y) for x,y in sets[b] if any(x<=u<=v<=y for u,v in S)]
     trace=sum(x+y==1 for x,y in image)
     ck(trace==(1 if all(v%2 for v in V) and a<=F(1,2) else 0),'image trace')
   rows+=1
 # Retain exact mixed empty/point/interval and AP13 controls.
 mixed=[good((1,2),a) for a in (F(2,5),F(1,3),F(1,4))]
 ck(not mixed[0] and all(a==b for a,b in mixed[1]) and all(a<b for a,b in mixed[2]),'mixed dimension hostile')
 ap=good(tuple(range(1,14)),F(1,14))
 ck(ap==[(F(k,14),F(k,14)) for k in(1,3,5,9,11,13)],'AP13 exact six points')
 print(json.dumps(dict(status='PASS',checks=N,velocity_sets=rows,ap13=[str(a) for a,b in ap]),indent=2))
if __name__=='__main__':main()
