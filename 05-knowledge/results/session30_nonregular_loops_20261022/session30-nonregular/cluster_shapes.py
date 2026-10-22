from fractions import Fraction
from itertools import permutations
from collections import defaultdict,Counter
from math import factorial

def partitions(n):
 a=[0]*n
 def rec(i,m):
  if i==n: yield tuple(a);return
  for t in range(m+2):
   a[i]=t;yield from rec(i+1,max(m,t))
 if n: yield from rec(1,0)

def comp(n):
 if n==0:yield ();return
 for x in range(2,n+1,2):
  for rest in comp(n-x):yield (x,)+rest

def canon(edges,v):
 best=None; signs=set()
 for p in permutations(range(v)):
  es=[];sgn=1
  for a,b in edges:
   a,b=p[a],p[b]
   if a>b:a,b=b,a;sgn=-sgn
   es.append((a,b))
  key=tuple(sorted(es))
  if best is None or key<best:best=key;signs={sgn}
  elif best==key:signs.add(sgn)
 return best, next(iter(signs)) if len(signs)==1 else 0

def cluster(L):
 ans=defaultdict(Fraction)
 for ls in comp(L):
  ed=[];o=0
  for l in ls: ed += [(o+i,o+i+1) for i in range(l)];o+=l+1
  for p in partitions(o):
   v=max(p)+1
   if v>L+1:continue
   es=[(p[a],p[b]) for a,b in ed]
   if any(a==b for a,b in es):continue
   seen={0}
   while True:
    z=seen|{b for a,b in es if a in seen}|{a for a,b in es if b in seen}
    if z==seen:break
    seen=z
   if len(seen)<v:continue
   mu=1
   for s in Counter(p).values():mu*=(-1)**(s-1)*factorial(s-1)
   key,sign=canon(es,v)
   if sign:ans[key]+=Fraction(mu*sign,factorial(len(ls)))
 return {k:v for k,v in ans.items() if v}
if __name__=='__main__':
 for L in [2,4]:
  for es,c in cluster(L).items():
   d=Counter(a for e in es for a in e);v=len(d)
   print(L,'e',L-v,'leaves',sum(x==1 for x in d.values()),'coefficient',c,'edges',es)
