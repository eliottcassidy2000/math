#!/usr/bin/env python3
"""Exact finite controls. No dependencies; explicit gates survive python -O."""
from itertools import product, combinations
from math import gcd
import json
checks=0
def check(v,msg):
 global checks
 checks+=1
 if not v: raise RuntimeError(msg)
def phi(n):return sum(gcd(k,n)==1 for k in range(n))
def compose(p,q):return tuple(p[q[i]] for i in range(len(p)))
def perm(r,m,t,z):return tuple((i+t)%r+r*((j+z[i])%m) for j in range(m) for i in range(r))
def order(p):
 q=tuple(range(len(p))); k=0
 while True:
  q=compose(p,q);k+=1
  if q==tuple(range(len(p))):return k

def subgroup(p):
 q=tuple(range(len(p)));out=[]
 while q not in out:
  out.append(q);q=compose(p,q)
 return frozenset(out)
def edge_orbits(parents,n):
 edges=list(combinations(range(n),2));idx={e:i for i,e in enumerate(edges)};uf=list(range(len(edges)))
 def root(x):
  while uf[x]!=x:uf[x]=uf[uf[x]];x=uf[x]
  return x
 for H in parents:
  p=next(p for p in H if order(p)==n)
  for k,(x,y) in enumerate(edges):
   a,b=root(k),root(idx[tuple(sorted((p[x],p[y])))]);uf[a]=b
 return {e:root(k) for k,e in enumerate(edges)}
def alpha(t):
 s=t//2;w=[]
 for k in range(s):w += [k,t-1-k]
 w += [s]
 return [2*x%t for x in w]
rows=[]
for r,m in [(3,3),(3,5),(3,7),(5,3),(3,9),(5,5)]:
 n=r*m;g=gcd(r,m);parents={};allw=[]
 for t in range(r):
  for z in product(range(m),repeat=r):
   p=perm(r,m,t,z); allw.append(p)
   x=0;seen=set()
   while x not in seen:seen.add(x);x=p[x]
   cyclic=len(seen)==n
   check(cyclic==(gcd(t,r)==1 and gcd(sum(z),m)==1),'cycle criterion')
   if cyclic:
    H=subgroup(p)
    rho=0 if g==1 else sum(z)*pow(t,-1,g)%g
    if H in parents:check(parents[H]==rho,'slope generator invariant')
    parents[H]=rho
 K=m**(r-1)*phi(g)//g
 check(len(parents)==K,'parent count')
 for rho in range(g):
  if gcd(rho,g)==1:check(sum(v==rho for v in parents.values())==m**(r-1)//g,'slope count')
 orbits=edge_orbits(parents,n);d=(r+m-2)//2
 check(len(set(orbits.values()))==d,'common pair rank')
 natural=list(range(n));selector=[i+r*j for i in alpha(r) for j in alpha(m)]
 for label,P,qexpected in [('natural',natural,1),('selector',selector,d)]:
  check(sorted(P)==list(range(n)) and (P[0],P[-1])==(0,n-1),'path endpoints')
  touched={orbits[tuple(sorted((a,b)))] for a,b in zip(P,P[1:])}
  check(len(touched)==qexpected,'pin rank')
  # Enumerate product banks, orient each distance pair by its lower positive representative.
  legal=0
  for bits in product((0,1),repeat=d):
   def arc(x,y):
    i,j=x%r,x//r;ii,jj=y%r,y//r
    t,delta,offset=(r,(ii-i)%r,0) if i!=ii else (m,(jj-j)%m,(r-1)//2)
    k=min(delta,t-delta)-1
    return bool(bits[offset+k])==(delta<=t//2)
   good=all(arc(a,b) for a,b in zip(P,P[1:]))
   if good:
    legal+=1
    for H in parents:
     p=next(p for p in H if order(p)==n)
     check(all(arc(x,y)==arc(p[x],p[y]) for x,y in combinations(range(n),2)),'all-parent owner')
  check(legal==2**(d-qexpected),'intersection size')
 rows.append({'r':r,'m':m,'parents':K,'slope_classes':phi(g),'parents_per_slope':m**(r-1)//g,'natural_intersection':2**(d-1),'selector_intersection':1})
print(json.dumps({'status':'pass','checks':checks,'rows':rows},indent=2))
