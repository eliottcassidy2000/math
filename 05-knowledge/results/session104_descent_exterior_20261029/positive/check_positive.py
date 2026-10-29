#!/usr/bin/env python3
"""Exact K4 and endpoint birth/filter controls. No external dependencies."""
from fractions import Fraction as F
from itertools import combinations
import json
CHECKS=0
def need(x,msg):
 global CHECKS
 CHECKS+=1
 if not x: raise RuntimeError(msg)
def rank(A):
 A=[[F(x) for x in row] for row in A];r=0
 for c in range(len(A[0])):
  p=next((i for i in range(r,len(A)) if A[i][c]),None)
  if p is None:continue
  A[r],A[p]=A[p],A[r];q=A[r][c];A[r]=[x/q for x in A[r]]
  for i in range(len(A)):
   if i!=r:
    q=A[i][c];A[i]=[x-q*y for x,y in zip(A[i],A[r])]
  r+=1
 return r
def tags(V,d):
 return {(v,k,e):((F(k)+e*d)/v)%1 for v in V for k in range(v) for e in (-1,1)}
def safe(t,V,d):
 return all(min((v*t)%1,1-(v*t)%1)>=d for v in V)
def accepted(V,d):return {a:t for a,t in tags(V,d).items() if safe(t,V,d)}
def main():
 dag=cyclic=0;edges=list(combinations(range(4),2))
 for bits in range(64):
  arcs=[(u,v) if bits>>i&1 else (v,u) for i,(u,v) in enumerate(edges)]
  B=[[int(v==i)-int(u==i) for u,v in arcs] for i in range(4)]
  need(rank(B)==3,'K4 incidence rank')
  remaining=set(range(4));order=[]
  while remaining:
   sources=[v for v in remaining if not any(b==v and a in remaining for a,b in arcs)]
   if not sources:break
   v=min(sources);order.append(v);remaining.remove(v)
  triangle=next(((a,b,c) for a,b,c in __import__('itertools').permutations(range(4),3) if (a,b) in arcs and (b,c) in arcs and (c,a) in arcs),None)
  if not remaining:
   dag+=1;h={v:i for i,v in enumerate(order)}
   need(all(h[v]-h[u]>=1 for u,v in arcs),'strict positive separator')
   need(triangle is None,'DAG has no triangle')
  else:
   cyclic+=1;need(triangle is not None,'tournament cyclic implies triangle')
   a,b,c=triangle;z=[int(e in [(a,b),(b,c),(c,a)]) for e in arcs]
   need(sum(z)==3 and all(sum(x*y for x,y in zip(row,z))==0 for row in B),'positive circulation erased')
 need((dag,cyclic)==(24,40),'orientation census')
 cases=0;increments=0
 for q in [4,7,14]:
  d=F(1,q)
  for m in range(1,5):
   for V in combinations(range(1,9),m):
    direct=accepted(V,d);old={}
    for j,w in enumerate(V):
     old={a:t for a,t in old.items() if safe(t,[w],d)}
     born={a:t for a,t in tags([w],d).items() if safe(t,V[:j],d)}
     need(not(set(old)&set(born)),'tag disjointness')
     old.update(born);need(old==accepted(V[:j+1],d),'incremental equals direct');increments+=1
    need(old==direct,'final')
    E=sorted(set(tags(V,d).values())|{F(0),F(1)})
    samples=E[:-1]+[(a+b)/2 for a,b in zip(E,E[1:])]
    need(bool(direct)==any(safe(t,V,d) for t in samples),'endpoint witness equivalence')
    cases+=1
 V=tuple(range(1,14));d=F(1,14);a=accepted(V,d);locations=set(a.values())
 need(len(tags(V,d))==182 and len(a)==12 and len(locations)==6,'tight AP atoms')
 need(locations=={F(j,14) for j in range(14) if __import__('math').gcd(j,14)==1},'coprime locations')
 need(all(sum(t==s for t in a.values())==2 for s in locations),'two owners per location')
 up=[len(accepted(V,F(1,4))) for V in [(1,),(1,3)]]
 need(up==[2,6],'birth term blocks count descent')
 print(json.dumps({'status':'PASS','require_checks':CHECKS,'K4_orientations':64,'positive_incidence_faithful':dag,'positive_triangle_erasure':cyclic,'endpoint_rows':cases,'incremental_steps':increments,'tight_AP':{'tags':182,'accepted_tags':12,'locations':6},'nonmonotone_atom_counts':up},indent=2,sort_keys=True))
if __name__=='__main__':main()
