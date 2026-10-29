#!/usr/bin/env python3
"""Integer-coordinate Euler repair; all checks active under -O."""
import json,math
from itertools import permutations
CHECKS=0
def need(v,msg='failed'):
 global CHECKS
 CHECKS+=1
 if not v:raise AssertionError(msg)
def padding(beta,k):
 b=(beta-1).bit_length();L=3*b+12
 return [(i,0) for i in range(L+1)]+[(3,1)]+[(6+3*j,1) for j in range(b) if k>>j&1]
def graph(points):
 g=[set() for _ in points]
 for i,p in enumerate(points):
  for j in range(i):
   q=points[j];d=sum((a-b)**2 for a,b in zip(p,q))
   need(d>=1,'minimum separation')
   if d==1:g[i].add(j);g[j].add(i)
 return g
def ahu(g):
 alive=set(range(len(g)));degree=list(map(len,g))
 while len(alive)>2:
  leaves=[i for i in alive if degree[i]==1];need(bool(leaves),'tree leaves')
  alive.difference_update(leaves)
  for i in leaves:
   for j in g[i]:degree[j]-=1
 def branch(i,p):
  data=[branch(j,i) for j in g[i] if j!=p];s=tuple(sorted(x for x,a in data));aut=math.prod(a for x,a in data)
  for x in set(s):aut*=math.factorial(s.count(x))
  return s,aut
 if len(alive)==1:
  return ('v',*branch(next(iter(alive)),-1))
 u,v=alive;a,x=branch(u,v);b,y=branch(v,u)
 return ('e',tuple(sorted((a,b))),x*y*(1+(a==b)))
def main():
 total=0
 for beta in range(1,65):
  sig=set()
  for k in range(beta):
   pts=padding(beta,k);g=graph(pts);need(sum(map(len,g))//2==len(g)-1,'tree Euler')
   s=ahu(g);need(s[-1]==1,'rigid');need(s[:-1] not in sig,'distinct');sig.add(s[:-1]);total+=1
  need(len(sig)==beta)
 # Contact convention: all-boundary square contact adds diagonal corners.
 p=padding(1,0);g=graph(p);side=sum(map(len,g))//2
 boundary=sum(max(abs(x-u),abs(y-v))==1 for i,(x,y) in enumerate(p) for u,v in p[:i])
 need((len(p),side,boundary)==(14,13,15),'corner-contact hostile')
 # A square core (C4) has beta1. Translate the rigid tree by5 to ensure gap2.
 core=[(0,0),(1,0),(1,1),(0,1)];pts=core+[(x+3,y) for x,y in p];g=graph(pts)
 need(sum(map(len,g))//2==len(g)-1,'repaired Euler')
 need(all(j<4 for i in range(4) for j in g[i]),'no cross-component contacts')
 # Graph symmetries of the square survive, but only identity preserves all distances.
 autos=[];isometries=[]
 for perm in permutations(range(4)):
  h=tuple(perm)+tuple(range(4,len(g)))
  if all((j in g[i])==(h[j] in g[h[i]]) for i in range(len(g)) for j in range(i)):
   autos.append(h)
   if all(sum((pts[i][z]-pts[j][z])**2 for z in (0,1))==sum((pts[h[i]][z]-pts[h[j]][z])**2 for z in (0,1)) for i in range(len(g)) for j in range(i)):
    isometries.append(h)
 need(len(autos)==8,'abstract D8 survives')
 need(isometries==[tuple(range(len(g)))],'global Euclidean symmetry trivial')

 print(json.dumps(dict(status='PASS',checks=CHECKS,padding_trees=total,max_beta=64,corner_control=[len(p),side,boundary]),indent=2))
if __name__=='__main__':main()
