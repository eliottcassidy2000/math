#!/usr/bin/env python3
"""Independent metric adjacency and rooted-tree symmetry audit of padding output."""
import importlib.util,json,math
from pathlib import Path
N=0
def ck(x,m):
 global N
 N+=1
 if not x:raise RuntimeError(m)
def canonical_tree(g):
 alive=set(range(len(g)));deg=list(map(len,g))
 while len(alive)>2:
  leaf=[i for i in alive if deg[i]<=1];ck(bool(leaf),'acyclic leaf peeling')
  alive.difference_update(leaf)
  for i in leaf:
   for j in g[i]:deg[j]-=1
 def rooted(v,p):
  cs=[rooted(w,v) for w in g[v] if w!=p];s=tuple(sorted(c[0] for c in cs));a=math.prod(c[1] for c in cs)
  for c in set(s):a*=math.factorial(s.count(c))
  return s,a
 if len(alive)==1:
  s,a=rooted(next(iter(alive)),-1);return ('v',s),a
 u,v=alive;s,a=rooted(u,v);t,b=rooted(v,u)
 return ('e',tuple(sorted((s,t)))),a*b*(2 if s==t else 1)
def main():
 path=Path(__file__).resolve().parents[1]/'geometry'/'check_geometric_repair.py';sp=importlib.util.spec_from_file_location('producer',path);p=importlib.util.module_from_spec(sp);sp.loader.exec_module(p)
 total=0
 for beta in range(1,97):
  signatures=set()
  for k in range(beta):
   pts=p.padding(beta,k);g=[set() for _ in pts]
   ck(len(pts)==len(set(pts)),'distinct integer points')
   for i in range(len(pts)):
    for j in range(i):
     dx=pts[i][0]-pts[j][0];dy=pts[i][1]-pts[j][1];D=dx*dx+dy*dy
     ck(D>=1,'separated points')
     if D==1:g[i].add(j);g[j].add(i)
   ck(sum(map(len,g))==2*(len(g)-1),'Euler tree count')
   code,a=canonical_tree(g);ck(a==1,'independent symmetry recursion')
   ck(code not in signatures,'pairwise nonisomorphic');signatures.add(code);total+=1
 pts=p.padding(1,0)
 euclidean=sum((x-u)**2+(y-v)**2==1 for i,(x,y) in enumerate(pts) for u,v in pts[:i])
 square=sum(max(abs(x-u),abs(y-v))==1 for i,(x,y) in enumerate(pts) for u,v in pts[:i])
 ck((len(pts),euclidean,square)==(14,13,15),'side versus corner predicate')
 print(json.dumps(dict(status='PASS',checks=N,padding_trees=total,max_beta=96,corner_hostile=[14,13,15]),indent=2))
if __name__=='__main__':main()
