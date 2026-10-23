#!/usr/bin/env python3
import itertools,json
from pathlib import Path
import sympy as s

def require(b,msg):
 if not b:raise RuntimeError(msg)
checks=0;universes=0
for n in range(2,5):
 allE=list(itertools.combinations(range(n),2))
 for mask in range(1<<len(allE)):
  E=[e for k,e in enumerate(allE) if mask>>k&1];m=len(E)
  x=s.symbols(f'x0:{n}');y=s.symbols(f'y0:{m}');L=s.zeros(n)
  for k,(a,b) in enumerate(E):
   w=y[k]*x[a]*x[b];L[a,a]+=w;L[b,b]+=w;L[a,b]-=w;L[b,a]-=w
  tau=s.Poly(s.expand(L[:n-1,:n-1].det()),*x)
  edgeindex={e:k for k,e in enumerate(E)}
  actual={pair:[] for pair in itertools.combinations(range(n),2)}
  for p in itertools.permutations(range(n)):
   if p[0]>p[-1]:continue
   es=[tuple(sorted((a,b))) for a,b in zip(p,p[1:])]
   if all(e in edgeindex for e in es):actual[(p[0],p[-1])].append(set(edgeindex[e] for e in es))
  for (a,b),paths in actual.items():
   mon=tuple(1 if v in (a,b) else 2 for v in range(n))
   pol=tau.coeff_monomial(mon);ones={z:1 for z in y}
   require(pol.subs(ones)==len(paths),'endpoint coefficient');checks+=1
   for k in range(m):
    require(s.diff(pol,y[k]).subs(ones)==sum(k in p for p in paths),'first edge derivative');checks+=1
   for k,l in itertools.combinations(range(m),2):
    require(s.diff(pol,y[k],y[l]).subs(ones)==sum(k in p and l in p for p in paths),'mixed edge derivative');checks+=1
  universes+=1
out={'universes':universes,'checks':checks,'scope':'all simple graphs n=2,3,4; all endpoint pairs; all first and mixed second edge derivatives','status':'PASS'}
Path(__file__).with_name('matrix-tree-projection-receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(out)
