#!/usr/bin/env python3
from itertools import permutations,combinations,product
from collections import Counter
import json

def ck(x,m):
 if not x:raise RuntimeError(m)
def rank(M,eqs):
 p=list(range(M))
 def root(a):
  while p[a]!=a:p[a]=p[p[a]];a=p[a]
  return a
 r=0
 for a,b in eqs:
  a,b=root(a),root(b)
  if a!=b:p[a]=b;r+=1
 return r
ps=[(0,)+p+(4,) for p in permutations((1,2,3))]
E=list(combinations(range(5),2));ei={e:i for i,e in enumerate(E)}
out={}
for kind in ('vertex','edge'):
 words=ps if kind=='vertex' else [tuple(ei[tuple(sorted(e))] for e in zip(p,p[1:])) for p in ps]
 M=5 if kind=='vertex' else 10
 events=[list(zip(a,b)) for a,b in combinations(words,2)]
 coefficients=Counter()
 for mask in range(1<<len(events)):
  eqs=[e for i,E in enumerate(events) if mask>>i&1 for e in E]
  coefficients[M-rank(M,eqs)]+=(-1)**mask.bit_count()
 counts={}
 for q in (2,3):
  theory=sum(c*q**d for d,c in coefficients.items())
  direct=sum(len({tuple(cols[x] for x in w) for w in words})==len(words) for cols in product(range(q),repeat=M))
  ck(theory==direct,'inclusion-exclusion mismatch')
  counts[str(q)]=direct
 if kind=='vertex':ck(counts=={'2':0,'3':54},'vertex hostile')
 out[kind]=dict(objects=M,paths=len(words),collision_events=len(events),coloring_count_polynomial={str(d):c for d,c in sorted(coefficients.items()) if c},counts=counts)
print(json.dumps(dict(status='PASS',carrier='all six Hamiltonian paths from0 to4 in K5',results=out),indent=2))
