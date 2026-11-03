#!/usr/bin/env python3
"""Independent anti-involution census using adjacency matrices and centralizers."""
from itertools import combinations, permutations, product
from math import factorial

def check(n):
 r=tuple(i^1 if i<2*(n//2) else i for i in range(n)); m=n//2; vertices=tuple(range(n))
 edges=list(combinations(vertices,2)); fixed=[]
 for bits in product((False,True),repeat=len(edges)):
  M=[[False]*n for _ in vertices]
  for (a,b),v in zip(edges,bits):M[a][b]=v; M[b][a]=not v
  if all(M[r[a]][r[b]] != M[a][b] for a,b in edges):
   fixed.append(tuple(int(M[a][b]) for a,b in edges))
 assertlen=2**(m*m+(n%2)*m)
 if len(fixed)!=assertlen:raise ValueError(('wrong fixed count', n, len(fixed),assertlen))
 central=[p for p in permutations(vertices) if all(p[r[i]]==r[p[i]] for i in vertices)]
 if len(central)!=2**m*factorial(m):raise ValueError('wrong centralizer size')
 def key(bits,p):
  M=[[False]*n for _ in vertices]
  for (a,b),v in zip(edges,bits):M[a][b]=v; M[b][a]=not v
  return tuple(M[p[a]][p[b]] for a,b in edges)
 reps={min(key(b,p) for p in central) for b in fixed}
 return len(fixed),len(central),len(reps)

for n in range(0,6):
 f,c,o=check(n); expected=(1,1,1,2,2,8)[n]
 if o!=expected:raise ValueError(('wrong class count',n,o,expected))
 print(f'n={n}: reversed={f}, centralizer={c}, self-converse classes={o}')
