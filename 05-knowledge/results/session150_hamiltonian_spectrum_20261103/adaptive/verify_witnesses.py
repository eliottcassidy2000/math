#!/usr/bin/env python3
"""Independently check the witness root histories printed by search.py."""
import ast, functools, sys
p=sys.argv[1] if len(sys.argv)>1 else 'check.out'
rows=open(p).read().splitlines()
items=ast.literal_eval(next(x[5:] for x in rows if x.startswith('hits ')))
for want, roots in items:
 a=[[0,1,0],[0,0,1],[1,0,0]]
 for v in roots:
  n=len(a)
  if not 0<=v<n: raise ValueError('invalid root')
  a=[row+[int(i==v)] for i,row in enumerate(a)]+[[int(i!=v) for i in range(n)]+[0]]
 n=len(a)
 @functools.lru_cache(None)
 def ways(mask, last):
  if not mask: return 1
  return sum(ways(mask^(1<<v),v) for v in range(n) if mask>>v&1 and a[last][v])
 actual=sum(ways(((1<<n)-1)^(1<<i),i) for i in range(n))
 if actual != want: raise ValueError((want,roots,actual))
print('Verified',len(items),'witnesses in',p)
