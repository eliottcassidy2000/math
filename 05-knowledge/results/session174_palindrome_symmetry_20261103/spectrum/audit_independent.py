#!/usr/bin/env python3
"""Independent permutation-and-cut enumeration, unlike subset cover recurrence."""
from itertools import permutations,combinations
from math import factorial,comb
import json
import common

def require(ok,s):
 if not ok:raise RuntimeError(s)
def literal_profile(g):
 n=len(g);ordered=[0]*(n+1);H=0
 for v in permutations(range(n)):
  bad=sum(not(g[v[i]]>>v[i+1]&1) for i in range(n-1))
  H+=bad==0
  for k in range(bad+1,n+1):ordered[k]+=comb(n-1-bad,k-1-bad)
 p=[ordered[k]//factorial(k) for k in range(n+1)]
 require(all(ordered[k]%factorial(k)==0 for k in range(n+1)),'integral cover quotient')
 require(p[1]==H,'path count')
 return p
rows=[]
for b in common.B:
 p=literal_profile(b);require(p==common.cover(b),'seed profile')
 rows.append({'type':'seed','module':b,'profile':p})
for bits in range(8):
 g=[0]*3
 for k,(i,j) in enumerate(combinations(range(3),2)):
  if bits>>k&1:g[i]|=1<<j
  else:g[j]|=1<<i
 for r in [1,2]:
  a=common.wrap(g,r);p=literal_profile(a)
  require(p==common.cover(a),'wrapped profile')
  n=len(a);sigma=[3,4,5,0,1,2]+list(reversed(range(6,n)))
  fixed=0
  for v in permutations(range(n)):
   if v==tuple(sigma[x] for x in reversed(v)) and all(a[v[i]]>>v[i+1]&1 for i in range(n-1)):fixed+=1
  base=literal_profile(g)
  F=lambda u:sum(base[k]*factorial(u)//factorial(u-k) for k in range(min(len(base),u+1)))
  expected=F(2)-F(1) if r==1 else 2*F(3)-3*F(2)+3*F(1)
  require(fixed==expected,'fixed formula')
  rows.append({'type':'wrap','bits':bits,'r':r,'profile':p,'fixed':fixed})
print(json.dumps({'result':'PASS','independent_cases':rows},indent=2))
