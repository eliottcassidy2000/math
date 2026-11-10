#!/usr/bin/env python3
"""Independent fixed-permutation Burnside certificate via full color matrices."""
from itertools import combinations, permutations, product
from math import factorial
from pathlib import Path
import sys,json
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'burnside'))
from palette_engines import palette_series

def require(x,label):
    if not x: raise RuntimeError(label)

results=[]
for q in range(5):
 for t in range(q+1):
  if (q-t)%2:continue
  dp,_=palette_series(15,q,t); brute=[]
  rev=tuple(c if c<t else t+((c-t)^1) for c in range(q))
  for n in range(5):
   edges=list(combinations(range(n),2)); ps=list(permutations(range(n)))
   fixed=[0]*len(ps)
   for colors in product(range(q),repeat=len(edges)):
    matrix=[[None]*n for _ in range(n)]
    for (a,b),c in zip(edges,colors): matrix[a][b]=c;matrix[b][a]=rev[c]
    for i,p in enumerate(ps):
     fixed[i]+=all(matrix[a][b]==matrix[p[a]][p[b]] for a,b in edges)
   require(sum(fixed)%factorial(n)==0,'Burnside integrality')
   count=sum(fixed)//factorial(n);brute.append(count)
   require(count==dp[n],'matrix exhaustive vs DP')
  for n in [0,1,2,7]:require(palette_series(n,q,t)[0]==dp[:n+1],'horizon consistency')
  results.append({'q':q,'t':t,'independent_brute':brute})
for args in [(-1,2,2),(3,2,1),(3,1,2),(3,-1,0),(3,3,-1),(3,2.0,0)]:
 try:palette_series(*args)
 except AssertionError:pass
 else:raise RuntimeError('invalid arguments accepted')
print(json.dumps({'status':'PASS','method':'Every coloring represented as full directed color matrix; fixedness checked for every permutation, then averaged','palettes':results},indent=2))
