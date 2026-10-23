from functools import lru_cache
from itertools import combinations
from math import comb
import json

def ordinary_matrices(n,k):
 @lru_cache(None)
 def rec(caps):
  total=sum(caps)
  if total==0:return 1
  rows=total//k
  if total%k or max(caps)>rows:return 0
  answer=0
  for pos in combinations([i for i,c in enumerate(caps) if c],k):
   cs=list(caps)
   for i in pos:cs[i]-=1
   answer+=rec(tuple(sorted(cs,reverse=True)))
  return answer
 return rec((k,)*n)

def symmetric_matrices(n,k):
 @lru_cache(None)
 def rec(caps):
  m=len(caps)
  if not m:return 1
  if min(caps)<0 or max(caps)>m:return 0
  d=caps[0];rest=caps[1:];answer=0
  for diag in [0,1]:
   q=d-diag
   if q<0 or q>m-1:continue
   for pos in combinations([i for i,c in enumerate(rest) if c],q):
    cs=list(rest)
    for i in pos:cs[i]-=1
    answer+=rec(tuple(sorted(cs,reverse=True)))
  return answer
 return rec((k,)*n)
rows=[]
for p in [3,5,7]:
 layers=[]
 for k in range(1,(p-1)//2+1):
  A=ordinary_matrices(p,k);S=symmetric_matrices(p,k)
  numerN=A+(p-1)**2*comb(p,k);numerF=S+(p-1)*comb((p-1)//2,k//2)
  if numerN%(p*p) or numerF%p:raise ValueError('nonintegral Burnside')
  N=numerN//(p*p);F=numerF//p
  if (N+F)%2:raise ValueError('merged parity')
  layers.append(dict(k=k,matrices=A,symmetric_matrices=S,ordinary=N,self_converse=F,merged=(N+F)//2))
 rows.append(dict(p=p,layers=layers,total_ordinary=1+sum(r['ordinary'] for r in layers),total_self_converse=1+sum(r['self_converse'] for r in layers),total_merged=1+sum(r['merged'] for r in layers)))
if rows[1]['total_ordinary']!=97 or rows[1]['total_self_converse']!=31 or rows[1]['total_merged']!=64:raise ValueError('five-cell control')
print(json.dumps(dict(status='PASS',method='independent remaining-column-capacity DP and symmetric-row DP, exact integers',rows=rows),indent=2))
