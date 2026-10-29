"""Exact conditional laws under all binary bracketings, independent producer-free test."""
from fractions import Fraction as F
from itertools import product
from collections import defaultdict
import json

def require(b,msg):
 if not b:raise RuntimeError(msg)
def ds(n):return [d for d in range(1,n+1) if n%d==0]
def trees(ids):
 if len(ids)==1:return [ids[0]]
 return [(a,b) for k in range(1,len(ids)) for a in trees(ids[:k]) for b in trees(ids[k:])]
def counttree(t,factors):
 if isinstance(t,int):return {d:{(d,):F(1)} for d in ds(factors[t])}
 a,b=map(lambda u:counttree(u,factors),t);counts={d:len(law) for d,law in a.items()};counts2={d:len(law) for d,law in b.items()}
 out=defaultdict(dict);total=defaultdict(int)
 for x in a:
  for y in b:total[x*y]+=counts[x]*counts2[y]
 for x,px in a.items():
  for y,py in b.items():
   weight=F(counts[x]*counts2[y],total[x*y])
   for u,pu in px.items():
    for v,pv in py.items():out[x*y][u+v]=weight*pu*pv
 return out
checks=0
for factors in product((2,3,4,6),repeat=4):
 direct=defaultdict(list)
 for tup in product(*(ds(n) for n in factors)):
  d=1
  for x in tup:d*=x
  direct[d].append(tup)
 for t in trees(tuple(range(4))):
  actual=counttree(t,factors)
  require(set(actual)==set(direct),'support')
  for d,tuples in direct.items():
   require(actual[d]=={tup:F(1,len(tuples)) for tup in tuples},'uniform conditional law');checks+=1
for p in (2,3,5,7,11):
 require(F(p*p)-F(p,2)**2==F(3*p*p,4),'non-tensor determinant');checks+=1
print(json.dumps({'factor_rows':256,'binary_shapes':5,'fibre_checks':checks,'status':'pass'}))
