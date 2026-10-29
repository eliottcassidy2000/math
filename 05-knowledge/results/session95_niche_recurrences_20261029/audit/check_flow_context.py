"""No-Schur independent cross-matrix enumeration for the degree12 context."""
from collections import defaultdict
from functools import lru_cache
import json
@lru_cache(None)
def weak(total,n):
 if n==1:return ((total,),)
 return tuple((j,)+tail for j in range(total+1) for tail in weak(total-j,n-1))
edges=[(1,0),(2,0),(3,0),(4,0),(2,1),(3,1),(4,1),(3,2),(2,4),(4,3)]
kernels=[]
for degree in range(10):
 out=defaultdict(int)
 for weights in weak(degree,10):
  div=[0]*5
  for w,(a,b) in zip(weights,edges):div[a]+=w;div[b]-=w
  out[tuple(div)]+=1
 kernels.append(out)
answers=[];cross_triples=0
for reverse in (False,True):
 coefficients=[]
 for D in range(13):
  value=int(D%3==0) # t0: one C3 internal circulation; T2 internalzero.
  for t in range(1,D//3+1):
   remain=D-3*t
   for matrix in weak(t,10):
    rows=(sum(matrix[:5]),sum(matrix[5:]));cols=tuple(matrix[i]+matrix[i+5] for i in range(5))
    for a in weak(t,5):
     div=tuple(c-x for c,x in zip(cols,a))
     if reverse:div=tuple(-x for x in div)
     for b in weak(t,2):
      j=b[0]-rows[0]
      if j<0 or j>remain:continue
      value+=kernels[remain-j].get(div,0);cross_triples+=1
  coefficients.append(value)
 answers.append(coefficients)
if answers[0][:12]!=answers[1][:12] or [v[12] for v in answers]!=[9470,9462]:raise RuntimeError(answers)
print(json.dumps({'status':'PASS','checks':2,'through_degree':12,'forward':answers[0],'reverse':answers[1],'cross_triples_evaluated':cross_triples,'producer_imports':False},sort_keys=True))
