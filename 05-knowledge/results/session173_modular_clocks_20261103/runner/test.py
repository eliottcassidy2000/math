from fractions import Fraction as F
from math import gcd,ceil
from itertools import combinations
import json,random

def require(ok, label):
 if not ok: raise RuntimeError(label)

def norm(x):
 r=x%1
 return min(r,1-r)
def cert(d,core,exc,delta,t0):
 require(all(v%d==0 and norm(v*t0)>=delta for v in core), "core validity")
 B=[{q for q in range(d) if norm(v*(t0+F(q,d)))<delta} for v in exc]
 weights=sorted(((len(B[i]&B[j]),i,j) for i,j in combinations(range(len(B)),2)),reverse=True)
 p=list(range(len(B)))
 def root(i):
  while p[i]!=i:i=p[i]
  return i
 edges=[]
 for w,i,j in weights:
  a,b=root(i),root(j)
  if a!=b:p[a]=b;edges.append((i,j,w))
 union=set().union(*B)
 bound=sum(map(len,B))-sum(w for i,j,w in edges)
 require(len(union)<=bound, "tree bound")
 good=[q for q in range(d) if q not in union]
 for q in good:require(all(norm(v*(t0+F(q,d)))>=delta for v in core+exc), "witness validity")
 caps=[gcd(v,d)*ceil(2*delta*(d//gcd(v,d))) for v in exc]
 require(all(len(b)<=c for b,c in zip(B,caps)), "gcd capacities")
 return dict(d=d,core=core,exceptions=exc,delta=str(delta),t0=str(t0),bad_sets=[sorted(b) for b in B],crude_bound=sum(caps),actual_sum=sum(map(len,B)),tree_bound=bound,union_size=len(union),good_shifts=good,witnesses=[str(t0+F(q,d)) for q in good],tree_edges=edges)
if __name__ == "__main__":
 examples=[]
 for d in [8,9,13,44,527]:
  exc=[1,3,5,7,13,44,527]
  exc=[v for v in exc if v%d]
  m=14-len(exc)
  core=[d*i for i in range(1,m+1)]
  examples.append(cert(d,core,exc,F(1,15),F(1,d*(m+1))))
 # Systematic family, with exact rational arithmetic.
 checks=0; improved=0
 for d in range(2,25):
  for den in range(3,18):
   core=[d]
   t=F(1,d*den)
   delta=F(1,den)
   exc=[i for i in range(1,25) if i%d][:min(den-2,10)]
   if not exc:continue
   c=cert(d,core,exc,delta,t);checks+=1
   improved+=c['crude_bound']>=d and c['tree_bound']<d
 print(json.dumps({'examples':examples,'systematic_checks':checks,'strict_tree_certificates_beyond_crude':improved},indent=2))
