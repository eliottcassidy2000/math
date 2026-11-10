#!/usr/bin/env python3
"""Permutation/split path covers, independent of producer subset DP."""
from itertools import permutations,product
import json,sys

def need(x,m):
 if not x:raise RuntimeError(m)
def path(a,p):return all(a[u]>>v&1 for u,v in zip(p,p[1:]))
def counts(a):
 n=len(a);H=0;ends=[0]*n;starts=[0]*n;W=[[0]*n for _ in a]
 for p in permutations(range(n)):
  if path(a,p):H+=1;ends[p[-1]]+=1
  for j in range(n):
   if path(a,p[:j]) and path(a,p[j:]):
    starts[p[j]]+=1
    if j:W[p[j-1]][p[j]]+=1;W[p[j]][p[j-1]]+=1
 return H,ends,starts,W
if '--fail' in sys.argv:need(False,'intentional sentinel')
a=[10,28,49,36,41,3];h,_,A,W=counts(a)
need(h==43 and A==[55,64,72,55,78,57],'base profiles')
expected=[[0,26,14,20,26,28],[26,0,20,18,20,20],[14,20,0,26,28,24],[20,18,26,0,24,22],[26,20,28,24,0,22],[28,20,24,22,22,0]]
need(W==expected,'interaction matrix')
def F(c):return h+2*sum(A[i] for i in range(6) if c>>i&1)-sum(W[i][j] for i in range(6) for j in range(i+1,6) if c>>i&1 and c>>j&1)
def closure(c):
 b=[x|(1<<6) for x in a]+[1<<7,0]
 for i in range(6):
  if c>>i&1:b[7]|=1<<i
  else:b[i]|=1<<7
 return b
cert=[]
for c in (37,49,13,25):
 H,ends,starts,_=counts(closure(c));need(H==F(c) and ends[7]==h,'closure count/end')
 need(starts[7]==h+2*sum(A[i] for i in range(6) if c>>i&1),'closure U')
 cert.append([c,H,starts[7],ends[7]])
need(cert==[[37,345,411,43],[49,347,423,43],[13,347,407,43],[25,349,419,43]],'certificate')
# Inequality (5), recounting every deletion graph by permutations.
for i in range(6):
 for j in range(i+1,6):
  vs=[v for v in range(6) if v not in (i,j)]
  hdel=sum(path(a,p) for p in permutations(vs))
  need(W[i][j]>=2*hdel,'positive lower bound')
# Bounded-digit iff, exhaustive finite hostile controls, zero weights included.
weight_tests=0
for weights in product(range(5),repeat=4):
 s={sum(w*b for w,b in zip(weights,bits)) for bits in product(range(2),repeat=4)}
 t=0;ok=True
 for w in sorted(w for w in weights if w):
  if w>t+1:ok=False
  t+=w
 need((s==set(range(t+1)))==ok,'complete sequence iff');weight_tests+=1
print(json.dumps(dict(status='PASS',core_path_count=h,closure_certificates=cert,weight_tests=weight_tests,scope='literal stated core and four physical eight-vertex closures; no coverage expansion'),indent=2,sort_keys=True))
