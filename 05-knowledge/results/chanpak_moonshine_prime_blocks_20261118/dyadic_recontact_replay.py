from itertools import combinations
from collections import defaultdict

def h(out):
 n=len(out);D=[[0]*n for _ in range(1<<n)]
 for i in range(n):D[1<<i][i]=1
 for S,row in enumerate(D):
  for i,val in enumerate(row):
   todo=out[i]&~S
   if not val:continue
   while todo:
    b=todo&-todo;D[S|b][b.bit_length()-1]+=val;todo^=b
 return sum(D[-1])
def out(n,z):
 A=[0]*n
 for j,(a,b) in enumerate(combinations(range(n),2)):
  if z>>j&1:A[a]|=1<<b
  else:A[b]|=1<<a
 return A

def twin(A,v,k):
 labels=[i for i in range(len(A)) if i!=v]+[v]*k
 return [[int(i<j) if a==b else A[a][b] for j,b in enumerate(labels)] for i,a in enumerate(labels)]

import json
records=[]
for k in range(1,9):
 A0=out(5,70);A=[[int(A0[i]>>j&1) for j in range(5)] for i in range(5)];B=twin(A,0,k);base=[sum(v<<j for j,v in enumerate(row)) for row in B];n=len(base);target=6*2**k;best=(0,0,0);anchored=(0,0,0)
 for apex in range(n):
  spokes=[v for v in range(n) if v!=apex];w=n-1
  def ev(mask):
   C=base.copy()
   for i,v in enumerate(spokes):
    if mask>>i&1:C[apex]^=1<<v;C[v]^=1<<apex
   return h(C)
  c=ev(0);b=[ev(1<<i)-c for i in range(w)];J={(i,j):ev(1<<i|1<<j)-c-b[i]-b[j] for i,j in combinations(range(w),2)}
  D=[c+sum(b[i] for i in range(w) if s>>i&1)+sum(v for (i,j),v in J.items() if s>>i&1 and s>>j&1) for s in range(1<<w)]
  for d in [1<<i for i in range(w)]:
   for fixed in [0,d]:
    live=[s for s in range(1<<w) if s&d==fixed and D[s^d]-D[s]==2];seen=set()
    for s in live:
     if s in seen:continue
     todo=[s];seen.add(s);C=[]
     while todo:
      u=todo.pop();C.append(u)
      for v in live:
       if v not in seen and (u^v).bit_count()<=2 and abs(D[u]-D[v])<=4:seen.add(v);todo.append(v)
     lo=min(D[u] for u in C);hi=max(D[u] for u in C)+2;rec=((hi-lo)//2+1,lo,hi)
     best=max(best,rec)
     if lo<=target-1 and hi>=target+1:anchored=max(anchored,rec)
 print('k',k,'n',n,'center',target,'best',best,'center-pair-containing',anchored,flush=True)
 records.append((k,n,target,best,anchored))
expected=[(4,9,15)]+[(2,6*2**k-1,6*2**k+1) for k in range(2,9)]
if [e[-1] for e in records] != expected: raise ValueError(records)
