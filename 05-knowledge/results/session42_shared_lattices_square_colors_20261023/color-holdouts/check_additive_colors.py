#!/usr/bin/env python3
from math import isqrt
from collections import deque
import json

def ck(x,m):
 if not x:raise RuntimeError(m)
def rank(A,q):
 A=[r[:] for r in A];i=0
 if not A:return 0
 for j in range(len(A[0])):
  k=next((k for k in range(i,len(A)) if A[k][j]%q),None)
  if k is None:continue
  A[i],A[k]=A[k],A[i];inv=pow(A[i][j],-1,q);A[i]=[v*inv%q for v in A[i]]
  for k in range(len(A)):
   if k!=i and A[k][j]:
    f=A[k][j];A[k]=[(a-f*b)%q for a,b in zip(A[k],A[i])]
  i+=1
 return i
out=[]
for n in range(1,81):
 E=[(a,b) for a in range(n) for b in range(a+1,n) if isqrt(a+b+2)**2==a+b+2]
 A=[[int(v==a or v==b) for v in range(n)] for a,b in E];G=[[] for _ in range(n)]
 for a,b in E:G[a].append(b);G[b].append(a)
 col={};components=0;bipcomponents=0
 for s in range(n):
  if s in col:continue
  components+=1;col[s]=0;todo=[s];bip=True
  for v in todo:
   for w in G[v]:
    if w not in col:col[w]=1-col[v];todo.append(w)
    elif col[w]==col[v]:bip=False
  bipcomponents+=bip
 r2=rank(A,2);r3=rank(A,3)
 ck(r2==n-components,'F2 kernel');ck(r3==n-bipcomponents,'F3 kernel')
 # Round-trip color transport on all directed edge pairs, any assigned vertex colors.
 for q in (2,3):
  for power in range(1,5):
   c=[pow(v+1,power,q) for v in range(n)]
   for a,b in E:ck(((c[a]+c[b])%q-c[a])%q==c[b],'color transport')
 out.append(dict(n=n,components=components,bipartite_components=bipcomponents,rank2=r2,rank3=r3))
print(json.dumps(dict(status='PASS',scope='unsigned incidence rank on Q1..Q80 and exact edge transport controls',rows=out),indent=2))
