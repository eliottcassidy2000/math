from itertools import combinations
from collections import defaultdict
import json

def strong(r):
 n=len(r);a=[x|1<<i for i,x in enumerate(r)]
 for k in range(n):
  for i in range(n):
   if a[i]>>k&1:a[i]|=a[k]
 return all(x==(1<<n)-1 for x in a)
def clone(r,a,m):
 old=[u for u in range(len(r)) if u!=a];labels=old+[a]*m;n=len(labels);out=[0]*n
 for i in range(n):
  for j in range(i+1,n):
   forward=labels[i]==labels[j] or r[labels[i]]>>labels[j]&1
   out[i if forward else j]|=1<<(j if forward else i)
 return out
def counts(r):
 n=len(r);N=1<<n;D=[[0]*n for _ in range(N)];Z=[[0]*n for _ in range(N)]
 for i in range(n):D[1<<i][i]=1
 Z[1][0]=1
 for S in range(1,N):
  for i in range(n):
   if not D[S][i]:continue
   Q=r[i]&~S
   while Q:
    b=Q&-Q;Q-=b;j=b.bit_length()-1
    D[S|b][j]+=D[S][i];Z[S|b][j]+=Z[S][i]
 return sum(D[-1]),sum(Z[-1][i] for i in range(n) if r[i]&1) if n>=3 else 0

def edges(r,a):return [(o,i) for o in range(len(r)) for i in range(len(r)) if r[a]>>o&1 and r[i]>>a&1 and r[o]>>i&1]
def two(E):return any(a!=c and b!=d for (a,b),(c,d) in combinations(E,2))
