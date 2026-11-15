#!/usr/bin/env python3
from collections import defaultdict
from itertools import combinations
import math
def rows(n,code):
  a=[0]*n
  for k,(i,j) in enumerate(combinations(range(n),2)):
    if code>>k&1:a[i]|=1<<j
    else:a[j]|=1<<i
  return a
def dp(a):
  n=len(a); f=(1<<n)-1; ends=[[0]*n for s in range(f+1)]; h=[1]+[0]*f
  for i in range(n):ends[1<<i][i]=1
  for s in range(1,f+1):
    for i,v in enumerate(ends[s]):
      h[s]+=v
      r=a[i]&~s
      while v and r:
        z=r&-r;r-=z;ends[s|z][z.bit_length()-1]+=v
  return h,ends
def packets(a):
  n=len(a);full=(1<<n)-1
  h,end=dp(a)
  _,start=dp([full^r^(1<<i) for i,r in enumerate(a)])
  Q=sum(h[s]*h[full^s] for s in range(full+1))
  S=[[sum(h[full^s]*m[i] for s,m in enumerate(start) if s.bit_count()%2==p) for i in range(n)] for p in (0,1)]
  T=[[sum(h[full^s]*m[i] for s,m in enumerate(end) if s.bit_count()%2==p) for i in range(n)] for p in (0,1)]
  W=[[sum(end[s][i]*start[full^s][j] for s in range(1,full)) for j in range(n)] for i in range(n)]
  for mask in range(1<<n):
    on=[i for i in range(n) if mask>>i&1];off=[i for i in range(n) if not mask>>i&1]
    se,so=[sum(v[i] for i in on)for v in S]
    te,to=[sum(v[i] for i in off)for v in T]
    N=h[-1]+se+so-sum(end[-1][i] for i in on)-sum(W[i][j] for i in on for j in on)
    yield (h[-1],Q,se+so,so-se,te+to,to-te),N,mask

def run(hosts):
  buckets=defaultdict(dict)
  for n,c in hosts:
    for k,N,m in packets(rows(n,c)):buckets[k].setdefault(N,(n,c,m))
  print("bucket count",len(buckets))
  v=sorted(((k,x)for k,x in buckets.items() if len(x)>1),key=lambda t:(-len(t[1]),t[0]))
  print("forks",len(v))
  for k,x in v[:15]:print(k,sorted(x.items()))
  cop=[]
  for k,x in v:
    h=k[0]
    if h>1 and math.gcd(h,k[1])==1:
      cop.append((len({n%h for n in x})/h,k,sorted(x.items())))
  for p,k,x in sorted(cop,reverse=True)[:15]:print("primitive",p,k,x)
  return buckets
if __name__=='__main__':
  import sys
  run([tuple(map(int,line.split())) for line in open(sys.argv[1])])
