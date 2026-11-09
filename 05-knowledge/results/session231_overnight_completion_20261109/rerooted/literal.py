# Literal subset DP copied from the sibling interval-placement verifier.
import json
from pathlib import Path

def need(c,m):
 if not c: raise RuntimeError(m)
def ordinal(blocks):
 n=sum(map(len,blocks)); out=[]; off=0
 for b in blocks:
  later=((1<<n)-1)^((1<<(off+len(b)))-1)
  out.extend((x<<off)|later for x in b);off+=len(b)
 return out

def close(rows,cut):
 n=len(rows); v=n;r=n+1
 return [x|(1<<v)|((0 if cut>>i&1 else 1)<<r) for i,x in enumerate(rows)]+[1<<r,cut]

def graft(rows,r,b):
 n=len(rows);m=len(b); add=((1<<m)-1)<<n
 return [x|(0 if i==r else add) for i,x in enumerate(rows)]+[(x<<n)|(1<<r) for x in b]

def state(rows,r):
 n=len(rows);N=1<<n; hp=[0]*N;hp[0]=1;end=[None]*N;end[0]=[0]*n;start=[0]*N
 for s in range(1,N):
  ds=[0]*n
  for v in range(n):
   if s>>v&1:
    t=s^(1<<v)
    ds[v]=sum(end[t][u] for u in range(n) if rows[u]>>v&1) if t else 1
  end[s]=ds;hp[s]=sum(ds)
 # Number of root-start paths on a set by first-step recurrence, using reversed endpoint DP.
 rev=[sum((rows[j]>>i&1)<<j for j in range(n)) for i in range(n)]
 se=[None]*N;se[0]=[0]*n
 for s in range(1,N):
  ds=[0]*n
  for v in range(n):
   if s>>v&1:
    t=s^(1<<v)
    ds[v]=sum(se[t][u] for u in range(n) if rows[v]>>u&1) if t else 1
  se[s]=ds;start[s]=ds[r]
 full=N-1
 return (hp[full],end[full][r],sum(start[s]*hp[full^s] for s in range(N)),hp[full^(1<<r)])

def hq(b):
 # add a sink root: its U is H, and direct subset HP reused below for Q.
 n=len(b);N=1<<n;d=[[0]*n for _ in range(N)];hs=[1]+[0]*(N-1)
 for s in range(1,N):
  for v in range(n):
   if s>>v&1:
    t=s^(1<<v);d[s][v]=sum(d[t][u] for u in range(n) if b[u]>>v&1) if t else 1
  hs[s]=sum(d[s])
 return hs[-1],sum(hs[s]*hs[N-1-s] for s in range(N))
