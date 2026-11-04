# Session180 standard-library endpoint DP, retained verbatim for reproducibility.
import itertools,json
from pathlib import Path

def response(rows):
 n=len(rows);N=1<<n;dp=[[0]*n for _ in range(N)];hs=[0]*N;hs[0]=1
 for mask in range(1,N):
  for v in range(n):
   if mask>>v&1:
    rest=mask^(1<<v)
    dp[mask][v]=1 if rest==0 else sum(dp[rest][u] for u in range(n) if rows[u]>>v&1)
  hs[mask]=sum(dp[mask])
 F=[sum(dp[s][w]*hs[(N-1)^s] for s in range(N)) for w in range(n)]
 return hs[-1],F

def strong(rows):
 n=len(rows)
 for reverse in [False,True]:
  seen=1
  while True:
   new=seen
   for v in range(n):
    if seen>>v&1:
     new|=sum(1<<u for u in range(n) if (rows[u]>>v&1)) if reverse else rows[v]
   if new==seen:break
   seen=new
  if seen!=(1<<n)-1:return False
 return True
