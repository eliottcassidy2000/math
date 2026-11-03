#!/usr/bin/env python3
"""Independently count one-apex witnesses via subset DP, compare prefix formula."""

def paths(bits):
 n=len(bits)+1; a=[[0]*n for _ in range(n)]
 for i in range(n-1):
  for j in range(i+1,n-1): a[i][j]=1
  if bits[i]: a[n-1][i]=1
  else: a[i][n-1]=1
 dp=[[0]*n for _ in range(1<<n)]
 for j in range(n): dp[1<<j][j]=1
 for s in range(1,1<<n):
  for v in range(n):
   if s>>v&1:
    p=s^(1<<v)
    if p: dp[s][v]=sum(dp[p][u]*a[u][v] for u in range(n) if p>>u&1)
 return sum(dp[-1])

def predicted(prefix):
 r=len(prefix); P=prefix[0]; C=0
 for i in range(1,r):
  b=prefix[i]
  C += P if b==0 else int(prefix[i-1]==0)
  P=2*P+b
 return 2+P*((1<<r)-1)+C, P, C
for r in range(1,7):
 seen=set()
 for z in range(1<<(r-1)):
  P=(1<<(r-1))+z
  prefix=[(P>>(r-i-1))&1 for i in range(r)]
  h,p,c=predicted(prefix)
  if not (p==P and 0<=c<(1<<r)-1): raise AssertionError('check failed')
  if not ((h-2)//((1<<r)-1)==P): raise AssertionError('check failed')
  if not (h<=4**r and h%2==1 and h not in seen): raise AssertionError('check failed')
  seen.add(h)
  got=paths(prefix+[0]*r)
  if not (got==h): raise AssertionError((prefix,h,got))
 print(f'r={r}, witnesses={len(seen)}, minimum={min(seen)}, maximum={max(seen)}, DP=PASS')
