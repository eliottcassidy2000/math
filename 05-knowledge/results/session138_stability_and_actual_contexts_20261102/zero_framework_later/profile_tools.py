"""Standalone actual ordered-cover profiles; backward H and anchored partitions."""
import math

def covers(a):
 n=len(a);N=1<<n;E=[[0]*n for _ in range(N)];hs=[0]*N;hs[0]=1
 for S in range(1,N):
  for i in range(n):
   if not(S>>i&1):continue
   R=S^(1<<i)
   E[S][i]=sum(E[R][j] for j in range(n) if R>>j&1 and a[j]>>i&1) if R else 1
  hs[S]=sum(E[S])
 G=[[0]*(n+1) for _ in range(N)];G[0][0]=1
 for S in range(1,N):
  z=S&-S;R=S^z;U=R
  while True:
   A=U|z;B=S^A;weight=hs[A]
   for r in range(1,S.bit_count()+1):G[S][r]+=weight*G[B][r-1]
   if U==0:break
   U=(U-1)&R
 return [G[-1][r]*math.factorial(r) for r in range(1,n+1)]

