#!/usr/bin/env python3
"""Verify independently via vertex-level start-DP and reachability, not ear formulas."""
from functools import lru_cache

def verify(expected, parent, cut, n=11):
 a=[[0]*n for _ in range(n)]; bit=0
 for i in range(n-1):
  for j in range(i+1,n-1):
   a[i][j]=(parent>>bit)&1; a[j][i]=1-a[i][j]; bit+=1
 for i in range(n-1):
  a[i][n-1]=(cut>>i)&1; a[n-1][i]=1-a[i][n-1]
 @lru_cache(None)
 def f(v,mask):
  if mask==0:return 1
  return sum(f(w,mask^(1<<w)) for w in range(n) if mask>>w&1 and a[v][w])
 h=sum(f(v,((1<<n)-1)^(1<<v)) for v in range(n))
 reach=[row[:] for row in a]
 for k in range(n):
  for i in range(n):
   for j in range(n): reach[i][j]|=reach[i][k]&reach[k][j]
 strong=all(reach[i][j] for i in range(n) for j in range(n))
 if h!=expected or not strong: raise ValueError((expected,parent,cut,h,strong))
 return h,strong
if __name__=='__main__':
 import sys
 for line in sys.stdin:
  h,c,cut=map(int,line.split()); print(h,c,cut,*verify(h,c,cut))
