#!/usr/bin/env python3
"""Finite independent checks for extension.md; checks remain active under -O."""
import itertools, math
from check import tournament, profile, cycle, hamilton_dp

def require(condition):
    if not condition: raise AssertionError('check failed')

def ordinal(A,B):
    n,m=len(A),len(B)
    E=[[False]*(n+m) for _ in range(n+m)]
    for i in range(n):
        for j in range(n): E[i][j]=A[i][j]
        for j in range(m): E[i][n+j]=True
    for i in range(m):
        for j in range(m): E[n+i][n+j]=B[i][j]
    return E

def bin_paths(E):
    n=len(E); count=0
    for assignment in range(1<<n):
        choices=[]
        for bit in (0,1):
            vertices=[v for v in range(n) if (assignment>>v)&1==bit]
            choices.append(sum(all(E[p[i]][p[i+1]] for i in range(len(p)-1))
                               for p in itertools.permutations(vertices)))
        count+=choices[0]*choices[1]
    return count

def prime(x):
    return x>=2 and all(x%d for d in range(2,math.isqrt(x)+1))

X=tournament(4,4)
require(profile(X)[1:3]==[5,22]); require(bin_paths(X)==32)
for n in range(1,4):
 for m in range(1,4):
  for x in range(1<<(n*(n-1)//2)):
   for y in range(1<<(m*(m-1)//2)):
    A=tournament(n,x); B=tournament(m,y); ab=ordinal(A,B)
    require(bin_paths(ab)==bin_paths(A)*bin_paths(B))
    require(hamilton_dp(ab)==hamilton_dp(A)*hamilton_dp(B))
unit=tournament(1,0)
for n in (0,1,2,4,8):
    B=ordinal(X,tournament(n,0)) if n else X
    value=hamilton_dp(cycle(unit,B,unit))
    require(value==5+2**(n+5))
require(prime(37) and prime(19) and prime(1171))
for t in range(40):
    p=(2**(3*t+7)+5)//7
    require(7*p==2**(3*t+7)+5)
    require(p%7==(5+2*t)%7)
print('PASS: two-bin state, strong 37 and 7·19, 7·1171 constructions, and 7-adic identity')
