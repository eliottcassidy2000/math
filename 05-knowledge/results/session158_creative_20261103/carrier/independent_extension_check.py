#!/usr/bin/env python3
"""Independent bitmask path counts and strong-connectivity checks for carrier witnesses."""
from itertools import permutations
from math import isqrt

def instance(t):
 n=3*t+2; N=n+6
 # 0 is the left singleton, 1..n+4 the middle block, and N-1 the right singleton.
 E=[[False]*N for _ in range(N)]
 for i in range(1,N-1):
  E[0][i]=True; E[i][N-1]=True
  for j in range(i+1,N-1): E[i][j]=True
 E[N-1][0]=True
 # Vertices 1..4 induce X, all higher middle vertices form T_n.
 E[1][4]=False; E[4][1]=True
 return E

def paths(E):
 n=len(E); D=[[0]*n for _ in range(1<<n)]
 for v in range(n): D[1<<v][v]=1
 for m in range(1,1<<n):
  for v in range(n):
   if (m>>v)&1:
    prev=m^(1<<v)
    if prev: D[m][v]=sum(D[prev][u] for u in range(n) if (prev>>u)&1 and E[u][v])
 return sum(D[-1])

def strong(E):
 n=len(E)
 for root in range(n):
  seen={root}; todo=[root]
  while todo:
   u=todo.pop()
   for v in range(n):
    if E[u][v] and v not in seen: seen.add(v); todo.append(v)
  if len(seen)!=n: return False
 return True

def prime(p): return p>1 and all(p%d for d in range(2,isqrt(p)+1))
for t,expected,p in [(0,133,19),(2,8197,1171)]:
 E=instance(t); value=paths(E)
 if value!=expected or expected!=7*p or not prime(p) or not strong(E):
  raise RuntimeError((t,len(E),value,expected,p))
 print('PASS: strong',len(E),'vertex witness; H =',value,'= 7 *',p)
