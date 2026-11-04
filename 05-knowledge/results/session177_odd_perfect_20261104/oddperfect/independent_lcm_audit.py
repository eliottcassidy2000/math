#!/usr/bin/env python3
from math import isqrt,lcm
import json
from check_lcm_closure import closure

def need(c):
    if not c:raise RuntimeError('Independent closure verification failed')
def prime(n):return n>=2 and all(n%d for d in range(2,isqrt(n)+1))
def order(p,r):return min(d for d in range(1,r) if (r-1)%d==0 and pow(p,d,r)==1)
def geom(p,a,r):return (a+1)%r if p%r==1 else ((pow(p,a+1,r)-1)*pow(p-1,-1,r))%r
P=[p for p in range(3,44,2) if prime(p)];nhood=0;subsets=0
for p in P:
 for role in [0,1]:
  b=1+role
  K={}
  for r in P:
   if p==r:continue
   d=order(p,r)
   if (d%2==1 if not role else d%4!=0):K[r]=r if d==1 else d
  R=list(K)
  for mask in range(1<<len(R)):
   T={r for j,r in enumerate(R) if mask>>j&1};L=lcm(b,*(K[r] for r in T));C={r for r in R if L%K[r]==0}
   need(closure(T,K,b)==(frozenset(C),L))
   for r in R:
    need(C <= set(closure(T|{r},K,b)[0]))
   if C==T:
    n=L
    if n==1:
     n=max([3]+list(K.values()))+1
     while not prime(n):n+=1
    need({r for r in P if r!=p and geom(p,n-1,r)==0}==T)
   subsets+=1
  for a in range(1,1001):
   if (a%2!=0 if not role else a%4!=1):continue
   T={r for r in P if r!=p and geom(p,a,r)==0}
   need(T<=set(K));need(set(closure(T,K,b)[0])==T);nhood+=1
need(prime(4561));need([order(3,r) for r in [11,13,4561]]==[5,3,15])
print(json.dumps({'status':'PASS','subsets':subsets,'actual_neighborhoods':nhood,'forced_prime':4561},sort_keys=True))
