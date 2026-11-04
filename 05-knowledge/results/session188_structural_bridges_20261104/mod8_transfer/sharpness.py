"""Literal word/cut and ambient enumeration; no producer imports."""
from itertools import permutations
from math import comb
import json
A=[48,17,59,19,32,10];B=[50,16,59,19,0,26]
def need(x,msg):
 if not x:raise RuntimeError(msg)
def tournament(a):
 return all(not (a[i]>>i&1) for i in range(len(a))) and all(((a[i]>>j)^(a[j]>>i))&1 for i in range(len(a)) for j in range(i+1,len(a)))
def q(a):
 n=len(a);r=[0]*(n+1)
 for p in permutations(range(n)):
  bad=sum(not (a[u]>>v&1) for u,v in zip(p,p[1:]));good=n-1-bad
  for k in range(1,n+1):
   if 0<=k-1-bad<=good:r[k]+=comb(good,k-1-bad)
 return r
def lift(a):
 n=len(a);return [r|(1<<n) for r in a]+[1<<(n+1),(1<<n)-1]
def hp(a):return sum(all(a[u]>>v&1 for u,v in zip(p,p[1:])) for p in permutations(range(len(a))))
need(tournament(A) and tournament(B),'physical tournament masks')
p=q(A);r=q(B);need(p==[0,11,178,828,1752,1800,720] and r==[0,3,90,612,1608,1800,720],'literal full profiles')
x=lift(A);y=lift(B);f=hp(x);g=hp(y);need((f,g)==(211,99),'literal ambient counts');need(f-g==3*(p[1]-r[1])+p[2]-r[2],'exact substitution formula');need((p[1]-r[1])%16==8 and (f-g)%16==0 and f!=g,'mod16 boundary without exact collision')
print(json.dumps({'status':'PASS','A':A,'B':B,'q_A':p,'q_B':r,'ambient_A':x,'ambient_B':y,'ambient_counts':[f,g],'core_gap':p[1]-r[1],'ambient_gap':f-g,'scope':'Failure of universal even-order mod16 odd-unit transfer, not an exact response or scalar collision.'},sort_keys=True,indent=2))
