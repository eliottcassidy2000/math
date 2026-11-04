from itertools import combinations
from math import prod,isqrt
from pathlib import Path
import json

def check(v,msg):
 if not v:raise RuntimeError(msg)
def prime(n):return n>1 and all(n%d for d in range(2,isqrt(n)+1))
def divs(n):return [d for d in range(1,n+1) if n%d==0]
def ordmod(p,r):return next(d for d in divs(r-1) if pow(p,d,r)==1)
def good(p,q,n):return n>1 and (n%4==2 if p==q else n%2==1)
def factor(n):
 ans={};d=2
 while d*d<=n:
  while n%d==0:ans[d]=ans.get(d,0)+1;n//=d
  d+=1
 if n>1:ans[n]=ans.get(n,0)+1
 return ans
P=[p for p in range(3,25,2) if prime(p)]
count=closed=0
for size in (2,3,4):
 for S in combinations(P,size):
  for q in S:
   if q%4!=1:continue
   count+=1
   for p in S:
    O={r:ordmod(p,r) for r in S if r!=p}
    C={n for n in range(2,max(S)) if good(p,q,n) and all(any(o==d for o in O.values()) for d in divs(n) if d>1)}
    check(len(C)<=size-1,'count bound')
    for n in range(2,35):
     if not good(p,q,n):continue
     v=sum(p**k for k in range(n));rem=v
     for r in (2,)+S:
      while rem%r==0:rem//=r
     if rem==1:
      closed+=1;check(n in C,'closed row omitted')
      check(len(divs(n))<=size,'divisor bound')
check((count,closed)==(189,124),'independent census mismatch')
R=json.loads(Path(__file__).with_name('RESULTS.json').read_text())
for q,report in R['concrete_support'].items():
 S=(3,5,7,11,13);q=int(q)
 for p in S:
  O={ordmod(p,r) for r in S if r!=p}
  C=[n for n in range(2,14) if good(p,q,n) and set(divs(n)[1:])<=O]
  check(C==report['domains'][str(p)],'concrete domain mismatch')
 edges={(p,r) for p in S for r in S if p!=r and (ordmod(p,r)%4!=0 if p==q else ordmod(p,r)%2==1)}
 # Boolean Warshall independently checks reachability.
 reach={(p,p) for p in S}|edges
 for k in S:
  reach|={(i,j) for i in S for j in S if (i,k) in reach and (k,j) in reach}
 check(len(reach)==25,'strong certificate failed')
for p in (3,5,7,11,13):
 for d in range(2,7):
  F=factor(p**d-1)
  primitive=[r for r in F if r!=p and next(k for k in range(1,d+1) if pow(p,k,r)==1)==d]
  exc=d==2 and ((p+1)&p)==0
  check(bool(primitive)!=exc,'primitive divisor exception')
check(factor(sum(5**k for k in range(9)))=={19:1,31:1,829:1},'n9 factor')
check(factor(sum(5**k for k in range(5)))=={11:1,71:1},'n5 factor')
print(json.dumps({'status':'PASS','independent_support_cases':count,'independent_closed_rows':closed,'trial_factored_primitive_levels':25,'source':'Gyory-Smyth 2010 Theorem 13, page 325 (PDF page 7), checked directly'},indent=2))
