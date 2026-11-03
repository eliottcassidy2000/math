from math import isqrt
from itertools import combinations
from collections import defaultdict

def square(n): return n>=0 and isqrt(n)**2==n

def construct(k):
 for q in range(1,10000):
  x=[z*q*q for z in (2,39,119,254,303)]
  a,b=321*q*q,322*q*q
  chosen=[]
  for t in range(25*q+1,isqrt(643*q*q-1)+1):
   y=t*t-b
   if not all(not square(y+v) for v in [*x,*chosen,a]):continue
   chosen.append(y)
   if len(chosen)==k: break
  if k==0:chosen=[]
  if len(chosen)==k:
   vs=x+chosen+[a,b]; n=len(vs)
   actual={(i,j) for i,j in combinations(range(n),2) if square(vs[i]+vs[j])}
   expect={(i,n-1) for i in range(n-2)}|{(0,2),(0,3)}
   if actual!=expect:raise Exception((q,actual,expect))
   return q,vs
 raise RuntimeError('not found')

def graph(vs,j):
 n=len(vs); b=n-1
 back={(i,b) for i in range(n-2)}
 if j>=1: back.add((0,2))
 if j>=2: back.add((0,3))
 return {(i,h) for i in range(n) for h in range(n) if i!=h and ((i<h and (i,h) not in back) or (i>h and (h,i) in back))}

def count(vs,j):
 n=len(vs); b=n-1
 back={(i,b) for i in range(n-2)}
 if j>=1: back.add((0,2))
 if j>=2: back.add((0,3))
 dp=[defaultdict(int) for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for mask in range(1,1<<n):
  for i,c in dp[mask].items():
   for h in range(n):
    if mask>>h&1:continue
    if (i<h and (i,h) not in back) or (i>h and (h,i) in back):dp[mask|1<<h][h]+=c
 return sum(dp[-1].values())
for k in range(5):
 q,vs=construct(k)
 if len(vs)!=len(set(vs)) or vs!=sorted(vs):raise Exception('labels not distinct and ordered')
 edges=[graph(vs,j) for j in range(3)]
 if edges[0]^edges[1]!={(0,2),(2,0)} or edges[1]^edges[2]!={(0,3),(3,0)}:raise Exception('wrong edge flips')
 for ed in edges:
  for start in range(len(vs)):
   reached={start}
   while True:
    more=reached|{v for u,v in ed if u in reached}
    if more==reached:break
    reached=more
   if len(reached)!=len(vs):raise Exception('not strongly connected')
 vals=[count(vs,j) for j in range(3)]
 expected=[1+32*2**k+j*(2+16*2**k) for j in range(3)]
 if vals!=expected:raise Exception((k,vals,expected))
 print('k=',k,'q=',q,'vertices=',vs,'counts=',vals)
