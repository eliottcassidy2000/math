#!/usr/bin/env python3
"""Independent C5 leakage blocks, determinantal divisors and modular powers."""
from itertools import combinations,permutations,product
from math import gcd
from independent_guard_leakage import rank_mod
import json
N=0
def ck(x,m):
 global N
 N+=1
 if not x:raise RuntimeError(m)
def mm(A,B):return [[sum(a*b for a,b in zip(row,col)) for col in zip(*B)] for row in A]
def det(A):
 n=len(A);ans=0
 for p in permutations(range(n)):
  inv=sum(p[i]>p[j] for i in range(n) for j in range(i+1,n));v=(-1)**inv
  for i in range(n):v*=A[i][p[i]]
  ans+=v
 return ans

def smith(A):
 prev=1;out=[]
 for k in range(1,len(A)+1):
  g=0
  for rows in combinations(range(len(A)),k):
   for cols in combinations(range(len(A)),k):g=gcd(g,det([[A[i][j] for j in cols] for i in rows]))
  ck(g%prev==0,'determinantal divisors');out.append(g//prev);prev=g
 return out

def main():
 edges=[(i,(i+1)%5) for i in range(5)]
 W=[s for s in range(32) if all(not(s>>i&1 and s>>j&1) for i,j in edges)]
 outside=[s for s in range(32) if s not in W]
 C=[[int((s^t).bit_count()==1) for t in outside] for s in W]
 D=mm(C,list(zip(*C)));ck(rank_mod(C,5)==10,'incidence rank5');ck(rank_mod(D,5)==8,'Gram rank5');ck(rank_mod(D,2)==9,'Gram rank2')
 blocks=[]
 for size,diag in [(1,2),(2,3)]:
  ix=[i for i,s in enumerate(W) if s.bit_count()==size];B=[[D[i][j] for j in ix] for i in ix]
  ck(all(B[i][i]==diag and sum(B[i][j] for j in range(5) if i!=j)==2 for i in range(5)),'cyclic block degrees')
  adj=[[B[i][j] if i!=j else 0 for j in range(5)] for i in range(5)]
  # Connected2-regular simple graph on5 vertices is C5.
  seen={0};todo=[0]
  while todo:
   i=todo.pop()
   for j in range(5):
    if adj[i][j] and j not in seen:seen.add(j);todo.append(j)
  ck(len(seen)==5 and all(x in(0,1) for row in adj for x in row),'actual five-cycle')
  S=smith(B);ck(S==([1,1,1,1,4] if size==1 else[1,1,1,5,25]),'exact Smith')
  blocks.append(B)
 B=blocks[1];B2=mm(B,B);B3=mm(B2,B)
 ck(all(x%5==1 for row in B2 for x in row),'square allones mod5')
 ck(all(x%5==0 for row in B3 for x in row),'cube zero mod5')
 ck([rank_mod(X,5) for X in(B,B2,B3)]==[3,1,0],'rank filtration')
 rho={s:sum(pow(4,i,11) for i in range(5) if s>>i&1)%11 for s in W}
 ck(len(set(rho.values()))==11,'golden reader')
 insertions=[]
 for s in W:
  for i in range(5):
   t=s|(1<<i)
   if t!=s and t in rho:
    ck((rho[t]-rho[s])%11 in {1,3,4,5,9},'positive Paley arc');insertions.append((s,t))
 ck(len(insertions)==15,'guarded arcs not55')
 # Nonzero weights gauge each forbidden-state equation exactly.
 for weights in[(1,-2,3,-4,5),(2,3,5,7,11)]:
  diagS=[__import__('math').prod(weights[i] for i in range(5) if s>>i&1) for s in W]
  for j,t in enumerate(outside):
   factors=[]
   for i,s in enumerate(W):
    if C[i][j]:
     v=(s^t).bit_length()-1;factors.append(weights[v]*diagS[i])
   ck(len(set(factors))<=1,'weighted column factorization')
 print(json.dumps(dict(status='PASS',checks=N,states=len(W),arcs=len(insertions),smith=[[1,1,1,1,4],[1,1,1,5,25]],mod5_powers=[3,1,0]),indent=2))
if __name__=='__main__':main()
