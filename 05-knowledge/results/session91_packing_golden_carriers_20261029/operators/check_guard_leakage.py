#!/usr/bin/env python3
"""Exact hard-constraint leakage ranks and the golden Lucas-five operator.
Standard library only; runtime checks survive -O.
"""
from fractions import Fraction
from math import gcd
from itertools import combinations
from random import Random
import json
CHECKS=0
def need(x,msg='check failed'):
 global CHECKS
 CHECKS+=1
 if not x:raise AssertionError(msg)
def mm(A,B):
 if not A:return []
 return [[sum(a*b for a,b in zip(r,c)) for c in zip(*B)] for r in A]
def transpose(A,cols=0):return list(map(list,zip(*A))) if A else [[] for _ in range(cols)]
def rank(A,p=None):
 if not A:return 0
 a=[[Fraction(x) if p is None else x%p for x in row] for row in A];r=0
 for j in range(len(a[0])):
  k=next((k for k in range(r,len(a)) if a[k][j]),None)
  if k is None:continue
  a[r],a[k]=a[k],a[r];v=a[r][j]
  a[r]=[x/v if p is None else x*pow(int(v),-1,p)%p for x in a[r]]
  for i in range(r+1,len(a)):
   v=a[i][j]
   if v:a[i]=[x-v*y if p is None else (x-v*y)%p for x,y in zip(a[i],a[r])]
  r+=1
  if r==len(a):break
 return r
def parts(n,edges):
 adj=[set() for _ in range(n)]
 for u,v in edges:adj[u].add(v);adj[v].add(u)
 unseen=set(range(n));b=0;out=[]
 while unseen:
  s=min(unseen);color={s:0};todo=[s];ok=True;unseen.remove(s)
  while todo:
   u=todo.pop()
   for v in adj[u]:
    if v not in color:color[v]=1-color[u];todo.append(v);unseen.remove(v)
    elif color[v]==color[u]:ok=False
  out.append((color,ok));b+=ok
 return b,out,adj
def matrices(n,edges):
 V=[s for s in range(1<<n) if all(not(s>>u&1 and s>>v&1) for u,v in edges)]
 W=[s for s in range(1<<n) if s not in set(V)]
 A=[[int((s^t).bit_count()==1) for t in V] for s in V]
 L=[[int((s^t).bit_count()==1) for t in W] for s in V]
 D=mm(L,transpose(L,len(W)))
 if not W:D=[[0]*len(V) for _ in V]
 return V,W,A,L,D
def cycle(n):return sorted({tuple(sorted((i,(i+1)%n))) for i in range(n) if n>1})
def fib(n):
 a,b=0,1
 for _ in range(n):a,b=b,a+b
 return a
def audit_graph(n,edges,extra=False):
 V,W,A,L,D=matrices(n,edges);b,components,adj=parts(n,edges)
 need(rank(L)==len(V)-(1<<b),'boundary rank')
 need(rank(D)==len(V)-(1<<b),'Gram rank')
 if n<=4:
  weights=[(-1 if i%2 else 1)*(i+1) for i in range(n)]
  WL=[[weights[(s^t).bit_length()-1] if (s^t).bit_count()==1 else 0 for t in W] for s in V]
  need(rank(WL)==len(V)-(1<<b),'nonzero named weight rank')
  for j,t in enumerate(W):
   entries=[(i,WL[i][j]) for i in range(len(V)) if WL[i][j]]
   if len(entries)==2:
    vals=[]
    for i,c in entries:
     product=1
     for h in range(n):
      if V[i]>>h&1:product*=weights[h]
     vals.append(c*product)
    need(vals[0]==vals[1],'source-weight gauge')

 # Independent direct full-cube two-step moment formula.
 A2=mm(A,A)
 for i,s in enumerate(V):
  for j,t in enumerate(V):
   full=n if s==t else 2 if (s^t).bit_count()==2 else 0
   need(full-A2[i][j]==D[i][j],'moment defect')
 # Every nonzero forbidden column has one or two predecessors.
 for col in transpose(L,len(W)):
  need(sum(col) in (0,1,2),'conflict star predecessor bound')
 if extra:
  for p in (3,5,7):need(rank(L,p)==len(V)-(1<<b),'odd characteristic incidence rank')
 return len(V),rank(D)

def determinant(A):
 a=[list(map(Fraction,r)) for r in A];out=Fraction(1)
 for j in range(len(a)):
  k=next((k for k in range(j,len(a)) if a[k][j]),None)
  if k is None:return 0
  if k!=j:a[j],a[k]=a[k],a[j];out=-out
  v=a[j][j];out*=v
  for i in range(j+1,len(a)):
   c=a[i][j]/v
   for h in range(j,len(a)):a[i][h]-=c*a[j][h]
 need(out.denominator==1,'integer determinant');return int(out)
def smith(A):
 n=len(A);div=[1]
 for k in range(1,n+1):
  d=0
  for rows in combinations(range(n),k):
   for cols in combinations(range(n),k):d=gcd(d,abs(determinant([[A[i][j] for j in cols] for i in rows])))
  div.append(d)
 return [div[k]//div[k-1] for k in range(1,n+1)]

def main():
 census=[]
 for n in range(6):
  pairs=list(combinations(range(n),2));num=0
  for bits in range(1<<len(pairs)):
   edges=[e for j,e in enumerate(pairs) if bits>>j&1];audit_graph(n,edges);num+=1
  census.append(num)
 rng=Random(9107)
 for n in (6,7,8):
  pairs=list(combinations(range(n),2))
  for _ in range(24):audit_graph(n,[e for e in pairs if rng.randrange(2)],True)
 cyclic=[]
 for n in range(3,11):
  V,W,A,L,D=matrices(n,cycle(n));E=sum(map(sum,A))//2
  need(len(V)==fib(n-1)+fib(n+1),'Lucas dimension')
  need(E==n*fib(n-1),'Lucas edges')
  need(sum(D[i][i] for i in range(len(V)))==n*(len(V)-2*fib(n-1)),'Lucas leakage trace')
  need(rank(D)==len(V)-1-(n%2==0),'Lucas leakage rank')
  if n%2==0:
   alt=[sum(1<<i for i in range(n) if i%2==a) for a in (0,1)]
   for mask in alt:
    row=D[V.index(mask)]
    need(row[V.index(mask)]==n//2 and sum(x!=0 for x in row)==1,'alternating amplitude recovery')
   need(all(x==0 for x in D[V.index(0)]),'vacuum Gram kernel')
   v=[sum((-1)**i for i in range(n) if mask==1<<i) for mask in V]
   need(all(sum(a*b for a,b in zip(row,v))==0 for row in D),'alternating singleton kernel')
  cyclic.append(dict(n=n,states=len(V),edges=E,rank=rank(D)))
 V,W,A,L,D=matrices(5,cycle(5));levels=[[i for i,s in enumerate(V) if s.bit_count()==k] for k in (0,1,2)]
 blocks=[]
 for k,inds in enumerate(levels):
  block=[[D[i][j] for j in inds] for i in inds];blocks.append(block)
  if k:
   diag=k+1
   need(all(block[i][i]==diag for i in range(5)),'golden block diagonal')
   off=[[block[i][j] if i!=j else 0 for j in range(5)] for i in range(5)]
   need(all(sum(r)==2 for r in off),'golden block cycle degree')
   # Minimal polynomial of C5 adjacency: (A-2I)(A²+A-I)=0.
   I=[[int(i==j) for j in range(5)] for i in range(5)];O2=mm(off,off)
   U=[[off[i][j]-2*I[i][j] for j in range(5)] for i in range(5)]
   Z=[[O2[i][j]+off[i][j]-I[i][j] for j in range(5)] for i in range(5)]
   need(all(x==0 for row in mm(U,Z) for x in row),'C5 minimal polynomial')
 need(all(D[i][j]==0 for a in range(3) for b in range(3) if a!=b for i in levels[a] for j in levels[b]),'level decomposition')
 C=blocks[2];C2=mm(C,C);C3=mm(C2,C)
 need(all(x%5==1 for row in C2 for x in row),'mod5 square all ones')
 need(all(x%5==0 for row in C3 for x in row),'mod5 nilpotent depth3')
 need(rank(C)==5 and rank(C,5)==3 and rank(C2,5)==1,'mod5 rank profile')
 need(rank(L,5)==10 and rank(D,5)==8,'incidence versus Gram rank')
 need(rank(D,2)==9,'mod2 Gram rank')
 snf=[smith(blocks[k]) for k in (1,2)]
 need(snf==[[1,1,1,1,4],[1,1,1,5,25]],'exact Smith minors')
 residues={s:sum(pow(4,i,11) for i in range(5) if s>>i&1)%11 for s in V}
 need(len(set(residues.values()))==11,'golden reader bijection')
 QR={1,3,4,5,9};arcs=[]
 for s in V:
  for i in range(5):
   t=s^(1<<i)
   if t in residues and not(s>>i&1):
    need((residues[t]-residues[s])%11 in QR,'positive insertion Paley arc')
    arcs.append((residues[s],residues[t]))
 need(len(arcs)==15,'guarded arc count')
 need(QR|{(-x)%11 for x in QR}==set(range(1,11)),'ambient Cayley complete')
 need(sum(A[V.index(0)])==5 and sum(A[V.index(5)])==2,'degree hostile')
 print(json.dumps(dict(status='PASS',checks=CHECKS,all_graph_counts=census,random_graphs=72,cycles=cyclic,golden5=dict(states=11,edges=15,leakage_rank=10,leakage_trace=25,mod5_Gram_rank=8,block_smith=snf,mod5_boundary_rank=10,pair_block_mod5_ranks=[3,1,0],residues=sorted(residues.items()),positive_arcs=arcs)),indent=2))
if __name__=='__main__':main()
