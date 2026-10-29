#!/usr/bin/env python3
from fractions import Fraction as Q
from itertools import combinations
import json
checks=0

def ck(v,msg):
 global checks
 checks+=1
 if not v:raise RuntimeError(msg)

def rank(A,p=None):
 a=[[Q(v) if p is None else v%p for v in row] for row in A]
 if not a:return 0
 r=0
 for c in range(len(a[0])):
  k=next((i for i in range(r,len(a)) if a[i][c]),None)
  if k is None:continue
  a[r],a[k]=a[k],a[r];x=a[r][c];iv=1/x if p is None else pow(x,-1,p)
  a[r]=[v*iv if p is None else v*iv%p for v in a[r]]
  for i in range(len(a)):
   if i==r:continue
   z=a[i][c]
   a[i]=[v-z*w if p is None else (v-z*w)%p for v,w in zip(a[i],a[r])]
  r+=1
  if r==len(a):break
 return r
bank=0
for n in range(1,6):
 edges=list(combinations(range(n),2))
 for bits in range(1<<len(edges)):
  A=[[0]*n for _ in range(n)]
  for k,(i,j) in enumerate(edges):A[i][j]=1 if bits>>k&1 else -1;A[j][i]=-A[i][j]
  ck(rank(A)==n-(n%2),('ordinary rational corank',n,bits));bank+=1
weighted=0
for n in range(2,13):
 for r in range(1,n//2+1):
  A=[[(j-i)*sum((i*j)**(2*k) for k in range(r)) for j in range(1,n+1)] for i in range(1,n+1)]
  ck(all(A[i][j]>0 for i in range(n) for j in range(i+1,n)),('same orientation',n,r))
  V=[[i**j for j in range(2*r)] for i in range(1,n+1)]
  reconstructed=[[sum(V[i][2*k]*V[j][2*k+1]-V[i][2*k+1]*V[j][2*k] for k in range(r)) for j in range(n)] for i in range(n)]
  ck(A==reconstructed,('factorization',n,r));ck(rank(V)==2*r,('Vandermonde',n,r));ck(rank(A)==2*r,('weighted rank',n,r));weighted+=1
M=[[0,1,1,1],[-1,0,-1,1],[-1,1,0,-1],[-1,-1,1,0]]
ck(rank(M)==4 and rank(M,3)==2,'same matrix different field')
A=[[j-i for j in range(4)] for i in range(4)];U=[[int(j>i)-int(j<i) for j in range(4)] for i in range(4)]
ck(rank(A)==2 and rank(U)==4,'same signs different rational corank')
print(json.dumps({'status':'PASS','checks':checks,'ordinary_tournaments':bank,'weighted_examples':weighted,'scope':'rational matrix rank, not elliptic-curve rank'},sort_keys=True))
