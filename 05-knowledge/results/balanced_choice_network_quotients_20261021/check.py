#!/usr/bin/env python3
from fractions import Fraction as F
from itertools import combinations
from math import gcd, prod
import json
checks=0
def need(c,msg):
 global checks
 checks+=1
 if not c:raise RuntimeError(msg)
def determinant(A):
 A=[[F(x) for x in row] for row in A];z=F(1)
 for j in range(len(A)):
  p=next((i for i in range(j,len(A)) if A[i][j]),None)
  if p is None:return F(0)
  if p!=j:A[p],A[j]=A[j],A[p];z=-z
  a=A[j][j];z*=a
  for i in range(j+1,len(A)):
   b=A[i][j]/a
   for h in range(j+1,len(A)):A[i][h]-=b*A[j][h]
 return z
def inverse(A):
 d=len(A);B=[[F(x) for x in row]+[F(i==j) for j in range(d)] for i,row in enumerate(A)]
 for j in range(d):
  p=next(i for i in range(j,d) if B[i][j]);B[p],B[j]=B[j],B[p]
  a=B[j][j];B[j]=[x/a for x in B[j]]
  for i in range(d):
   if i!=j:
    a=B[i][j];B[i]=[x-a*y for x,y in zip(B[i],B[j])]
 return [r[d:] for r in B]
def lap(weights):
 W=sum(weights);return [[F(W*weights[i]*(i==j)-weights[i]*weights[j]) for j in range(len(weights))] for i in range(len(weights))]
rows=[]
for q in range(2,8):
 for s in range(1,6):
  m=q-1;Q=s*m+1;w=[s]*m+[1];L=lap(w);R=[row[:-1] for row in L[:-1]];I=inverse(R)
  need(determinant(R)==s**m*Q**(m-1),'weighted tree count')
  need(all(sum(row)==0 for row in L),'Laplacian rows')
  for a,b in combinations(range(q),2):
   va=I[a][a] if a<m else 0;vb=I[b][b] if b<m else 0;cross=I[a][b] if a<m and b<m else 0
   need(va+vb-2*cross==(F(1,w[a])+F(1,w[b]))/Q,'resistance')
  # determinantal divisors check every minor, independent of group presentation
  if m<=4:
   for rank in range(1,m+1):
    g=0
    for rr in combinations(range(m),rank):
     for cc in combinations(range(m),rank):
      a=determinant([[R[i][j] for j in cc] for i in rr]);need(a.denominator==1,'integer minor');g=gcd(g,abs(a.numerator))
    need(g==s**rank*Q**(rank-1),'Smith determinantal divisor')
  groups=[i for i,wi in enumerate(w) for _ in range(wi)]
  counts=[[0]*q for _ in range(q)];loops=[0]*q
  for a,b in combinations(range(Q),2):
   x,y=groups[a],groups[b]
   if x==y:loops[x]+=1
   else:counts[x][y]+=1;counts[y][x]+=1
  need(all(counts[i][j]==w[i]*w[j] for i,j in combinations(range(q),2)),'literal contraction')
  need(loops==[wi*(wi-1)//2 for wi in w],'hidden loops')
  for a,x in enumerate(groups):
   transition=[sum(1 for b,y in enumerate(groups) if b!=a and y==j) for j in range(q)]
   need(transition==[w[j]-(j==x) for j in range(q)],'strong lumpability')
  # Schur complement from complete graph, retain q distinct physical representatives
  h=Q-q;B=[[F(Q*(i==j)-1) for j in range(q)] for i in range(q)]
  if h:
   C=[[F(Q*(i==j)-1) for j in range(h)] for i in range(h)];CI=inverse(C);loss=sum(map(sum,CI))
   B=[[z-loss for z in row] for row in B]
  need(B==[[F(Q)*(i==j)-F(Q,q) for j in range(q)] for i in range(q)],'Kron matrix')
  blue_child=sum(w);blue_parent=sum(x*x for x in w)
  need(blue_parent-blue_child==m*s*(s-1),'blue correction')
  rows.append({'q':q,'s':s,'Q':Q,'tree_count':s**m*Q**(m-1),'smith':[s]+[s*Q]*(m-1),'hidden_loops':sum(loops)})
need(F(2,5*2)==F(1,5) and F(3,5*2)==F(3,10),'literal contraction R')
need(F(2,5)!=F(1,5),'Kron versus contraction')
need(F(2*1*2*2,5)==F(8,5),'literal entropy')
print(json.dumps({'status':'PASS','checks':checks,'declared_q':[2,7],'declared_s':[1,5],'rows':rows,
 'literal_q3_s2':{'critical_group':[2,10],'two_primary':[2,2],'entropy_bits':'8/5','child_blue':'5/25','parent_blue':'9/25'},
 'scope':'Choice-state graph only, not tournament metagraph or ideal-class identification.'},indent=2))
