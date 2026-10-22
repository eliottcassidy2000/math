#!/usr/bin/env python3
"""Exact rational tests for polynomial moment realizations."""
from fractions import Fraction as F
from math import comb, factorial
import json

def require(condition):
 if not condition: raise RuntimeError('Exact check failed')
def tr(A): return list(map(list,zip(*A)))
def mm(A,B): return [[sum(x*y for x,y in zip(row,col)) for col in zip(*B)] for row in A]
def ident(n): return [[F(i==j) for j in range(n)] for i in range(n)]
def det(A):
 A=[list(map(F,r)) for r in A];v=F(1);n=len(A)
 for k in range(n):
  p=next((i for i in range(k,n) if A[i][k]),None)
  if p is None:return F(0)
  if p!=k:A[p],A[k]=A[k],A[p];v=-v
  v*=A[k][k]
  for i in range(k+1,n):
   q=A[i][k]/A[k][k]
   for j in range(k,n):A[i][j]-=q*A[k][j]
 return v

def inertia(A):
 A=[list(map(F,r)) for r in A];p=n=z=0
 while A:
  k=len(A);diag=next((i for i in range(k) if A[i][i]),None)
  if diag is None:
   pair=next(((i,j) for i in range(k) for j in range(i+1,k) if A[i][j]),None)
   if pair is None:z+=k;break
   i,j=pair;T=ident(k);T[j][i]=1;A=mm(tr(T),mm(A,T));diag=i
  if diag:
   A[0],A[diag]=A[diag],A[0]
   for row in A:row[0],row[diag]=row[diag],row[0]
  pivot=A[0][0];p+=pivot>0;n+=pivot<0
  A=[[A[i][j]-A[i][0]*A[0][j]/pivot for j in range(1,k)] for i in range(1,k)]
 return [p,n,z]

out=[]
for d in range(1,9):
 for variant in range(4):
  coef=[F(((-1)**(j+variant))*(j+variant+1),variant+1) for j in range(d)]+[F(variant+1)]
  P=lambda n:sum(a*F(n)**j for j,a in enumerate(coef))
  dif=[P(n) for n in range(d+1)];c=[]
  while dif:c.append(dif[0]);dif=[b-a for a,b in zip(dif,dif[1:])]
  require(c[-1]==factorial(d)*coef[-1]);G=[[c[i+j] if i+j<=d else F(0) for j in range(d+1)] for i in range(d+1)]
  A=ident(d+1)
  for i in range(d):A[i+1][i]=1
  require(mm(tr(A),G)==mm(G,A))
  V=[[F(comb(i,j)) if j<=i else F(0) for i in range(d+1)] for j in range(d+1)]
  H=[[P(i+j) for j in range(d+1)] for i in range(d+1)]
  require(mm(tr(V),mm(G,V))==H)
  expected=[(d+2)//2,(d+1)//2,0]
  require(inertia(G)==expected and inertia(H)==expected)
  require(det(G)==(-1)**(d*(d+1)//2)*c[-1]**(d+1))
  v=[[F(i==0)] for i in range(d+1)]
  for n in range(2*d+7):
   require(mm(G,v)[0][0]==P(n));v=mm(A,v)
  # Full Hankel extension has no additional positive or negative squares.
  Hbig=[[P(i+j) for j in range(d+3)] for i in range(d+3)]
  require(inertia(Hbig)==expected[:2]+[2])
  out.append({'degree':d,'variant':variant,'newton_coefficients':list(map(str,c)),'inertia':expected,'determinant':str(det(G)),'moments_checked':2*d+7})
print(json.dumps({'cases':out,'case_count':len(out),'arithmetic':'Exact Fraction; failures enforced under python -O'},indent=2,sort_keys=True))
