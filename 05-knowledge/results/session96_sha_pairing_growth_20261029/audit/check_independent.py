#!/usr/bin/env python3
"""Independent finite pairing and tournament controls. No producer imports."""
from itertools import product
from fractions import Fraction as Q
from math import isqrt
import json
CHECKS=0

def ck(ok,c):
 global CHECKS
 CHECKS+=1
 if not ok:raise RuntimeError(c)

def parity(x):return x.bit_count()%2

def image(M,x):return sum(parity(row&x)<<i for i,row in enumerate(M))

def dot(x,y):return parity(x&y)

counts=[]
for n in range(1,5):
 pairs=[(i,j) for i in range(n) for j in range(i,n)];types=[0,0,0];perfect=0
 for bits in range(1<<len(pairs)):
  M=[0]*n
  for k,(i,j) in enumerate(pairs):
   if bits>>k&1:M[i]|=1<<j;M[j]|=1<<i
  ims=[image(M,x) for x in range(1<<n)]
  if len(set(ims))<1<<n:continue
  perfect+=1;diag=sum(((M[i]>>i)&1)<<i for i in range(n));c=ims.index(diag)
  eps=dot(c,image(M,c));types[0 if c==0 else 1+eps]+=1
  f=lambda x:x^(c if dot(diag,x) else 0)
  fs=[f(x) for x in range(1<<n)]
  ck(all(dot(x,image(M,x))==dot(x,image(M,c)) for x in range(1<<n)),('characteristic',n,bits))
  ck(all(f(f(x))==(x if eps==0 else f(x)) for x in range(1<<n)),('repair power',n,bits))
  ck(len(set(fs))==(1<<n if eps==0 else 1<<(n-1)),('repair image',n,bits))
  ck(all(dot(x,image(M,f(x)))==0 for x in range(1<<n)),('alternating repair',n,bits))
  if eps:
   H=[x for x in range(1<<n) if dot(x,image(M,c))==0]
   ck(set(H)==set(fs),('orthogonal projection',n,bits))
   ck(all(any(dot(x,image(M,y)) for y in H) for x in H if x),('perfect complement',n,bits))
  order=1<<n;R=isqrt(order//(2**eps))
  ck(order==2**eps*R*R,('order shape',n,bits))
 counts.append({'dimension':n,'perfect_forms':perfect,'alternating_isotropic_anisotropic_characteristic':types})
# Same square-order group and different pairing data.
ck(dot(1,image([1,2],1))==1,'I2 nonalternating')
ck(all(dot(x,image([2,1],x))==0 for x in range(4)),'H2 alternating')
ck(image([1,2],3)==3 and dot(3,3)==0,'I2 nonzero isotropic characteristic')
for e,f,R,S in product(range(2),range(2),range(1,8),range(1,8)):
 ck((2**e*R*R)*(2**f*S*S)==2**(e^f)*(2**(e*f)*R*S)**2,'binary carry')

def pf(M):
 n=len(M)
 if n==0:return 1
 return sum(((-1)**(j+1))*M[0][j]*pf([[M[a][b] for b in range(1,n) if b!=j] for a in range(1,n) if a!=j]) for j in range(1,n))

def det(M):
 A=[[Q(x) for x in row] for row in M];out=Q(1)
 for i in range(len(A)):
  pivot=next((j for j in range(i,len(A)) if A[j][i]),None)
  if pivot is None:return 0
  if pivot!=i:A[pivot],A[i]=A[i],A[pivot];out=-out
  v=A[i][i];out*=v
  for j in range(i+1,len(A)):
   t=A[j][i]/v
   for k in range(i+1,len(A)):A[j][k]-=t*A[i][k]
 return out

def inverse(M):
 n=len(M);A=[[Q(x) for x in row]+[Q(i==j) for j in range(n)] for i,row in enumerate(M)]
 for i in range(n):
  j=next(j for j in range(i,n) if A[j][i]);A[i],A[j]=A[j],A[i]
  v=A[i][i];A[i]=[x/v for x in A[i]]
  for j in range(n):
   if j!=i:
    v=A[j][i];A[j]=[x-v*y for x,y in zip(A[j],A[i])]
 return [row[n:] for row in A]

tcounts=[]
for n in (2,4,6):
 pairs=[(i,j) for i in range(n) for j in range(i+1,n)];maxmask=1<<len(pairs);dets=set()
 for mask in range(maxmask):
  M=[[0]*n for _ in range(n)]
  for k,(i,j) in enumerate(pairs):M[i][j]=1 if mask>>k&1 else -1;M[j][i]=-M[i][j]
  p=pf(M);ck(p%2==1,('odd pfaffian',n,mask));dets.add(p*p)
  if n<=4 or mask%257==0:
   ck(det(M)==p*p,('det square',n,mask));I=inverse(M)
   ck(all(I[i][j]==-I[j][i] for i in range(n) for j in range(n)),('inverse skew',n,mask))
   for x in (list(range(n)),[1]*n):
    ck(sum(Q(x[i])*I[i][j]*x[j] for i in range(n) for j in range(n))==0,'discriminant alternating')
 tcounts.append({'vertices':n,'labelled_tournaments':maxmask,'determinants':sorted(dets)})
def mm(A,B):return [[sum(x*y for x,y in zip(row,col)) for col in zip(*B)] for row in A]
def rankmod(A,p):
 A=[[x%p for x in row] for row in A];rank=0
 for j in range(len(A[0])):
  i=next((i for i in range(rank,len(A)) if A[i][j]),None)
  if i is None:continue
  A[rank],A[i]=A[i],A[rank];v=pow(A[rank][j],-1,p);A[rank]=[x*v%p for x in A[rank]]
  for i in range(len(A)):
   if i!=rank:
    v=A[i][j];A[i]=[(x-v*y)%p for x,y in zip(A[i],A[rank])]
  rank+=1
 return rank
for n in (2,4):
 pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
 for mask in range(1<<len(pairs)):
  M=[[0]*n for _ in range(n)]
  for k,(i,j) in enumerate(pairs):M[i][j]=1 if mask>>k&1 else -1;M[j][i]=-M[i][j]
  D=[M[i]+[M[i][j]+int(i==j) for j in range(n)] for i in range(n)]
  D += [[M[i][j]-int(i==j) for j in range(n)]+[-v for v in M[i]] for i in range(n)]
  M2=mm(M,M);S=[[int(i==j)-2*M2[i][j] for j in range(n)] for i in range(n)]
  ck(det(D)==det(S),('doubling determinant',n,mask))
  ID=inverse(D);IS=inverse(S);link=mm(M,IS)
  ck(all(ID[i][j]==-link[i][j] for i in range(n) for j in range(n)),('doubling transported pairing',n,mask))
  B=[M[i]+[-1,1] for i in range(n)]+[[1]*n+[0,1],[-1]*n+[-1,0]]
  ck(pf(B)==pf(M),('plus two pfaffian',n,mask))
  IB=inverse(B)
  ck(all(IB[i][j]==inverse(M)[i][j] for i in range(n) for j in range(n)),('plus two pairing',n,mask))
M=[[0,1,1,1],[-1,0,-1,1],[-1,1,0,-1],[-1,-1,1,0]]
J=[row+[1]*4 for row in M]+[[-1]*4+row for row in M]
ck(det(M)==9 and rankmod(M,3)==2,'base C3 squared')
ck(det(J)==81 and rankmod(J,3)==6,'selfjoin C9 squared not C3 fourth')
IJ=inverse(J)
ck(all((9*v).denominator==1 for row in IJ for v in row),'selfjoin exponent divides9')
ck(any((3*v).denominator!=1 for row in IJ for v in row),'selfjoin exponent not3')
# A finite-level pairing is perfect but the standard divisible injections are not isometries.
for p in (2,3,5,7):
 for n in range(1,5):
  q=p**n
  ck(Q(p*p,p*q)%1==(p*Q(1,q))%1,('scaled direct limit pairing',p,n))
  ck(Q(p*p,p*q)%1!=Q(1,q)%1,('not compatible',p,n))
print(json.dumps({'status':'PASS','checks':CHECKS,'binary_pairings':counts,'tournament_universes':tcounts},sort_keys=True))
