#!/usr/bin/env python3
"""Exact standard-library checks; no numerical eigenvalues or Smith package."""
from fractions import Fraction as F
from itertools import product
from math import gcd,lcm
import json,random
COUNT=0

def check(test,msg):
 global COUNT
 COUNT+=1
 if not test:raise RuntimeError(msg)
def eye(n):return [[int(i==j) for j in range(n)] for i in range(n)]
def trans(A):return list(map(list,zip(*A)))
def neg(A):return [[-x for x in row]for row in A]
def add(A,B):return [[x+y for x,y in zip(a,b)]for a,b in zip(A,B)]
def mul(A,B):return [[sum(x*y for x,y in zip(a,b)) for b in zip(*B)]for a in A]
def scale(c,A):return [[c*x for x in row]for row in A]
def block(A,B,C,D):return [a+b for a,b in zip(A,B)]+[c+d for c,d in zip(C,D)]
def inv(A):
 n=len(A);T=[[F(x)for x in row]+list(map(F,e))for row,e in zip(A,eye(n))]
 for j in range(n):
  i=next(i for i in range(j,n)if T[i][j]);T[j],T[i]=T[i],T[j]
  z=T[j][j];T[j]=[x/z for x in T[j]]
  for i in range(n):
   if i!=j:
    z=T[i][j];T[i]=[x-z*y for x,y in zip(T[i],T[j])]
 return [row[n:]for row in T]
def det(A):
 n=len(A);T=[row[:]for row in A];old=1;sgn=1
 for k in range(n-1):
  if not T[k][k]:
   opts=[i for i in range(k+1,n)if T[i][k]]
   if not opts:return 0
   i=opts[0];T[k],T[i]=T[i],T[k];sgn=-sgn
  pivot=T[k][k]
  for i in range(k+1,n):
   for j in range(k+1,n):T[i][j]=(T[i][j]*pivot-T[i][k]*T[k][j])//old
  for i in range(k+1,n):T[i][k]=0
  old=pivot
 return sgn*T[-1][-1] if n else 1

def pf(A):
 n=len(A)
 if n==0:return 1
 return sum((-1)**(j+1)*A[0][j]*pf([[A[a][b]for b in range(1,n)if b!=j]for a in range(1,n)if a!=j])for j in range(1,n))
def tournament(n,mask):
 A=[[0]*n for _ in range(n)];b=0
 for i in range(n):
  for j in range(i+1,n):A[i][j]=1 if mask>>b&1 else -1;A[j][i]=-A[i][j];b+=1
 return A

def rankmod(A,p):
 T=[[x%p for x in row]for row in A];r=0
 for j in range(len(T[0])):
  opts=[i for i in range(r,len(T))if T[i][j]]
  if not opts:continue
  i=opts[0];T[i],T[r]=T[r],T[i];z=pow(T[r][j],-1,p);T[r]=[(x*z)%p for x in T[r]]
  for i in range(len(T)):
   if i!=r:
    z=T[i][j];T[i]=[(x-z*y)%p for x,y in zip(T[i],T[r])]
  r+=1
  if r==len(T):break
 return r

def denom(A):return lcm(*(x.denominator for row in A for x in row))
# Full even tournament Pfaffian universe through six vertices.
for n in [2,4,6]:
 for mask in range(1<<(n*(n-1)//2)):
  A=tournament(n,mask);p=pf(A)
  check(p%2==1,('odd pf',n,mask))
  check(det(A)==p*p,('det square',n,mask))
  check(rankmod(A,2)==n,('mod2 rank',n,mask))
# All tournaments through four vertices; deterministic samples at five/six.
rng=random.Random(960447);cases=[]
for n in range(1,5):cases += [tournament(n,m)for m in range(1<<(n*(n-1)//2))]
for n in [5,6]:cases += [tournament(n,rng.randrange(1<<(n*(n-1)//2)))for _ in range(16)]
for A in cases:
 n=len(A);I=eye(n);Z=[[0]*n for _ in range(n)];J2=[[0,1],[-1,0]]
 U=[[-1,1]for _ in range(n)];C=mul(neg(J2),trans(U))
 P=block(A,U,neg(trans(U)),J2)
 Q=block(I,[[0,0]for _ in range(n)],C,eye(2))
 check(mul(trans(Q),mul(P,Q))==block(A,[[0,0]for _ in range(n)],[[0]*n for _ in range(2)],J2),'padding congruence')
 D=block(A,add(A,I),add(A,neg(I)),neg(A));B=add(I,scale(-2,mul(A,A)))
 check(trans(D)==neg(D),'double skew')
 check(mul(D,D)==block(neg(B),Z,Z,neg(B)),'double square')
 check(det(D)==det(B),'double determinant')
 Pi=[I[i]+add(scale(2,A),I)[i]for i in range(n)]
 check(mul(Pi,D)==[[-x for x in B[i]]+B[i]for i in range(n)],'quotient relations')
 # Original block columns map -B and B, so their sum eliminates y.
 R=inv(D);link=scale(-1,mul(A,inv(B)))
 check([row[:n]for row in R[:n]]==link,'transported pairing')
 check(trans(link)==neg(link),'link alternating')
 check(all(x.denominator==1 for row in mul(B,link)for x in row),'link left integral')
 if n%2==0:
  Ai=inv(A);ones=[[1]*n for _ in range(n)];K=block(A,ones,neg(ones),A)
  q=mul(Ai,[[1]for _ in range(n)]);qq=mul(q,trans(q))
  check(inv(K)==block(Ai,qq,neg(qq),Ai),'join inverse')
  check(det(K)==det(A)**2,'join order square')
  d=denom(Ai);e=denom(inv(K));check(e%d==0 and d*d%e==0,'exponent bounds')
  check(pf(P)==pf(A),'padding pf')
  for _ in range(3):
   u=[[rng.choice([-1,1])]for _ in range(n)];v=[[rng.choice([-1,1])]for _ in range(n)];eps=rng.choice([-1,1]);UV=[a+b for a,b in zip(u,v)]
   Border=block(A,UV,neg(trans(UV)),scale(eps,J2))
   contraction=mul(trans(u),mul(Ai,v))[0][0]
   check(pf(Border)==pf(A)*(eps+contraction),'general border')
M=[[0,1,1,1],[-1,0,-1,1],[-1,1,0,-1],[-1,-1,1,0]]
K=block(M,[[1]*4 for _ in range(4)],[[-1]*4 for _ in range(4)],M)
check((det(M),denom(inv(M)),4-rankmod(M,3))==(9,3,2),'C3 square')
check((det(K),denom(inv(K)),8-rankmod(K,3))==(81,9,2),'C9 square hostile')
print(json.dumps({'status':'PASS','checks':COUNT,'full_even_tournaments':sum(1<<(n*(n-1)//2)for n in [2,4,6]),'transport_cases':len(cases),'self_join_hostile':{'old_order':9,'old_exponent':3,'new_order':81,'new_exponent':9,'new_3_rank':2}},sort_keys=True))
