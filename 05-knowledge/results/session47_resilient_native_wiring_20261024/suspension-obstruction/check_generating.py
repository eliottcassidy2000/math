#!/usr/bin/env python3
"""Exact formal-series Schur determinant and cyclic-four controls."""
import argparse,itertools,json
from fractions import Fraction
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--output',type=Path);args=p.parse_args();checks=0

def req(x):
 global checks
 checks+=1
 if not x:raise ValueError('control '+str(checks))
def mm(A,B):return [[sum(x*y for x,y in zip(row,col)) for col in zip(*B)] for row in A]
def trace(A):return sum(A[i][i] for i in range(len(A)))
def ident(n):return [[int(i==j) for j in range(n)] for i in range(n)]
def Tmat(n,b):
 A=[[0]*n for _ in range(n)]
 for i in range(n-1):A[i][i+1]=1
 for k,(i,j) in enumerate((i,j) for i in range(n) for j in range(i+2,n)):A[i][j]=(b>>k)&1;A[j][i]=1-A[i][j]
 return A
def matrices(T):
 n=len(T);K=[r[:] for r in T]
 for i in range(n):
  for j in range(i+2,n):K[i][j]^=1;K[j][i]^=1
 u=[int(i!=0) for i in range(n)];v=[int(i==n-1) for i in range(n)];w=[1-x for x in v]
 H=[r+[u[i],v[i]] for i,r in enumerate(K)]+[[int(i==0) for i in range(n)]+[0,0],w+[1,0]]
 return K,H,u,v,w
def powertr(A,r):
 P=ident(len(A));out=[len(A)]
 for _ in range(r):P=mm(P,A);out.append(trace(P))
 return out
def detseries(A,R):
 tr=powertr(A,R);d=[1]
 for k in range(1,R+1):
  num=-sum(tr[i]*d[k-i] for i in range(1,k+1));req(num%k==0);d.append(num//k)
 req(all(x==0 for x in d[len(A)+1:]));return d
stats={'triangle_four_cases':0,'series_cases':0}
for n in range(2,7):
 for bits in range(1<<((n-1)*(n-2)//2)):
  T=Tmat(n,bits);K,H,u,v,w=matrices(T);K2=mm(K,K);H2=mm(H,H)
  boundary=trace(mm(H2,H2))-trace(mm(K2,K2));rhs=sum(K2[0][j]*u[j] for j in range(n))+sum(w[i]*K2[i][n-1] for i in range(n))+K[0][n-1]
  req(boundary==4*rhs);stats['triangle_four_cases']+=1
  if n>5:continue
  R=n+6;P=ident(n);aa=[];bb=[];cc=[];dd=[]
  for j in range(R+1):
   aa.append(sum(P[0][i]*u[i] for i in range(n)));bb.append(P[0][n-1]);cc.append(sum(w[i]*P[i][j]*u[j] for i in range(n) for j in range(n)));dd.append(sum(w[i]*P[i][n-1] for i in range(n)));P=mm(P,K)
  q=[1]+[0]*R
  for k in range(2,R+1):q[k]-=aa[k-2]+dd[k-2]
  for k in range(3,R+1):q[k]-=bb[k-3]
  for k in range(4,R+1):q[k]+=sum(aa[i]*dd[k-4-i]-bb[i]*cc[k-4-i] for i in range(k-3))
  dk=detseries(K,R);dh=detseries(H,R)
  req([sum(dk[j]*q[k-j] for j in range(k+1)) for k in range(R+1)]==dh)
  # -z d(log q) is the difference in closed-walk traces.
  tw=powertr(H,R);tk=powertr(K,R)
  req(all(-k*q[k]==sum(q[j]*(tw[k-j]-tk[k-j]) for j in range(k)) for k in range(1,R+1)))
  stats['series_cases']+=1
out={'status':'PASS','checks':checks,**stats,'scope':'All1099 fixed-forward-path tournaments n2..6 for fourth-cycle boundary; all75 n2..5 for exact formal determinant and closed-walk series through n+6.'}
s=json.dumps(out,indent=2,sort_keys=True)+'\n'
if args.output:args.output.write_text(s)
else:print(s,end='')
