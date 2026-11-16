"""Exact finite calibration + constructive 2^s·5^d controller.
Generate return_automaton.json first with the earlier eight/braid lattice.py.
A plan uses pairs (return-block index, repetition count), not a literal expansion.
"""
from pathlib import Path
from itertools import permutations,combinations,product
from collections import deque
import argparse,json
I=tuple(int(i==j) for i in range(4) for j in range(4))
def need(x):
 if not x:raise AssertionError("certificate failure")
def mm(A,B,m):
 return tuple(sum(A[4*i+k]*B[4*k+j] for k in range(4))%m for i in range(4) for j in range(4))
def apply(A,x,m):
 return tuple(sum(a*b for a,b in zip(A[4*i:4*i+4],x))%m for i in range(4))
def dot(a,b,m):return sum(x*y for x,y in zip(a,b))%m
def powm(A,k,m):
 Y=I;A=tuple(a%m for a in A)
 while k:
  if k&1:Y=mm(Y,A,m)
  A=mm(A,A,m);k>>=1
 return Y
def rank(rows,p):
 rows=[list(row) for row in rows]
 if not rows:return 0
 k=0
 for j in range(len(rows[0])):
  z=next((i for i in range(k,len(rows)) if rows[i][j]%p),None)
  if z is None:continue
  rows[k],rows[z]=rows[z],rows[k]
  pivot=pow(rows[k][j]%p,-1,p);rows[k]=[v*pivot%p for v in rows[k]]
  for i in range(k+1,len(rows)):
   c=rows[i][j]
   if c%p:rows[i]=[(a-c*b)%p for a,b in zip(rows[i],rows[k])]
  k+=1
  if k==len(rows):break
 return k
def charpoly(A):
 C=[1]
 for k in range(1,5):
  d=0
  for ii in combinations(range(4),k):
   for P in permutations(range(k)):
    a=(-1)**sum(P[i]>P[j] for i in range(k) for j in range(i+1,k))
    for i,j in enumerate(P):a*=A[4*ii[i]+ii[j]]
    d+=a
  C.append((-1)**k*d)
 return C
def ord5(X):
 t=96720
 need(powm(X,t,5)==I)
 for p in (2,3,5,13,31):
  while t%p==0 and powm(X,t//p,5)==I:t//=p
 return t
def all5(X,start,goal=None):
 start=tuple(a%5 for a in start)
 todo=deque([start]); path={start:()}
 while todo:
  a=todo.popleft()
  if a==goal:return path[a]
  for i,B in X:
   b=apply(B,a,5)
   if b not in path:path[b]=path[a]+(i,);todo.append(b)
 return path if goal is None else None
def calibrate(D):
 A=[tuple(M) for M in D['matrices']]
 h=tuple(D['readout']);z=tuple(D['start'])
 iR=D['words'].index([4,3,5]);T=A[iR]
 need(charpoly(T)==[1,-905,1960,9424,3456] and powm(T,620,25)==I)
 U=[i for i,B in enumerate(A) if rank([B[4*j:4*j+4] for j in range(4)],5)==4]
 need(len(U)==139 and iR in U)
 R=powm(T,16,200);v=(4,0,0,3)
 need(dot(h,v,8)==1 and dot(h,apply(R,z,8),8)%2==1)
 need(all(R[4*i+j]%8==v[i]*sum(h[k]*R[4*k+j] for k in range(4))%8 for i in range(4) for j in range(4)))
 scal={dot(h,apply(R,apply(A[i],v,8),8),8):i for i in U}
 need(set(scal)=={1,3,5,7})
 G={i:powm(mm(R,A[i],25),16,25) for i in U}
 X=list(G.items())
 orbit=all5(X,apply(R,z,5))
 need(len(orbit)==624 and (0,0,0,0) not in orbit)
 orders={i:ord5(G[i]) for i in U}
 C=[]
 for i in U:
  W=powm(G[i],orders[i],25)
  need(all((a-b)%5==0 for a,b in zip(W,I)))
  C.append(tuple((a-b)//5%5 for a,b in zip(W,I)))
 need(rank(C,5)==16)
 return A,h,z,iR,U,scal,orders
def plan(D,s,d,N,start=None,C=None):
 if C is None:C=calibrate(D)
 A,h,z,iR,U,scal,orders=C;T=A[iR]
 s=max(3,s);q=2**s;f=5**d;M=q*f
 need(N%2 and d>=1 and d<=30 and s<=30)
 x=tuple(z if start is None else start)
 K=16+620*max(0,(3*s-16+619)//620)
 R=powm(T,K,M);ans=[(iR,K)]
 x=apply(R,x,M)
 ne=tuple(R[4*i+j]%q for i in range(4) for j in range(4))
 cols=[tuple(ne[4*i+j] for i in range(4)) for j in range(4)]
 v=next(y for y in cols if dot(h,y,q)%2)
 v=tuple(y*pow(dot(h,v,q),-1,q)%q for y in v)
 c=dot(h,apply(R,v,q),q)
 need(c%2 and apply(R,v,q)==tuple(c*a%q for a in v) and any(a%5 for a in x))
 g={i:mm(R,A[i],M) for i in U}
 def do(i,k):
  nonlocal x
  if k:
   # compressed macro (i, k, iR, K): repeat [i, reset^K] k times
   ans.append((i,k,iR,K));x=apply(powm(g[i],k,M),x,M)
 i3,i5=scal[3],scal[5]
 al,be=(dot(h,apply(g[i],v,q),q) for i in (i3,i5))
 goal=N*pow(dot(h,x,q),-1,q)%q
 if goal%4!=1:do(i3,1);goal=goal*pow(al,-1,q)%q
 val=1
 for k in range(q//4):
  if val==goal:do(i5,k);break
  val=val*be%q
 else:raise AssertionError("dyadic unit")
 E=2**max(s-2,4)
 G={i:powm(g[i],E,M) for i in U}
 w=(int((N-h[3])*pow(h[0],-1,f)%f),0,0,1)
 route=all5(list(G.items()),x,tuple(a%5 for a in w))
 need(route is not None)
 for i in route:do(i,E)
 P={i:powm(G[i],orders[i],25) for i in U}
 for j in range(1,d):
  b=5**j
  delta=tuple((a-c)//b%5 for a,c in zip(w,x))
  need(all((a-c)%b==0 for a,c in zip(w,x)))
  if not any(delta):continue
  rows=[];idx=[]
  for i,T0 in P.items():
   F=tuple((a-c)//5%5 for a,c in zip(T0,I))
   y=apply(F,x,5)
   if rank(rows+[y],5)>len(rows):rows.append(y);idx.append(i)
   if len(rows)==4:break
  need(len(rows)==4)
  a=next(c for c in product(range(5),repeat=4) if tuple(sum(c[k]*rows[k][l] for k in range(4))%5 for l in range(4))==delta)
  for i,mul in zip(idx,a):do(i,E*orders[i]*5**(j-1)*mul)
 need(dot(h,x,M)==N%M)
 return ans,x
def load():
 here=Path(__file__).resolve().parent
 paths=[here/'return_automaton.json',here.parent/'eight_braid_root'/'return_automaton.json',here.parent/'chanpak_eight_braid_20261116'/'return_automaton.json']
 for p in paths:
  if p.exists():return json.loads(p.read_text())
 raise FileNotFoundError('Run lattice.py in the eight/braid packet first.')
def main():
 D=load();C=calibrate(D)
 print('624-state primitive orbit; rank-16 tangent check passed.')
 x=None
 for s,d,N in [(10,3,1023),(10,3,1001),(10,3,1707)]:
  ans,x=plan(D,s,d,N,x,C)
  print('Target',N,'mod',2**s*5**d,'compressed instructions',len(ans),'output',dot(C[1],x,2**s*5**d))
if __name__=='__main__':main()
