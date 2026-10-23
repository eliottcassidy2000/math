"""Independent direct-arc and rational-matrix audit. No subject imports."""
from itertools import product
from fractions import Fraction as F
from collections import Counter
import json
checks=0

def need(c,msg):
 global checks
 checks+=1
 if not c:raise RuntimeError(msg)

def arc(u,v):
 if u==v:return 0
 if u==0 or v==19:return 1
 if v==0 or u==19:return 0
 ca=(u-1)//9;cb=(v-1)//9;a=(u-1)%9;b=(v-1)%9
 if ca==cb:return int((b-a)%9 in (1,2,3,4))
 return int((b-a)%3==0) if ca==0 else int((b-a)%3!=0)
A=[[arc(u,v) for v in range(20)] for u in range(20)]
need(all(A[u][v]+A[v][u]==1 for u in range(20) for v in range(u+1,20)),'actual tournament')

def path(a,b):return (0,)+tuple(1+(a+i)%9 for i in range(9))+tuple(10+(b+2+i)%9 for i in range(9))+(19,)
def K(P):
 kept={frozenset(e) for e in zip(P,P[1:])}
 return [[0 if u==v else A[u][v] if frozenset((u,v)) in kept else A[v][u] for v in range(20)] for u in range(20)]

def g(a,b,u):return u if u in (0,19) else 1+(u-1+a)%9 if u<10 else 10+(u-10+b)%9
states=[];K0=K(path(0,0))
for a,b in product(range(9),repeat=2):
 P=path(a,b);valid=all(A[u][v] for u,v in zip(P,P[1:]))
 need(valid==((a-b)%3==0),'complete chart validity')
 if not valid:continue
 states.append((a,b));KP=K(P)
 need(tuple(g(a,b,u) for u in path(0,0))==P,'literal observer transport')
 need(all(A[g(a,b,u)][g(a,b,v)]==A[u][v] for u,v in product(range(20),repeat=2)),'native source automorphism')
 need(all(KP[g(a,b,u)][g(a,b,v)]==K0[u][v] for u,v in product(range(20),repeat=2)),'native target transporter')
 need(P[0]==0 and P[-1]==19,'fixed endpoint owners')
 rev=P[::-1];kept={frozenset(e) for e in zip(rev,rev[1:])}
 need(all((A[v][u] if frozenset((u,v)) in kept else A[u][v])==KP[v][u] for u,v in product(range(20),repeat=2) if u!=v),'converse square with actual reversed clocks')
need(len(states)==27,'all retained lifts')

def comps(moves):
 unseen=set(states);out=[]
 while unseen:
  todo=[unseen.pop()];C=set(todo)
  for a,b in todo:
   for da,db in moves:
    p=((a+da)%9,(b+db)%9)
    need(p in states,'full-domain move')
    if p in unseen:unseen.remove(p);todo.append(p);C.add(p)
  out.append(C)
 return out
C=comps([(1,1),(-1,-1)]);need(sorted(map(len,C))==[9,9,9],'three carry components')
need(all(len({((b-a)//3)%3 for a,b in c})==1 for c in C),'actual separating carry')
need(list(map(len,comps([(1,1),(-1,-1),(3,0),(-3,0)])))==[27],'kernel move repairs connectivity')

def mm(A,B):return [[sum(a*b for a,b in zip(row,col)) for col in zip(*B)] for row in A]
def charpoly(A):
 n=len(A);power=[[F(i==j) for j in range(n)] for i in range(n)];tr=[0]
 for k in range(1,n+1):power=mm(power,A);tr.append(sum(power[i][i] for i in range(n)))
 c=[F(1)]
 for k in range(1,n+1):c.append(-sum(c[k-j]*tr[j] for j in range(1,k+1))/k)
 return c

def pm(a,b):
 c=[F(0)]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return c
M=[[F(11 if i==j else 8 if (i-j)%9 in (1,8) else 0) for j in range(9)] for i in range(9)]
want=pm(pm([1,-27],pm([1,-3],[1,-3])),pm([1,-33,171,1293],[1,-33,171,1293]))
need(charpoly(M)==want,'exact native C9 polynomial')
for e in (F(1,27),F(1,81),F(1,243)):
 P=[[F(0) for _ in range(6)] for _ in range(6)]
 for u,a in product(range(2),range(3)):
  i=3*u+a;P[i][i]=F(5,9)-e;P[i][3*(1-u)+a]=F(4,9)
  P[i][3*u+(a+1)%3]+=e/2;P[i][3*u+(a-1)%3]+=e/2
  need(sum(P[i])==1 and all(z>=0 for z in P[i]),'stochastic positive support')
 roots=[F(1),F(1,9),1-3*e/2,1-3*e/2,F(1,9)-3*e/2,F(1,9)-3*e/2]
 poly=[1]
 for root in roots:poly=pm(poly,[1,-root])
 need(charpoly(P)==poly,'all lifted eigenvalues')
 need(min(roots)>0 and 1-sorted(roots)[-2]==3*e/2,'PSD gap arbitrarily small')
print(json.dumps({'status':'PASS','checks':checks,'native_paths':len(states),'unrepaired_components':sorted(map(len,C)),'native_scaled_characteristic_polynomial':[str(c) for c in want],'scope':'Independent native and spectral hostile audit; general holonomy proof reviewed separately.'},indent=2))
