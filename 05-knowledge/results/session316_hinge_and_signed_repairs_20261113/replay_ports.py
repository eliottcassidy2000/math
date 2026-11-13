from itertools import combinations,permutations
from math import isqrt
from collections import Counter

def add(A,B,s=1):
 D=A.copy()
 for x,y in B.items():D[x]=D.get(x,0)+s*y
 return {x:y for x,y in D.items() if y}
def mul(A,B):
 D=Counter()
 for a,s in A.items():
  for b,t in B.items():
   if not a&b:D[a|b]+=s*t
 return {x:y for x,y in D.items() if y}
def const(t):return {0:t} if t else {}
def rec(V,arc,I,sign=1):
 E={(u,v):const(sign*arc(u,v)) for u in V for v in V}
 U=set(V);d={0:1}
 for k in I:
  b={1<<k:1};U.remove(k)
  d=mul(d,add({0:1},mul(E[k,k],b),-1))
  for i in U:
   for j in U:E[i,j]=add(E[i,j],mul(mul(E[i,k],b),E[k,j]))
 return {(i,j):E[i,j] for i in U for j in U},d
def sp(x,y):return int(x!=y and isqrt(x+y)**2==x+y)
def rp(R):return lambda i,j:(R[i]>>j)&1
def naive(u,v,I,arc):
 D=Counter()
 for k in range(len(I)+1):
  for J in permutations(I,k):
   p=(u,)+J+(v,)
   if all(arc(i,j) for i,j in zip(p,p[1:])):D[sum(1<<j for j in J)]+=1
 return dict(D)
def test(s,t,I,arc,sgn):
 V=list(s+t+I);E,_=rec(V,arc,I)
 assert E==rec(V,arc,I[::-1])[0]
 assert all(E[u,v]==naive(u,v,I,arc) for u in s for v in t)
 p=mul(E[s[0],t[0]],E[s[1],t[1]])
 q=mul(E[s[0],t[1]],E[s[1],t[0]])
 Z=sum(1<<i for i in I)
 assert p.get(Z,0)+q.get(Z,0)==1 and p.get(Z,0)-q.get(Z,0)==sgn
 print('signed exchange',s,t,I,sgn)
test((6,8),(13,15),(1,3),sp,1)
test((6,8),(13,15),(1,3,10),sp,-1)
R=(82,48,107,67,108,9,34); arc=rp(R)
assert all(arc(u,v)+arc(v,u)==1 for u,v in combinations(range(7),2))
assert all(arc(u,v) for u in (2,3) for v in (0,1))
test((0,1),(2,3),(4,5),arc,1)
test((0,1),(2,3),(4,5,6),arc,-1)

def inverse(D,n):
 C=add({0:1},D,-1);Q=I={0:1}
 for _ in range(n):
  Q=mul(Q,C);I=add(I,Q)
 return I
for n in range(1,6):
 V=list(range(n));pairs=list(combinations(V,2));Z=(1<<n)-1
 for b in range(1<<len(pairs)):
  R=[0]*n
  for k,(u,v) in enumerate(pairs):
   if b>>k&1:R[u]|=1<<v
   else:R[v]|=1<<u
  A=rp(R)
  _,D=rec(V,A,V);_,Dneg=rec(V,A,V,-1)
  assert (D,Dneg)==(rec(V,A,V[::-1])[1],rec(V,A,V[::-1],-1)[1])
  # product (1+x_v) has coefficient one at every subset
  H=mul(mul(Dneg,inverse(D,n)),dict.fromkeys(range(1<<n),1)).get(Z,0)
  H0=sum(all(A(u,v) for u,v in zip(P,P[1:])) for P in permutations(V))
  assert H==H0,(n,b,H,H0)
 print('ratio and both elimination orders verified',n,1<<len(pairs))
from math import factorial
for k in range(1,8):
 for p in range(k+1):
  result=0
  for pi in permutations(range(k)):
   seen=set();ncycles=0
   for u in range(k):
    if u in seen:continue
    j=u;s=0
    while j not in seen:seen.add(j);s+=(j<p);j=pi[j]
    if s%2==0:break
    ncycles+=1
   else:result+=1<<ncycles
  assert result==2*p*factorial(k-1)
print('weighted odd seam formula through k=7')
