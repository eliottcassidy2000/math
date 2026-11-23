from itertools import combinations
from functools import lru_cache
from random import Random

def orient(T,u,v):
 T[u]|=1<<v;T[v]&=~(1<<u)

def converse(T):
 return [z^((1<<len(T))-1)^(1<<a) for a,z in enumerate(T)]

def add(T,S):
 n=len(T);D=list(T)+[0]
 for i in range(n):orient(D,*( (n,i) if S>>i&1 else (i,n)))
 return D

def directed_cycles(T,L):
 n=len(T);M=(1<<n)-1
 def walk(root,end,S,j):
  if j==L:
   if T[end]>>root&1:yield S
   return
  a=T[end]&M&~((1<<(root+1))-1)&~S
  while a:
   bit=a&-a;a-=bit
   yield from walk(root,bit.bit_length()-1,S|bit,j+1)
 for u in range(n):yield from walk(u,u,1<<u,1)

def feature(T):
 tris=list(directed_cycles(T,3))
 return 2*sum(1 for _ in directed_cycles(T,5))+4*sum((s&t)==0 for s,t in combinations(tris,2))

def experiment(T,i,j,k):
 n=len(T);U=add(add(T,0),0)
 path=(n+1,i,n,j)
 for e,(u,v) in enumerate(zip(path,path[1:])):orient(U,*((u,v) if k>>e&1 else (v,u)))
 return U
def sign(k):return (-1)**(3-k.bit_count())
def pm(T,u,v):return 2*(T[u]>>v&1)-1
def tensors(T):
 n=len(T);M=(1<<n)-1;H=[1]*(M+1)
 @lru_cache(None)
 def ep(A,i,st):
  B=A^(1<<i)
  if not B:return 1
  return sum(ep(B,j,st) for j in range(n) if B>>j&1 and (T[i]>>j&1 if st else T[j]>>i&1))
 for A in range(1,M+1):H[A]=sum(ep(A,i,False) for i in range(n) if A>>i&1)
 @lru_cache(None)
 def q(A,m):
  if not m:return int(not A)
  ans=0;X=A
  while True:
   ans+=H[X]*q(A^X,m-1)
   if not X:return ans
   X=(X-1)&A
 def D(i,m):
  return sum((ep(A,i,True)-ep(A,i,False))*q(M^A,m-1) for A in range(1,M+1) if A>>i&1)
 return (lambda m:q(M,m)),D

r=Random(26);n0=n1=0
for n in range(3,7):
 for _ in range(6):
  T=[0]*n
  for u,v in combinations(range(n),2):orient(T,*((v,u) if r.randrange(2) else (u,v)))
  i,j=r.sample(range(n),2)
  read=sum(sign(k)*feature(experiment(T,i,j,k)) for k in range(8))
  assert read==sum(pm(T,j,jj)+1 for jj in range(n) if jj not in (i,j)) # y loses all untested old vertices
  n0+=1
  if n<=5:
   V=converse(T);_,D=tensors(T)
   for S in range(1<<n):
    q1,_=tensors(add(T,S));q2,_=tensors(add(V,S))
    for m in range(1,4):
     assert q1(m)-q2(m)==m*sum(D(i,m)+D(i,m+1) for i in range(n) if S>>i&1)
     n1+=1
T=[74,28,25,112,97,7,38]
assert all(v.bit_count()==3 for v in T)
assert [tensors(T)[1](i,1) for i in range(7)]==[0,-1,1,-1,1,1,-1]
# Cayley coordinates, left invariant orientation.
D={1,2,3,14,8,9,10,11,12,13}
def times(g,h):
 a,b=g%7,g//7;c,d=h%7,h//7
 return (a+(1<<b)*c)%7+7*((b+d)%3)
T=[sum(1<<times(g,h) for h in D) for g in range(21)]
V=converse(T)
assert all((T[i]>>j&1)+(T[j]>>i&1)==1 for i,j in combinations(range(21),2))
def sig(U):
 N=[i for i in range(21) if U[0]>>i&1]
 d={i:sum(U[i]>>j&1 for j in N) for i in N}
 a=next(i for i in N if d[i]==7)
 return a,sorted(d.values()),sorted(d[j] for j in N if U[a]>>j&1)
assert sig(T)==(1,[3,3,4,4,4,4,5,5,6,7],[3,4,4,4,5,5,6])
assert sig(V)==(6,[3,3,4,4,4,4,5,5,6,7],[3,4,4,4,4,5,6])
def all_features(U):return [feature(experiment(U,0,1,k)) for k in range(8)]
a,b=all_features(T),all_features(V)
assert sum(sign(k)*(a[k]-b[k]) for k in range(8))==2
print(n0,"two-apex mixed probes,",n1,"one-apex divergence checks.")
print("Chiral Cayley order 21, 8 feature differences:", [x-y for x,y in zip(a,b)])
