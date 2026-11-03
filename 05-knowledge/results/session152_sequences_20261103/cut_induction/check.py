def require(x):
 if not x: raise AssertionError("cut induction check failed")
from itertools import combinations
from collections import Counter

def require(condition):
 if not condition: raise AssertionError("verification failed")

def adj(n,code):
 a=[[0]*n for _ in range(n)]
 for k,(i,j) in enumerate(combinations(range(n),2)):
  a[i][j]=(code>>k)&1;a[j][i]=1-a[i][j]
 return a

def hp(a):
 n=len(a); d=[[0]*n for _ in range(1<<n)]
 for i in range(n): d[1<<i][i]=1
 for s in range(1,1<<n):
  for i in range(n):
   if s>>i&1: d[s][i]+=sum(d[s^(1<<i)][j]*a[j][i] for j in range(n) if j!=i and s>>j&1)
 return sum(d[-1])

def ear(a,s):
 n=len(a);return [row+[int(not(s>>i&1))] for i,row in enumerate(a)]+[[s>>i&1 for i in range(n)]+[0]]

def allocate(a,K):
 n=len(a); edges=list(K); rem=a[:]; c=[[0]*n for i in range(n)]
 def go(p):
  if p==len(edges): return all(v==0 for v in rem)
  i,j=edges[p]; w=K[i,j]
  for t in range(max(0,w-rem[j]), min(w,rem[i])+1):
   rem[i]-=t;rem[j]-=w-t
   c[i][j]=t;c[j][i]=w-t
   if go(p+1):return True
   rem[i]+=t;rem[j]+=w-t
  return False
 require(go(0))
 return c

def test(n,code,strong):
 a=adj(n,code);h=hp(a); F=[hp(ear(a,s)) for s in range(1<<n)]
 A=[F[1<<i]-h for i in range(n)];K={(i,j):F[1<<i]+F[1<<j]-F[(1<<i)|(1<<j)]-h for i,j in combinations(range(n),2)}
 require(all(x%2==0 and x>=0 for x in A) and all(x%2==0 and x>=0 for x in K.values()))
 c=allocate([x//2 for x in A],{e:x//2 for e,x in K.items()})
 for s in range(1<<n):
  cut=sum(c[i][j] for i in range(n) for j in range(n) if s>>i&1 and not s>>j&1)
  require(F[s]==h+2*cut)
  if strong and s not in (0,(1<<n)-1): require(F[s]>=h+2*((n)//2))
 print('n,code,H=',n,code,h,' A=',A,'K=',K,'c=',c,'image=',sorted(set(F[1:-1])))

test(3,2,True)
test(5,8,True)
test(5,759,True)
print('PASS: direct subset DP, integral cut reconstruction, strong lower bound')
