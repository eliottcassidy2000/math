from itertools import combinations,permutations
from collections import defaultdict,Counter

def stats(T):
 n=len(T);N=(1<<n)-1
 def tab(T):
  D=[[0]*n for _ in range(N+1)]
  for i in range(n): D[1<<i][i]=1
  for s in range(1,N+1):
   for i,d in enumerate(D[s]):
    if not d:continue
    a=T[i]&~s
    while a:
     b=a&-a;a-=b;D[s|b][b.bit_length()-1]+=d
  return D
 E=tab(T); R=tab([N^(1<<i)^r for i,r in enumerate(T)])
 M=[[sum(E[s][u]*R[N^s][v] for s in range(1,N)) if T[u]>>v&1 else 0 for v in range(n)]for u in range(n)]
 return sum(E[N]),M,E[N],R[N]
def rows(n,x):
 A=[0]*n
 for a,(i,j) in enumerate(combinations(range(n),2)):
  u,v=(i,j) if x>>a&1 else (j,i)
  A[u]|=1<<v
 return A
def strong(T):
 n=len(T);N=(1<<n)-1
 for i in range(n):
  S=1<<i
  for _ in T:
   for j,r in enumerate(T):
    if S>>j&1:S|=r
  if S!=N:return False
 return True
if __name__=='__main__':
 for n in range(2,6):
  dist=defaultdict(set);z=None
  for x in range(1<<(n*(n-1)//2)):
   A=rows(n,x);h,M,*_=stats(A);g=tuple(sorted([M[i][j] for i in range(n) for j in range(n) if A[i]>>j&1],reverse=True))
   dist[h].add(g)
   if 0 in g and strong(A) and z is None:z=A,h,g
  print(n,dict(dist) if n<=4 else {h:len(s) for h,s in dist.items()},'strong zero',z)
