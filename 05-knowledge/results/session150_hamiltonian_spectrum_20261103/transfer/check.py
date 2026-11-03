from itertools import permutations,combinations,product
from math import factorial,comb
from collections import defaultdict

def tournament(n,m):
 E=list(combinations(range(n),2))
 A=[[0]*n for _ in range(n)]
 for i,(a,b) in enumerate(E):
  if m>>i&1:a,b=b,a
  A[a][b]=1
 return A

def profile(A):
 n=len(A); c=[0]*(n+1)
 for p in permutations(range(n)):
  bad=sum(not A[a][b] for a,b in zip(p,p[1:]))
  for k in range(1,n+1):
   if bad<=k-1:c[k]+=comb(n-1-bad,k-1-bad)
 return c

def substitution(A,B):
 n=len(A); C=[[0]*(n+2) for _ in range(n+2)]
 # A module then x,y where A -> x -> y -> A
 for i in range(n):
  for j in range(n):C[i][j]=A[i][j]
  C[i][n]=1; C[n+1][i]=1
 C[n][n+1]=1
 return C

for n in range(1,6):
 groups=defaultdict(dict)
 for mask in range(1<<(n*(n-1)//2)):
  A=tournament(n,mask); c=profile(A)
  groups[c[1]].setdefault(c[2] if n>=2 else 0,(mask,c))
  # direct DP rather than factorial for composite
  C=substitution(A,None); N=len(C); dp=[[0]*N for _ in range(1<<N)]
  for i in range(N):dp[1<<i][i]=1
  for s in range(1,1<<N):
   for v in range(N):
    if dp[s][v]:
     for w in range(N):
      if not s>>w&1 and C[v][w]:dp[s|1<<w][w]+=dp[s][v]
  actual=sum(dp[-1]); pred=3*c[1]+(c[2] if n>=2 else 0)
  if actual != pred: raise AssertionError((n,mask,c,actual,pred))
 amb=[(h,list(d.items())) for h,d in groups.items() if len(d)>1]
 print('n',n,'checked',1<<(n*(n-1)//2),'ambiguous',amb[:1])

def transfer(c):
 n=len(c)-1; out=[0]*(n+3)
 for k in range(1,n+1):
  for p in permutations(range(k+2)):
   # unique words with k repeated T labels, quotient via canonical labels order
   if [z for z in p if z<k]!=list(range(k)):continue
   w=['T' if z<k else ('x' if z==k else 'y') for z in p]
   bad=sum((a,b) not in [('T','x'),('x','y'),('y','T')] for a,b in zip(w,w[1:]))
   for q in range(1,n+3):
    if bad<=q-1<=k+1:out[q]+=c[k]*comb(k+1-bad,q-1-bad)
 return out
for n in range(1,5):
 for mask in range(1<<(n*(n-1)//2)):
  A=tournament(n,mask)
  if transfer(profile(A)) != profile(substitution(A,None)): raise AssertionError((n,mask))
 print('full-profile transfer checked',n)
