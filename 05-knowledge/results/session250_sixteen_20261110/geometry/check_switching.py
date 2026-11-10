"""Exact switch-orbit Fourier and walk checks; no external dependencies."""
from itertools import combinations
from collections import Counter
import random,json

def h(A):
 n=len(A); d=[[0]*n for _ in range(1<<n)]
 for v in range(n):d[1<<v][v]=1
 for S in range(1,1<<n):
  for v in range(n):
   if S>>v&1:
    for w in range(n):
     if not S>>w&1 and A[v][w]:d[S|1<<w][w]+=d[S][v]
 return sum(d[-1])
def matrix(n,b):
 A=[[0]*n for _ in range(n)]
 for k,(i,j) in enumerate(combinations(range(n),2)):
  A[i][j]=b>>k&1; A[j][i]=1-A[i][j]
 return A
def orbit(A):
 n=len(A)
 return [h([[0 if i==j else A[i][j]^((s>>i^s>>j)&1) for j in range(n)]for i in range(n)])for s in range(1<<(n-1))]
def check(n,b):
 vals=orbit(matrix(n,b)); N=len(vals)
 assert sum(vals)==__import__('math').factorial(n)
 degrees=set()
 for c in range(N):
  f=sum(v*(-1 if (s&c).bit_count()%2 else 1)for s,v in enumerate(vals))
  wt=c.bit_count();wt+=wt%2 # restore nth coordinate so full character is even
  if f:degrees.add(wt)
 assert max(degrees)<=2*(n//3)
 def act(v):return [sum(v[s^(1<<i)] for i in range(n-1))+v[s^(N-1)] for s in range(N)]
 residual=vals[:]
 for j in range(n//3+1):
  nxt=act(residual);residual=[x-(n-4*j)*y for x,y in zip(nxt,residual)]
 assert not any(residual)
 if n==5:
  for s in range(N): assert sum(vals[s^(1<<i)]for i in range(n-1))+vals[s^(N-1)]==vals[s]+30
 return sorted(degrees),dict(sorted(Counter(vals).items()))
r=random.Random(250)
report={}
for n in range(2,8):
 bs=list(range(1<<(n*(n-1)//2))) if n<=5 else [0,(1<<(n*(n-1)//2))-1]+[r.getrandbits(n*(n-1)//2)for _ in range(12)]
 ds=set()
 for b in bs:
  deg,_=check(n,b);ds.update(deg)
 report[n]={'tournaments':len(bs),'fourier_degrees':sorted(ds),'operator_annihilator_residual':0}
print(json.dumps(report,indent=2))
print('regular5',check(5,sum(1<<k for k,(i,j) in enumerate(combinations(range(5),2)) if (j-i)%5<=2)))
