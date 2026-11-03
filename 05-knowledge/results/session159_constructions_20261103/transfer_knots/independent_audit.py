from itertools import combinations
from collections import Counter
from random import Random

def data(m):
 n=2*m+2; return [[i<j and (j-i!=3 or i%2!=0) or (i>j and i-j==3 and j%2==0) for j in range(n)] for i in range(n)]

def cycles(A):
 n=len(A); ret=[]
 for u in range(n):
  def dfs(path,mask):
   v=path[-1]
   if len(path)>=3 and len(path)%2 and A[v][u]:ret.append(tuple(path))
   for w in range(u+1,n):
    if not mask>>w&1 and A[v][w]:dfs(path+[w], mask | 1<<w)
  dfs([u],1<<u)
 return ret

def path_count(A):
 n=len(A); dp=[[0]*n for _ in range(1<<n)]
 for i in range(n): dp[1<<i][i]=1
 for mask in range(1,1<<n):
  for i in range(n):
   if dp[mask][i]:
    for j in range(n):
     if not mask & 1<<j and A[i][j]:dp[mask|1<<j][j]+=dp[mask][i]
 return sum(dp[-1])

def pack(cs):
 D={0:1}
 for C in cs:
  cm=sum(1<<i for i in C)
  for used,c in list(D.items()):
   if not used&cm:D[used|cm]=D.get(used|cm,0)+2*c
 return sum(D.values())

A=data(4); C=cycles(A); obs=Counter()
for cy in C:
 backward=sorted((v//2 for u,v in zip(cy,cy[1:]+cy[:1]) if u>v))
 if not backward or backward!=list(range(backward[0],backward[-1]+1)):raise ValueError(cy,backward)
 a,b=backward[0],backward[-1]
 ends={t for i in backward for t in (2*i,2*i+3)}
 cap=set(cy)-ends
 if len(cap)!=1 or next(iter(cap)) not in (2*a+1,2*b+2):raise ValueError(cy,cap)
 obs[a,b,next(iter(cap))]+=1
exp=Counter((a,b,c) for a in range(4) for b in range(a,4) for c in (2*a+1,2*b+2))
if obs!=exp:raise ValueError((obs,exp))
if pack(C)!=193:raise ValueError('packing')
print('m4 independent DFS: 20 odd cycles, correct interval/caps, weighted packing=193')
if path_count(data(7))!=7473:raise ValueError('m7 path count')
print('m7 independent subset DP =7473')
r=Random(714)
for n in range(2,8):
 for trial in range(12):
  A=[[False]*n for _ in range(n)]
  for i,j in combinations(range(n),2):
   if r.getrandbits(1):A[i][j]=True
   else:A[j][i]=True
  h=path_count(A); pc=pack(cycles(A))
  if h!=pc:raise ValueError((n,trial,h,pc))
print('independent odd-cycle identity: 72 random tournaments n=2..7')
prog=[(17,143),(133,85),(17,19),(23,17),(1015,69),(23,29),(31,23),(38073,62),(31,37),(41,31),(86,287),(41,43),(13,41)]
for A in range(4):
 for B in range(4):
  for C in range(4):
   n=13*11*2**A*3**B*5**C; steps=0
   while True:
    for u,v in prog:
     if n%v==0:n=n//v*u;steps+=1;break
    else:break
    if steps>1000:raise ValueError('nontermination', A,B,C)
   expect=13*2**(3*A+B+C)*3**A*5**B
   if n!=expect:raise ValueError((A,B,C,n,expect))
print('literal 13-fraction program: 64 independent register triples pass, including zero cases')
