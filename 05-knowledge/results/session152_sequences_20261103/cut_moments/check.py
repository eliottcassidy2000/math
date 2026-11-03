def require(x):
 if not x: raise AssertionError("moment check failed")
from itertools import combinations
from math import comb
from collections import Counter

def adj(n,code):
 a=[[0]*n for _ in range(n)]
 for k,(i,j) in enumerate(combinations(range(n),2)):
  a[i][j]=(code>>k)&1; a[j][i]=1-a[i][j]
 return a

def hp(a):
 n=len(a); d=[[0]*n for _ in range(1<<n)]
 for i in range(n): d[1<<i][i]=1
 for m in range(1,1<<n):
  for i in range(n):
   if m>>i&1: d[m][i]+=sum(d[m^(1<<i)][j]*a[j][i] for j in range(n) if j!=i and m>>j&1)
 return sum(d[-1])

def check(n,code):
 a=adj(n,code); h=hp(a); ds=[]
 for m in range(1<<n):
  child=[row+[int(not(m>>i&1))] for i,row in enumerate(a)]+[[m>>i&1 for i in range(n)]+[0]]
  f=hp(child)
  if (f-h)%2: raise AssertionError("parity")
  ds.append((f-h)//2)
 k=[[0]*n for _ in range(n)]
 for i,j in combinations(range(n),2): k[i][j]=k[j][i]=ds[1<<i]+ds[1<<j]-ds[(1<<i)|(1<<j)]
 W=sum(k[i][j] for i,j in combinations(range(n),2)); deg=list(map(sum,k)); b=[2*ds[1<<i]-deg[i] for i in range(n)]; D=sum(b[i]*deg[i] for i in range(n))
 require(sum(b)==0)
 for t in range(1,n):
  rank=[ds[m] for m in range(1<<n) if m.bit_count()==t]
  other=[ds[m] for m in range(1<<n) if m.bit_count()==n-t]
  require(sum(rank)==comb(n-2,t-1)*W==sum(other))
  require((n-2)*(sum(x*x for x in rank)-sum(x*x for x in other))==(n-2*t)*comb(n-2,t-1)*D)
 return h,W,D,ds

for n in (3,4):
 count=1<<(n*(n-1)//2); nonzero=[]
 for c in range(count):
  h,w,d,ds=check(n,c)
  if d: nonzero.append((c,h,w,d))
 print('n',n,'checked',count,'nonzero D example',nonzero[:1])
for c in (8,759):
 h,w,d,ds=check(5,c); print('n=5 code',c,'H W D',h,w,d,'rank1',Counter(ds[m] for m in range(32) if m.bit_count()==1),'rank4',Counter(ds[m] for m in range(32) if m.bit_count()==4))
print('PASS')

# Independently check the general ordered-arc moment formula for the published
# code-8 carrier (positive arcs only), including repeated arcs and conflicts.
from itertools import product
c=[[0,0,0,0,4], [3,0,0,0,0], [2,5,0,0,0], [2,5,5,0,0], [0,2,2,3,0]]
n=5; edges=[(i,j,c[i][j]) for i in range(n) for j in range(n) if c[i][j]]
_,_,_,ds=check(5,8)
require(all(ds[s]==sum(w for i,j,w in edges if s>>i&1 and not(s>>j&1)) for s in range(1<<n)))
for r in range(n+1):
 for m in range(5):
  rhs=0
  for es in product(edges,repeat=m):
   tails={e[0] for e in es}; heads={e[1] for e in es}
   if tails & heads: continue
   a,b=len(tails),len(heads); t=r-a
   if not 0<=t<=n-a-b: continue
   weight=1
   for e in es: weight*=e[2]
   rhs+=weight*comb(n-a-b,t)
  require(rhs==sum(ds[s]**m for s in range(1<<n) if s.bit_count()==r))
print('ordered-arc moments: code 8, all ranks, m=0..4 PASS')
