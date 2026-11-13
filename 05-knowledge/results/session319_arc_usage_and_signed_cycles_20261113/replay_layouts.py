from itertools import permutations
from collections import Counter
from arcscan import stats
A=(54,64,122,115,98,66,1)
B=(54,32,122,115,98,64,3)
t=3; ext=tuple(i for i in range(7) if i!=t)
def layouts(R,r):
 for p in permutations(range(7+r-1)):
  if r==2 and p.index(t)>p.index(7):continue
  v=[t if i==7 else i for i in p]
  if all(R[a]>>b&1 for a,b in zip(v,v[1:])):yield tuple(v)
def inflate(R,m):
 q=ext+(t,)*m
 return [sum(1<<j for j,b in enumerate(q) if i!=j and ((i<j) if a==b else (R[a]>>b&1))) for i,a in enumerate(q)]
for R,zeros in ((A,{(2,6)}),(B,{(0,5),(2,5),(2,6)})):
 L={r:list(layouts(R,r)) for r in (1,2)}
 if (len(L[1]),len(L[2]))!=(27,9):raise ValueError(R)
 support={e for r in L for v in L[r] for e in zip(v,v[1:])}
 two={e for v in L[2] for e in zip(v,v[1:])}
 universe={(i,j)for i in range(7)for j in range(7)if R[i]>>j&1}
 if universe-support!=zeros or not {e for e in universe if t in e}<=two:raise ValueError(R)
 for m in range(2,7):
  T=inflate(R,m); H,W,*_=stats(T)
  lit={(ext[i] if i<6 else t,ext[j] if j<6 else t) for i in range(m+6) for j in range(m+6) if T[i]>>j&1 and W[i][j]==0}
  if H!=9*(2**m+1) or lit!=zeros:raise ValueError((R,m,lit))
  for i in range(1,m+1):
   for j in range(i+1,m+1):
    b=j==i+1
    if W[i+5][j+5]!=27*b+18*(2**(m-j+i-1)-b):raise ValueError((R,m,i,j))
print('one- and two-run quotient supports and exact literal weights checked')
