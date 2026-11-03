from itertools import permutations
from math import prod

def f(n,mask):
 pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
 B={p for z,p in enumerate(pairs) if mask>>z&1}
 S=sorted(set(v for e in B for v in e)); U=set(range(n))-set(S)
 if not S:return 1
 def arc(x,y):return ((min(x,y),max(x,y)) in B)==(x>y)
 total=0
 for p in permutations(S):
  m=len(p)
  allowed={v: {j for j in range(m+1) if (j==0 or p[j-1]<v) and (j==m or v<p[j])} for v in U}
  bad=[j for j in range(1,m) if not arc(p[j-1],p[j])]
  for a in range(1<<len(bad)):
   closed={bad[z] for z in range(len(bad)) if a>>z&1}
   total+=(-1 if len(closed)%2 else 1)*prod(len(allowed[v]-closed) for v in U)
 return total

def brute(n,mask):
 pairs=[(i,j) for i in range(n) for j in range(i+1,n)]; B={p for z,p in enumerate(pairs) if mask>>z&1}
 def arc(x,y):return ((min(x,y),max(x,y)) in B)==(x>y)
 return sum(all(arc(x,y) for x,y in zip(p,p[1:])) for p in permutations(range(n)))
for n in range(1,6):
 for mask in range(1<<(n*(n-1)//2)):
  a,b=f(n,mask),brute(n,mask)
  if a!=b: raise AssertionError((n,mask,a,b))
 print(n,1<<(n*(n-1)//2),'passed',flush=True)
