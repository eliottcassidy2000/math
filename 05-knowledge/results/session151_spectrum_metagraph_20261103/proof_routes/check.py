from functools import lru_cache
from itertools import combinations

def paths(adj,vs):
 if not vs:return 1
 @lru_cache(None)
 def f(mask,last):
  if mask==1<<last:return 1
  return sum(f(mask^(1<<last),i) for i in vs if mask>>i&1 and adj[i][last])
 mask=sum(1<<i for i in vs)
 return sum(f(mask,i) for i in vs)
def F(adj,v):
 vs=[i for i in range(len(adj)) if i!=v]; k=len(vs)
 @lru_cache(None)
 def f(mask,last):
  if mask==1<<last:return 1
  return sum(f(mask^(1<<last),i) for i in vs if mask>>i&1 and adj[i][last])
 out=0
 for code in range(1<<k):
  s=[vs[i] for i in range(k) if code>>i&1];t=[vs[i] for i in range(k) if not code>>i&1]
  sm=sum(1<<i for i in s)
  first=1 if not s else sum(f(sm,i) for i in s if adj[i][v])
  out+=first*paths(adj,t)
 return out
def P(adj):
 n=len(adj)
 return sum(paths(adj,[i for i in range(n) if s>>i&1])*paths(adj,[i for i in range(n) if not s>>i&1]) for s in range(1<<n))
def mat(n,mask):
 a=[[0]*n for _ in range(n)]
 for bit,(i,j) in enumerate(combinations(range(n),2)):
  a[i][j]=(mask>>bit)&1;a[j][i]=1-a[i][j]
 return a
def combine(a,b,v,k):
 m,n=len(a),len(b);N=m+n+k
 t=[[0]*N for _ in range(N)]
 for i in range(m):
  for j in range(m):t[i][j]=a[i][j]
  for j in range(m,m+n):t[i][j]=1
 for i in range(n):
  for j in range(n):t[m+i][m+j]=b[i][j]
 for x in range(m+n,N):
  for y in range(m+n):t[y][x]=int(y==m+v);t[x][y]=int(y!=m+v)
  for y in range(m+n,x):t[y][x]=1
 return t
cnt=0
for m in range(1,4):
 for ai in range(1<<(m*(m-1)//2)):
  a=mat(m,ai);ha=paths(a,list(range(m)));pa=P(a)
  for n in range(1,4):
   for bi in range(1<<(n*(n-1)//2)):
    b=mat(n,bi);hb=paths(b,list(range(n)))
    for v in range(n):
     qb=paths(b,[i for i in range(n) if i!=v]);fb=F(b,v)
     for k in (1,2,3):
      t=combine(a,b,v,k); got=paths(t,list(range(len(t))))
      expected=ha*hb+(2-2**k)*ha*qb+(2**k-1)*pa*fb
      if got!=expected:raise Exception((m,ai,n,bi,v,k,got,expected))
      cnt+=1
print('checked',cnt,'cases')
for m in (1,2,3):
 for ai in range(1<<(m*(m-1)//2)):
  a=mat(m,ai);h=paths(a,list(range(m)));p=P(a)
  print('A',m,ai,'H',h,'P',p)
from math import isqrt
c=mat(4,4)
if paths(c,list(range(4)))!=5 or P(c)!=32: raise RuntimeError('seed counts')
roots=[F(c,v) for v in range(4)]
if roots!=[7,9,6,5]: raise RuntimeError(('root counts',roots))
print('seed root coefficients',roots)
for n in (31,313,6269,9341,197233):
 if any(n%d==0 for d in range(2,isqrt(n)+1)): raise RuntimeError(('not prime',n))
 print('trial division prime',n,'through',isqrt(n))
for v in range(4):
 for m in (1,2,3):
  z=5*5**m+roots[v]*32**m
  print('root family',v,m,z)
