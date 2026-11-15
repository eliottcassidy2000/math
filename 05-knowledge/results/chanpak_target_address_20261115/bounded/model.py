from itertools import combinations
import math

def host(n,c):
 a=[0]*n
 for k,(x,y) in enumerate(combinations(range(n),2)):
  if (c>>k)&1: a[x]|=1<<y
  else: a[y]|=1<<x
 return a

def prof(a):
 n=len(a);full=(1<<n)-1
 es=[[0]*n for _ in range(full+1)]; ss=[[0]*n for _ in range(full+1)]
 for k in range(n):es[1<<k][k]=ss[1<<k][k]=1
 for t in range(1,full):
  for k in range(n):
   bit=1<<k
   if t&bit:continue
   es[t|bit][k]=sum(es[t][i] for i in range(n) if a[i]&bit)
   ss[t|bit][k]=sum(ss[t][i] for i in range(n) if (a[k]>>i)&1)
 f=[1]+[sum(v) for v in ss[1:]]
 v=[[sum(es[t][i]*ss[full^t][j] for t in range(1,full)) for j in range(n)] for i in range(n)]
 s=[sum(f[t]*ss[full^t][i] for t in range(full)) for i in range(n)]
 return f[-1],sum(x*y for x,y in zip(f,f[::-1])),s,v

