from itertools import combinations
from collections import defaultdict
from math import factorial
def profile(ar):
 n=len(ar); full=(1<<n)-1
 dp=[[0]*n for _ in range(full+1)]
 for u in range(n):dp[1<<u][u]=1
 for s in range(1,full+1):
  for u,a in enumerate(dp[s]):
   if a:
    for v in range(n):
     if not (s>>v)&1 and ar[u][v]:dp[s|1<<v][v]+=a
 h=[sum(x) for x in dp]
 p=[[0]*(n+1) for _ in range(full+1)];p[0][0]=1
 for s in range(1,full+1):
  lo=s&-s;t=s
  while t:
   if t&lo and h[t]:
    for i in range(s.bit_count()-t.bit_count()+1):
     p[s][i+1]+=h[t]*p[s^t][i]
   t=(t-1)&s
 return [x*factorial(i) for i,x in enumerate(p[-1])]
def decode(n,c):
 ar=[[0]*n for _ in range(n)]
 for k,(i,j) in enumerate(combinations(range(n),2)):ar[i][j]=(c>>k)&1;ar[j][i]=1-ar[i][j]
 return ar
if __name__=='__main__':
 for n in range(3,7):
  seen={}
  for c in range(1<<(n*(n-1)//2)):
   p=profile(decode(n,c))
   k=tuple(p[1:3])
   if k in seen and p[3]!=seen[k][1][3]:
    print(n,c,p,seen[k]);break
   seen[k]=(c,p)
  else:continue
  break
