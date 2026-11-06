import json
from functools import lru_cache
from pathlib import Path
rows=json.loads((Path(__file__).resolve().parents[1] / 'RESULT.json').read_text())['rows']

def hp(a):
 n=len(a)
 @lru_cache(None)
 def f(s,v):
  if s==1<<v:return 1
  return sum(f(s^(1<<v),u) for u in range(n) if u!=v and (s>>u&1) and a[u][v])
 return sum(f((1<<n)-1,v) for v in range(n))
for n,m,k,h,g,gg in rows:
 base=[[False]*n for _ in range(n)];z=0
 for i in range(n):
  for j in range(i+1,n):
   b=True if j==i+1 else not(m>>z&1)
   if j!=i+1:z+=1
   base[i][j]=b;base[j][i]=not b
 N=n+2*k;a=[[i<j for j in range(N)]for i in range(N)]
 for i in range(n):
  for j in range(n):a[k+i][k+j]=base[i][j]
 comp=[[a[i][j] if abs(i-j)==1 else not a[i][j] if i!=j else False for j in range(N)]for i in range(N)]
 assert hp(a)==h and hp(comp)==gg
 # Every retained consecutive arc and the closure is present.
 assert all(comp[i][i+1] for i in range(N-1)) and comp[-1][0]
 if k==2:assert gg>=g+2
 if k==3:assert gg>=4*g+1
print('PASS all36 rows by independent memoized endpoint recurrence; spanning cycles and both growth inequalities')
