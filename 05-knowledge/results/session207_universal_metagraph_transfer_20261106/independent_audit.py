import itertools,json
from functools import lru_cache

def graph(n,c):
 a=[[False]*n for _ in range(n)];k=0
 for i in range(n):
  for j in range(i+1,n):a[i][j]=bool(c>>k&1);a[j][i]=not a[i][j];k+=1
 return a

def hp(a):
 n=len(a)
 if not n:return 1
 @lru_cache(None)
 def end(s,v):
  if s==1<<v:return 1
  return sum(end(s^(1<<v),u) for u in range(n) if u!=v and s>>u&1 and a[u][v])
 return sum(end((1<<n)-1,v) for v in range(n))

def tail(a,q,t):
 n=len(a);b=[[False]*(n+t) for _ in range(n+t)]
 for i in range(t):
  for j in range(i+1,t):b[i][j]=j==i+1;b[j][i]=not b[i][j]
 for i in range(n):
  for j in range(n):b[t+i][t+j]=a[i][j]
  for j in range(t):b[t+i][j]=True
 if t:b[t-1][t+q]=True;b[t+q][t-1]=False
 return b
cases=0
for n in range(1,4):
 for code in range(1<<(n*(n-1)//2)):
  a=graph(n,code)
  for q in range(n):
   weights=[0,0]
   others=[v for v in range(n) if v!=q]
   for length in range(n):
    for suffix in itertools.permutations(others,length):
     p=(q,)+suffix
     if not all(a[u][v] for u,v in zip(p,p[1:])):continue
     left=[v for v in range(n) if v not in p]
     weights[len(p)%2]+=hp([[a[u][v] for v in left]for u in left])
   actual=[hp(tail(a,q,t))for t in range(7)]
   want=[hp(a),hp(a)+2*weights[0],hp(a)+2*sum(weights)]
   for t in range(3,7):want.append(sum(want[-3:]))
   assert actual==want,(n,code,q,actual,want)
   cases+=len(actual)
print(json.dumps({'status':'PASS','independent_pointed_tail_counts':cases,'universe':'all labelled cores n1..3, every root, tail lengths0..6'}))
