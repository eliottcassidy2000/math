from itertools import combinations
from collections import deque
from math import factorial
import json

def nbrs(st):
 p,s=st;k=len(p)
 for i in range(k):yield p,s^(1<<i)
 for i,j in combinations(range(k),2):
  mask=(1<<i)|(1<<j)
  yield p,s^mask
  q=list(p);q[i],q[j]=q[j],q[i];q=tuple(q)
  for bits in (0,1<<i,1<<j,mask):yield q,(s&~mask)|bits

def dist_formula(x,y):
 p,s=x;t,z=y;k=len(p);ti={v:i for i,v in enumerate(t)}
 r=[ti[v] for v in p];seen=set();c=0
 for i in range(k):
  if i not in seen:
   c+=1;j=i
   while j not in seen:seen.add(j);j=r[j]
 f=sum(p[i]==t[i] and bool((s^z)>>i&1) for i in range(k))
 return k-c+(f+1)//2

def check(k,origin):
 seen={origin:0};q=deque([origin])
 while q:
  st=q.popleft()
  for nb in nbrs(st):
   if nb not in seen:seen[nb]=seen[st]+1;q.append(nb)
 if len(seen)!=factorial(k)*2**k:raise RuntimeError('state count')
 for y,d in seen.items():
  f=dist_formula(origin,y)
  if f!=d:raise RuntimeError((k,origin,y,d,f))
 return {'k':k,'origin':[list(origin[0]),origin[1]],'states_checked':len(seen),'eccentricity':max(seen.values())}
rows=[]
for k in range(1,6):rows.append(check(k,(tuple(range(k)),0)))
rows.append(check(5,((2,0,4,1,3),21)))
print(json.dumps({'status':'PASS','distance_checks':sum(x['states_checked'] for x in rows),'rows':rows},indent=2))
