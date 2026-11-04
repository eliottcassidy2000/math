"""Independent exact controls for fixed-quotient transitive blowup obstruction."""
from functools import lru_cache
from itertools import product
import json

def require(v,msg):
 if not v: raise RuntimeError(msg)
def hp(rows):
 n=len(rows)
 @lru_cache(None)
 def f(S,v):
  if S==1<<v:return 1
  return sum(f(S^(1<<v),u) for u in range(n) if S>>u&1 and rows[u]>>v&1)
 return sum(f((1<<n)-1,v) for v in range(n)) if n else 1

def strong(r):
 n=len(r)
 for root in range(n):
  S=1<<root
  while True:
   T=S
   for u in range(n):
    if S>>u&1:T|=r[u]
   if T==S:break
   S=T
  if S!=(1<<n)-1:return False
 return True

def blow(q,s):
 block=[i for i,a in enumerate(s) for _ in range(a)];r=[0]*len(block)
 for u,i in enumerate(block):
  for v,j in enumerate(block):
   if (i==j and u<v) or (i!=j and q[i]>>j&1):r[u]|=1<<v
 return r

tri=[2,4,1]
triangle=[]
for s in range(1,12):
 h=hp(blow(tri,[s,1,1]));require(h==2**s+1,'triangle formula');triangle.append([s,h])
checks=0;strong_count=0
for k in range(3,5):
 pairs=[(i,j) for i in range(k) for j in range(i+1,k)]
 for mask in range(1<<len(pairs)):
  q=[0]*k
  for a,(u,v) in enumerate(pairs):
   if mask>>a&1:q[u]|=1<<v
   else:q[v]|=1<<u
  if not strong(q):continue
  strong_count+=1
  for i in range(k):
   require(any(q[i]>>j&1 and q[j]>>l&1 and q[l]>>i&1 for j in range(k) for l in range(k)),'vertex triangle')
  for sizes in product(range(1,4),repeat=k):
   if sum(sizes)>10:continue
   r=blow(q,sizes);h=hp(r)
   require(strong(r),'strong blowup')
   require(all(h>=2**s+1 for s in sizes),'module lower bound')
   checks+=1
# hostile boundary: a nonstrong quotient has arbitrary modules yet H=1.
require(hp(blow([2,0],[7,7]))==1,'strongness boundary')
# Count-level depth uniqueness and derived scale bound.
for n in range(1,13):
 vals={5}
 for _ in range(n):vals={3*h+d for h in vals for d in (2,4)}
 require(len(vals)==2**n,'digit uniqueness')
 require(max(vals)==3**n*7-2,'max depth value')
receipt={'status':'passed','triangle_counts':triangle,'strong_quotients':strong_count,'blowup_controls':checks,'branch_depths_checked':12,'hostile_control':'TT2 with module lengths 7,7 has H=1; strong quotient essential'}
print(json.dumps(receipt,indent=2))
