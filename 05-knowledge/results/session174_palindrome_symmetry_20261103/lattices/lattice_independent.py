from itertools import product,combinations
import json
checks=0
def require(x):
 global checks
 checks+=1
 if not x:raise RuntimeError('independent audit failure')
ambient=[]
for r in product((-1,1),repeat=8):
 if sum(v==1 for v in r)%2==0:ambient.append(r)
for i,j in combinations(range(8),2):
 for a,b in product((-2,2),repeat=2):
  r=[0]*8;r[i]=a;r[j]=b;ambient.append(tuple(r))
R={r for r in ambient if sum(r)==0};H={r for r in R if all(abs(x)==1 for x in r)}
require(len(ambient)==240);require(len(R)==126);require(len(H)==70);new=set(H)
for r in R:
 for x in R:
  d=sum(a*b for a,b in zip(r,x));require(d%4==0)
  y=tuple(a-d//4*b for a,b in zip(x,r));require(y in R)
  if r in H and x in H:new.add(y)
require(new==R)
for n in (2,4,6):
 edges=list(combinations(range(n),2));seen=set()
 for mask in range(1<<len(edges)):
  scores=[0]*n
  for k,(u,v) in enumerate(edges):scores[u if mask>>k&1 else v]+=1
  if set(scores)<={n//2-1,n//2}:seen.add(tuple(2*d-(n-1) for d in scores))
 expected={r for r in product((-1,1),repeat=n) if sum(r)==0};require(seen==expected)
print(json.dumps({'status':'PASS','explicit_checks':checks,'near_regular_orders':[2,4,6],'root_count':len(R),'reflection_pairs':len(R)**2},sort_keys=True))
