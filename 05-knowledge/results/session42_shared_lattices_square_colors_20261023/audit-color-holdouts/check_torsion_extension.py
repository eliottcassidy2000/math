from itertools import combinations,product
from math import gcd
import json
checks=0
for n in range(1,5):
 E=list(combinations(range(n),2))
 for mask in range(1<<len(E)):
  edges=[e for i,e in enumerate(E) if mask>>i&1];adj=[[] for _ in range(n)]
  for a,b in edges:adj[a].append(b);adj[b].append(a)
  sides={};k=bip=0
  for root in range(n):
   if root in sides:continue
   k+=1;sides[root]=0;queue=[root];ok=True
   for a in queue:
    for b in adj[a]:
     if b not in sides:sides[b]=1-sides[a];queue.append(b)
     elif sides[a]==sides[b]:ok=False
   bip+=ok
  for q in (4,6):
   count=sum(all((c[a]+c[b])%q==0 for a,b in edges) for c in product(range(q),repeat=n))
   if count!=q**bip*gcd(q,2)**(k-bip):raise RuntimeError((n,mask,q))
   checks+=1
print(json.dumps({'status':'PASS','checks':checks,'cyclic_color_groups':[4,6],'all_simple_graph_orders':[1,2,3,4],'kernel_formula':'q^b * gcd(q,2)^(k-b)'},indent=2))
