import sys,json,itertools
from collections import deque
from exact_utils import enumerate_fiber

def graph(states):
 buckets={};adj=[set() for _ in states];r=len(states[0])
 # share r-2 coordinates
 for i,t in enumerate(states):
  for omitted in itertools.combinations(range(r),min(2,r)):
   key=(omitted,tuple(t[j] for j in range(r) if j not in omitted))
   for k in buckets.get(key,[]):adj[i].add(k);adj[k].add(i)
   buckets.setdefault(key,[]).append(i)
 return adj

def reach(adj,start,allowed=None):
 q=deque([start]);prev={start:None}
 while q:
  x=q.popleft()
  for y in adj[x]:
   if y not in prev and (allowed is None or y in allowed):prev[y]=x;q.append(y)
 return prev
A=[('a','b','c','e'),('a','b','e'),('b','c'),('c','a')]
states=[x for x in itertools.product(*A) if len(set(x))==4]
adj=graph(states);source=('e','a','b','c');target=('e','b','c','a');u=states.index(source);v=states.index(target)
prev=reach(adj,u);path=[v]
while path[-1]!=u:path.append(prev[path[-1]])
path.reverse();allowed={i for i,t in enumerate(states) if t[0]=='e'}
print(json.dumps({'multistar_states':states,'distance':len(path)-1,'path':[states[i] for i in path],'fixed_A1_reachable':[states[i] for i in reach(adj,u,allowed)]}),flush=True)
for m in range(3,9):
 alltrees=enumerate_fiber(m,m,{})
 ts=[tuple(t[d] for d in range(1,m+1)) for t in alltrees];adj=graph(ts)
 found=None
 for d in range(m):
  for e in sorted(set(t[d] for t in ts)):
   allowed={i for i,t in enumerate(ts) if t[d]==e};u=min(allowed);seen=reach(adj,u,allowed)
   if len(seen)<len(allowed):
    v=min(allowed-set(seen));found={'m':m,'fixed_color':d+1,'edge':e,'source':ts[u],'target':ts[v],'restricted_states':len(allowed),'source_component':len(seen)};break
  if found:break
 print(json.dumps({'m':m,'states':len(ts),'obstruction':found}),flush=True)
 if found:break
