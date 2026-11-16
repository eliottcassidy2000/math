from collections import Counter
from math import isqrt

N=2025
old={(u,v) for u in range(1,N) for s in (576,2025,2601)
     if u<(v:=s-u)<N}
new={(u,1024-u) for u in range(1,512)}
classes=lambda u:min(u%9,(-u)%9)

adj=[set() for _ in range(N)]
for u,v in old:
    adj[u].add(v);adj[v].add(u)
seen=set()
order=[]
for i in range(1,N):
    if i in seen: continue
    todos=[i]; comp={i}
    while todos:
        for j in adj[todos.pop()]-comp:
            comp.add(j);todos.append(j)
    order.append((classes(i),len(comp),min(comp),max(comp),
                  sorted(x for x in comp if len(adj[x])==1)))
    seen.update(comp)
assert sorted(order)==[(0,224,9,2016,[288,576])]+[
    (r,450,r,2025-r,[]) for r in range(1,5)]
assert len(old)==2023
assert Counter(tuple(sorted((classes(u),classes(v)))) for u,v in new)=={
    (0,2):113,(1,1):56,(1,3):114,(2,4):114,(3,4):114}
assert sum(u+v==2025 and v<=1023 for u,v in old)==11
assert adj[736]=={1289,1865}
assert adj[448]=={128,1577}
assert adj[896]=={1129,1705}

added={(288,241),(335,194),(382,147),(429,100)}
assert all(u+v==529 and isqrt(u+v)**2==u+v for u,v in added)
for u,v in [(241,335),(194,382),(147,429),(100,476)]:
    adj[u].remove(v);adj[v].remove(u)
for u,v in added:
    adj[u].add(v);adj[v].add(u)
path=[576]
prev=None
while nxt:=adj[path[-1]]-{prev}:
    assert len(nxt)==1
    prev=path[-1];path.append(nxt.pop())
assert len(path)==2024 and len(set(path))==2024
print("Five pieces:",sorted(order))
print("Square 529 four-edge Hamilton repair passed; endpoints",path[0],path[-1])
