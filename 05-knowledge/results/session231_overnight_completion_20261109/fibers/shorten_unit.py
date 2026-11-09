from pathlib import Path
from collections import deque
import json
root=Path(__file__).parent;w=json.loads((root/'UNIT_LOOP_WITNESS.json').read_text());h,q,A,G,T,S,E=w['core'];start=(A,G,A%2)
comp=next(c for c in json.loads((root/f'SAT_SCC_{h}_{q}.json').read_text()) if list(start) in c);nodes=sorted(map(tuple,comp));ids={v:i for i,v in enumerate(nodes)};r=ids[start];graph=[]
for a,g,p in nodes:
 x=g%2;edges=[]
 for y in (0,1):
  aa=2*a+h*y-q*x;gg=(g+q*y-h*x)//2;dest=(aa,gg,y)
  if dest not in ids:continue
  wt=a*y+h*p*y-x*gg-h*x*(gg%2);edges.append((ids[dest],wt,x,y))
 graph.append(edges)
pot={r:0};todo=deque([r])
while todo:
 v=todo.popleft()
 for z,wt,x,y in graph[v]:
  if z not in pot:pot[z]=pot[v]+wt;todo.append(z)
for i,edges in enumerate(graph):graph[i]=[(j,(wt+pot[i]-pot[j])//2,x,y) for j,wt,x,y in edges]
bound=500;origin=(r,0);pred={origin:None};todo=deque([origin]);target=None
while todo and target is None:
 v,g=todo.popleft()
 for z,wt,x,y in graph[v]:
  gg=g+wt;dest=(z,gg)
  if abs(gg)>bound or dest in pred:continue
  pred[dest]=((v,g),x,y)
  if dest==(r,1):target=dest;break
  todo.append(dest)
if target is None:print('NO FIND',len(pred));raise SystemExit
path=[];cur=target
while pred[cur] is not None:
 prev,x,y=pred[cur];path.append((x,y));cur=prev
path=path[::-1];out=dict(w);out['X']=''.join(str(x) for x,y in path);out['Y']=''.join(str(y) for x,y in path);out['search']={'gain_bound_half':bound,'visited':len(pred),'claim':'shortest only within bounded gauge-gain search'}
(root/'SHORT_UNIT_LOOP_WITNESS.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(dict(length=len(path),visited=len(pred),last_x=path[-1][0],first_y=path[0][1])))
