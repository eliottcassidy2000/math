import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from pathlib import Path
from collections import deque
import json
from core import states,compose
rootdir=Path(__file__).parent
nodes={tuple(v) for comp in json.loads((rootdir/'NONZERO_SCC.json').read_text()) for v in comp}
h,q=95,1276;graph={};reverse={v:[] for v in nodes}
for a,g,p in nodes:
 v=(a,g,p);x=g%2;edges=[]
 for y in (0,1):
  aa=2*a+h*y-q*x;gg=(g+q*y-h*x)//2;w=(aa,gg,y)
  if w not in nodes:continue
  wt=a*y+h*p*y-x*gg-h*x*(gg%2);edges.append((w,wt,x,y));reverse[w].append((v,wt,x,y))
 graph[v]=edges
start=(270,372,0);adj=[34,20,9,35,13,86,31];mask=17
parent={start:None};gain={start:0};dist={start:0};todo=deque([start])
while todo:
 v=todo.popleft()
 for w,wt,x,y in graph[v]:
  if w not in parent:parent[w]=(v,x,y);gain[w]=gain[v]+wt;dist[w]=dist[v]+1;todo.append(w)
nextnode={start:None};gainback={start:0};distback={start:0};todo=deque([start])
while todo:
 w=todo.popleft()
 for v,wt,x,y in reverse[w]:
  if v not in nextnode:nextnode[v]=(w,x,y);gainback[v]=gainback[w]+wt;distback[v]=distback[w]+1;todo.append(v)
choices=[]
for v in sorted(nodes):
 for w,wt,x,y in graph[v]:
  delta=gain[v]+wt+gainback[w]
  if delta:choices.append((dist[v]+1+distback[w],abs(delta),v,w,x,y,delta))
L,_,v,w,x,y,delta=min(choices);path=[];cur=v
while parent[cur] is not None:
 prev,px,py=parent[cur];path.append((px,py));cur=prev
path=path[::-1]+[(x,y)];cur=w
while nextnode[cur] is not None:
 nex,px,py=nextnode[cur];path.append((px,py));cur=nex
if len(path)!=L:raise ValueError('length')
core=states(adj)[mask];one=[(1,2,b,b,0,b,b) for b in (0,1)]
def fold(parts):
 s=parts[0]
 for p in parts[1:]:s=compose(s,p)
 return s
u=path[-1][0];v=path[0][1];rows=[]
for k in range(1,7):
 zs=[]
 for j in range(k+1):
  parts=([one[1]] if not u else [])+[one[u]]+[one[x] for x,y in path]*j+[core]+[one[y] for x,y in path]*(k-j)+[one[v]]
  zs.append(fold(parts))
 if len({s[:4]+s[5:] for s in zs})!=1:raise ValueError('interface')
 differences=[zs[j+1][4]-zs[j][4] for j in range(k)]
 if len(set(differences))!=1 or abs(differences[0])!=abs(delta):raise ValueError(('gain',differences,delta))
 rows.append(dict(k=k,T_differences=differences,interface=zs[0][:4]+zs[0][5:],T_values=[s[4] for s in zs]))
out=dict(adj=adj,mask=mask,core=core,start=start,cycle_length=L,cycle_gain=delta,first_y=v,last_x=u,path=path,checks=rows)
(rootdir/'PHYSICAL_CYCLE.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(dict(cycle_length=L,cycle_gain=delta,first_y=v,last_x=u,core=core,path=path,checks_k=list(range(1,7)))))
