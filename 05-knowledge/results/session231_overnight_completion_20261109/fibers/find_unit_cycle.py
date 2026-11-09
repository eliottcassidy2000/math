from pathlib import Path
from collections import deque
import json,math
rootdir=Path(__file__).parent
allrows=[]
for n in range(3,10):
 for row in map(json.loads,(rootdir/f'SATURATED{n}.jsonl').read_text().splitlines()):allrows.append((n,row))
for h,q in [(37,740),(45,1440),(37,1480)]:
 for ci,comp in enumerate(json.loads((rootdir/f'SAT_SCC_{h}_{q}.json').read_text())):
  nodes={tuple(v) for v in comp};physical=[(n,r,(r[2],r[3],r[2]%2)) for n,r in allrows if r[:2]==[h,q] and (r[2],r[3],r[2]%2) in nodes]
  if not physical:continue
  graph={};reverse={v:[] for v in nodes}
  for a,g,p in nodes:
   v=(a,g,p);x=g%2;edges=[]
   for y in (0,1):
    aa=2*a+h*y-q*x;gg=(g+q*y-h*x)//2;w=(aa,gg,y)
    if w not in nodes:continue
    wt=a*y+h*p*y-x*gg-h*x*(gg%2);edges.append((w,wt,x,y));reverse[w].append((v,wt,x,y))
   graph[v]=edges
  for n,row,start in physical:
   parent={start:None};gain={start:0};dist={start:0};todo=deque([start])
   while todo:
    v=todo.popleft()
    for w,wt,x,y in graph[v]:
     if w not in parent:parent[w]=(v,x,y);gain[w]=gain[v]+wt;dist[w]=dist[v]+1;todo.append(w)
   nxt={start:None};back={start:0};db={start:0};todo=deque([start])
   while todo:
    w=todo.popleft()
    for v,wt,x,y in reverse[w]:
     if v not in nxt:nxt[v]=(w,x,y);back[v]=back[w]+wt;db[v]=db[w]+1;todo.append(v)
   cycles={}
   for v in sorted(nodes):
    for w,wt,x,y in graph[v]:
     delta=gain[v]+wt+back[w]
     if not delta:continue
     L=dist[v]+1+db[w]
     if delta in cycles and len(cycles[delta])<=L:continue
     path=[];cur=v
     while parent[cur] is not None:
      prev,px,py=parent[cur];path.append((px,py));cur=prev
     path=path[::-1]+[(x,y)];cur=w
     while nxt[cur] is not None:
      nex,px,py=nxt[cur];path.append((px,py));cur=nex
     cycles[delta]=path
   print(json.dumps(dict(h=h,q=q,core=row,n=n,cycle_gains={d:len(p) for d,p in cycles.items()})),flush=True)
   choices=[]
   for delta,path in cycles.items():
    if abs(delta)==2:choices.append((len(path),delta,path))
   for pos,pp in cycles.items():
    if pos<=0:continue
    for neg,pn in cycles.items():
     if neg>=0:continue
     for a in range(1,-neg+2):
      num=a*pos-2
      if num>=0 and num%(-neg)==0:
       b=num//(-neg);choices.append((a*len(pp)+b*len(pn),2,pp*a+pn*b));break
   if choices:
    L,delta,path=min(choices,key=lambda v:v[0]);adj=[0]*n;code=row[7];bit=0
    for i in range(n):
     for j in range(i+1,n):
      if code>>bit&1:adj[i]|=1<<j
      else:adj[j]|=1<<i
      bit+=1
    out=dict(adj=adj,mask=row[8],core=row[:7],cycle_gain=delta,X=''.join(str(x) for x,y in path),Y=''.join(str(y) for x,y in path))
    (rootdir/'UNIT_LOOP_WITNESS.json').write_text(json.dumps(out,indent=2)+'\n');print('FOUND UNIT',h,q,L,flush=True);raise SystemExit
