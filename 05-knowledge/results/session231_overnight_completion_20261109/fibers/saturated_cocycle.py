import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from core import states
from collections import defaultdict,deque
import json
banks={'H5':[[8,1,3,6]],'H9':[[16,1,3,7,14]],'H15':[[24,17,3,6,12]],'H29':[[26,44,57,16,34,9]],'H45':[[48,49,3,7,12,28],[40,37,17,6,11,28]],'H95':[[64,121,35,37,13,81,28]],'H389':[[20,85,72,195,204,159,33,71]]}
from pathlib import Path
bankstates=defaultdict(list)
for n in range(3,9):
 for row in map(json.loads,Path(__file__).with_name(f'SATURATED{n}.jsonl').read_text().splitlines()):bankstates[row[0],row[1]].append(row)
for (h,q),rows in sorted(bankstates.items()):
 name=f'H{h}Q{q}';ss=[row[:7] for row in rows];starts={(s[2],s[3],b) for s in ss for b in (0,1)};todo=deque(starts);graph={};seen=set(starts);reverse=defaultdict(list)
 while todo:
  node=todo.popleft();a,g,prev=node;x=g%2;edges=[]
  for y in (0,1):
   ap=2*a+h*y-q*x;gp=(g+q*y-h*x)//2;dest=(ap,gp,y)
   if not(-h<ap<q and -h<gp<q):continue
   w=a*y+h*prev*y-x*gp-h*x*(gp%2)
   edges.append((dest,w));reverse[dest].append(node)
   if dest not in seen:seen.add(dest);todo.append(dest)
  graph[node]=edges
 co=set(starts);todo=deque(starts)
 while todo:
  for p in reverse[todo.popleft()]:
   if p not in co:co.add(p);todo.append(p)
 # Iterative Kosaraju.
 visit=set();order=[]
 for root in sorted(co):
  if root in visit:continue
  visit.add(root);stack=[(root,iter(graph[root]))]
  while stack:
   v,it=stack[-1]
   try:w,_=next(it)
   except StopIteration:order.append(v);stack.pop();continue
   if w in co and w not in visit:
    visit.add(w);stack.append((w,iter(graph[w])))
 components=[];which={}
 for root in reversed(order):
  if root in which:continue
  idx=len(components);part=[];stack=[root];which[root]=idx
  while stack:
   v=stack.pop();part.append(v)
   for w in reverse[v]:
    if w in co and w not in which:which[w]=idx;stack.append(w)
  components.append(part)
 nonzero=[];cyclic=0
 for idx,part in enumerate(components):
  if len(part)>1 or any(w==part[0] for w,_ in graph[part[0]]):cyclic+=1
  potential={part[0]:0};todo=[part[0]]
  while todo:
   v=todo.pop()
   for w,wt in graph[v]:
    if which.get(w)!=idx:continue
    if w not in potential:potential[w]=potential[v]+wt;todo.append(w)
    elif potential[w]!=potential[v]+wt:nonzero.append((len(part),sum(z in starts for z in part),v,w,potential[v]+wt-potential[w]))
 if nonzero:
  from pathlib import Path
  ids={which[tuple(row[2])] for row in nonzero}
  Path(__file__).with_name(f'SAT_SCC_{h}_{q}.json').write_text(json.dumps([sorted(components[i]) for i in sorted(ids)]))
 print(json.dumps(dict(bank=name,reachable=len(seen),useful=len(co),SCCs=len(components),cyclic=cyclic,nonzero=nonzero)),flush=True)
