from pathlib import Path
from collections import defaultdict,deque,Counter
import json
HERE=Path(__file__).resolve().parent;D=json.loads((HERE/'wiring-results.json').read_text())
rows=D['witness_rows'];targets=sorted({tuple(r['merged_key']) for r in rows});sids={'AA':0,'AB/BA':1,'BB':2};tid={k:i+3 for i,k in enumerate(targets)}
edges=[(sids[r['source']],tid[tuple(r['merged_key'])]) for r in rows];adj=defaultdict(list)
for z,(a,b) in enumerate(edges):adj[a].append((b,z));adj[b].append((a,z))
def path(a,b,exclude=None,black=False):
 seen={a:None};q=deque([a])
 while q:
  v=q.popleft()
  if v==b:
   out=[]
   while seen[v] is not None:v,e=seen[v];out.append(e)
   return out[::-1]
  for w,e in adj[v]:
   if e==exclude or (black and rows[e]['blue']):continue
   if w not in seen:seen[w]=(v,e);q.append(w)
 return None
def require(b,msg):
 if not b:raise RuntimeError(msg)
certs=[];bridges=[]
for z,(a,b) in enumerate(edges):
 p=path(a,b,z)
 if p is None:bridges.append(z)
 else:
  v=a
  for e in p:
   require(e!=z and v in edges[e],'cycle step');v=edges[e][0]^edges[e][1]^v
  require(v==b,'cycle closes');certs.append({'edge':z,'alternative_path':p})
require(len(bridges)==8 and all(rows[e]['blue'] for e in bridges),'eight blue restricted bridges')
blacknodes={v for e in range(len(edges)) if not rows[e]['blue'] for v in edges[e]};require(all(path(0,v,black=True) is not None for v in blacknodes),'connected black boundary core')
blackmult=Counter(edges[e] for e in range(len(edges)) if not rows[e]['blue']);require(all(m%2==0 for m in blackmult.values()),'black parallel pairing')
simple=Counter(edges);leaf_targets=[v for v in tid.values() if len({a for a,_ in adj[v]})==1]
result={'vertices':len(adj),'target_nodes':len(targets),'multigraph_edges':len(edges),'underlying_simple_edges':len(simple),'underlying_simple_leaf_targets':len(leaf_targets),'restricted_multigraph_bridges':[{'edge':z,'source':rows[z]['source'],'target_key':rows[z]['merged_key'],'H':rows[z]['H']} for z in bridges],'full_metagraph_certified_nonbridge_edges':len(certs),'black_core_vertices':len(blacknodes),'black_core_edges':sum(not r['blue'] for r in rows),'cycle_certificates':certs}
(HERE/'cycle-results.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='cycle_certificates'},indent=2))
