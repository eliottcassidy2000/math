from itertools import combinations,permutations,product
from collections import defaultdict,Counter,deque
from pathlib import Path
import json
D=Path(__file__).parent;P=list(combinations(range(5),2));IX={p:i for i,p in enumerate(P)}
PER=list(permutations(range(5)));ROOT=[p for p in PER if p[4]==4]
def require(b,s):
 if not b:raise RuntimeError(s)
def rel(m,p):
 z=0
 for j,(i,k) in enumerate(P):
  a,b=p[i],p[k];yes=(m>>j)&1
  if a>b:a,b=b,a;yes=1-yes
  z|=yes<<IX[a,b]
 return z
def can(m,ps):return min(rel(m,p) for p in ps)
def tau(m):return rel(m^1023,(4,3,2,1,0))
spine=sum(1<<IX[i,i+1] for i in range(4));free=1023^spine
T=[m for m in range(1024) if m&spine==spine]
row={t:can(t,PER) for t in T};col={t:can(t^free,ROOT) for t in T};ri=defaultdict(list);ci=defaultdict(list)
for t in T:ri[row[t]].append(t);ci[col[t]].append(t)
require(len(ri)==12 and len(ci)==33,'alphabet')
# An ordered pair (e,f) means actual spine letters (e,tau(f)).
def endpoints(e,f):return (('S',)+tuple(sorted((row[e],row[f]))),('Q',)+tuple(sorted((col[e],col[f]))))
edges={};adj=defaultdict(list)
for e,f in product(T,repeat=2):
 u,v=endpoints(e,f);edges[e,f]=(u,v);adj[u].append(((e,f),v));adj[v].append(((e,f),u))
def alternate(e):
 forbidden=(e,e);s,t=edges[forbidden];todo=deque([s]);prev={s:None}
 while todo and t not in prev:
  u=todo.popleft()
  for eid,v in adj[u]:
   if eid==forbidden or v in prev:continue
   prev[v]=(u,eid);todo.append(v)
 if t not in prev:return None
 path=[];u=t
 while prev[u] is not None:
  u,eid=prev[u];path.append(eid)
 return list(reversed(path))
old=json.loads((D.parent/'session29-target-wiring/wiring-results.json').read_text());oldbridges=json.loads((D.parent/'session29-target-wiring/cycle-results.json').read_text())['restricted_multigraph_bridges'];result=[]
for b in oldbridges:
 w=old['witness_rows'][b['edge']];e,u=w['tiles'];require(u==tau(e),'blue diagonal')
 path=alternate(e)
 result.append({'old_edge':b['edge'],'old_source':b['source'],'tile':e,'actual_original_word':[e,u],'root':col[e],'base_source_class':row[e],'row_degree':len(ri[row[e]]),'column_degree':len(ci[col[e]]),'column_tiles':ci[col[e]],'alternate_normalized_pairs':path,'alternate_actual_words':None if path is None else [[a,tau(bb)] for a,bb in path],'H':b['H']})
allbridges=[e for e in T if alternate(e) is None]
res={'tile_count':len(T),'source_classes':len(ri),'root_classes':len(ci),'boundary_source_nodes':len({a for a,b in edges.values()}),'boundary_target_nodes':len({b for a,b in edges.values()}),'boundary_edges':len(edges),'base_incidence':[{ 'tile':t,'source_class':row[t],'root_class':col[t]} for t in T],'original_eight':result,'all_boundary_bridge_tiles':allbridges}
(D/'results.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps({k:v for k,v in res.items() if k!='base_incidence'},indent=2))
