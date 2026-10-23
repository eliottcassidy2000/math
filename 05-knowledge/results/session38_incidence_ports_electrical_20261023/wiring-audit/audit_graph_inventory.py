"""Independent recursive graph inventory and all endpoint-capacity comparison."""
import collections,csv,gzip,json,pathlib,sys
sys.setrecursionlimit(10000)
P=pathlib.Path(__file__).resolve().parent;W=P.parent/'wiring'
if not W.exists():W=P.parent/'session38-wiring'
def need(q,m):
 if not q:raise ValueError(m)
N=3528;groups={k:collections.Counter() for k in ('all','black','blue','paired_black','fundamental')}
for r in csv.DictReader((P/'independent_undirected_edges.tsv').open(),delimiter='\t'):
 u,v,col,w=map(int,(r['source_node'],r['target_node'],r['blue'],r['literal_count']));groups['all'][u,v]+=w;groups['blue' if col else 'black'][u,v]+=w
 groups['fundamental'][u,v]+=w if col else w//2
 if not col:groups['paired_black'][u,v]+=w//2
out={}
for name,E in groups.items():
 adj=[set() for _ in range(N)];deg=[0]*N
 for (u,v),w in E.items():
  deg[u]+=w;deg[v]+=w
  if u!=v:adj[u].add(v);adj[v].add(u)
 disc=[-1]*N;low=[-1]*N;bridges=[];components=[];clock=0
 def dfs(v,parent,component):
  global clock
  disc[v]=low[v]=clock;clock+=1;component.append(v)
  for w in adj[v]:
   if disc[w]<0:
    dfs(w,v,component);low[v]=min(low[v],low[w])
    if low[w]>disc[v]:bridges.append(tuple(sorted((v,w))))
   elif w!=parent:low[v]=min(low[v],disc[w])
 for v in range(N):
  if disc[v]<0:
   comp=[];dfs(v,-1,comp);components.append(sorted(comp))
 bridges.sort();support=sum(u!=v for u,v in E)
 out[name]={'vertices':N,'components':len(components),'component_sizes':sorted(map(len,components),reverse=True),'isolated_vertices':[c[0] for c in components if len(c)==1],'support_leaves':[v for v in range(N) if len(adj[v])==1],'multigraph_leaves':[v for v in range(N) if deg[v]==1],'nonloop_support_edges':support,'nonloop_lines':sum(w for (u,v),w in E.items() if u!=v),'loop_support_nodes':sum(u==v for u,v in E),'loop_lines':sum(w for (u,v),w in E.items() if u==v),'support_cycle_rank':support-N+len(components),'multigraph_cycle_rank':sum(E.values())-N+len(components),'support_bridges':[{'source':u,'target':v,'multiplicity':E[u,v]} for u,v in bridges],'multigraph_bridges':[list(e) for e in bridges if E[e]==1]}
 if name=='all':full_degree=deg
receipt=json.loads((W/'receipt.json').read_text())
for k,v in out.items():need(v==receipt['graph_inventories'][k],'entire inventory '+k)
# Every literal port contributes once to degree; loops contribute twice.
fiber=[0]*N;bluefiber=[0]*N;HC=[0]*N
for r in csv.DictReader((P/'independent_classes.tsv').open(),delimiter='\t'):
 node=int(r['node_id']);fiber[node]+=int(r['ordinary_fiber']);bluefiber[node]+=int(r['blue_fiber']);HC[node]+=int(r['closed_fiber'])
need(full_degree==fiber,'all node degrees equal literal port fibers')
closed=collections.Counter();opened=collections.Counter()
for r in csv.DictReader((P/'independent_ports.tsv').open(),delimiter='\t'):
 col=int(r['blue']);count=int(r['count'])
 closed[int(r['source_cid']),int(r['start_score']),int(r['end_score']),col]+=count
 opened[int(r['target_cid']),int(r['target_start_score']),int(r['target_end_score']),col]+=count
with gzip.open(W/'endpoint_capacities.tsv.gz','rt') as f:
 seen=set()
 for r in csv.DictReader(f,delimiter='\t'):
  key=(int(r['class_id']),int(r['start_score']),int(r['end_score']),int(r['blue']=='True'))
  need(key not in seen,'duplicate capacity');seen.add(key)
  need((closed[key],opened[key])==(int(r['closed']),int(r['open'])),'every capacity entry')
need(seen==set(closed)|set(opened),'complete endpoint capacity support')
R={'status':'PASS','endpoint_capacity_entries_compared':len(seen),'node_degree_entries_compared':len(fiber),'graph_inventories':out}
(P/'independent_graph_receipt.json').write_text(json.dumps(R,indent=2)+'\n')
print(json.dumps({'status':'PASS','endpoint_capacity_entries_compared':len(seen),'graph_summary':{k:{x:v[x] for x in ['components','support_bridges','multigraph_bridges','support_leaves','multigraph_leaves']} for k,v in out.items()}},indent=2))
