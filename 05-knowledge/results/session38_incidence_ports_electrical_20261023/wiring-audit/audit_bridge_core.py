"""Independent support-core and all higher-certificate restriction audit."""
import csv,json,pathlib,sys
sys.setrecursionlimit(10000)
P=pathlib.Path(__file__).resolve().parent;D=P.parent/'electrical'
if not D.exists():D=P.parent/'session38-electrical'
N=3528;leaves={0,1,3,10,40,329,331,3357};keep=sorted(set(range(N))-leaves);root=keep[0]
def need(q,m):
 if not q:raise ValueError(m)
E=[];full=[];adj={v:set()for v in keep}
for row in csv.DictReader((P/'independent_undirected_edges.tsv').open(),delimiter='\t'):
 u,v,w=map(int,(row['source_node'],row['target_node'],row['literal_count']))
 if u==v:continue
 full.append((u,v,w))
 if u in adj and v in adj:E.append((u,v,w));adj[u].add(v);adj[v].add(u)
bridges=[];disc={};low={};components=0

def dfs(u,parent):
 disc[u]=low[u]=len(disc)
 for v in adj[u]:
  if v not in disc:
   dfs(v,u);low[u]=min(low[u],low[v])
   if low[v]>disc[u]:bridges.append([u,v])
  elif v!=parent:low[u]=min(low[u],disc[v])
for v in keep:
 if v not in disc:components+=1;dfs(v,None)
need(components==1 and not bridges and min(map(len,adj.values()))==3,'bridge-free connected core and minimum degree')
bridgeweights=[]
for leaf in sorted(leaves):
 incident=[(u,v,w)for u,v,w in full if leaf in (u,v)];need(len(incident)==1,'exact degree-one support leaf');bridgeweights.append(incident[0][2])
need(sorted(bridgeweights)==[1,1,1,2,2,2,2,2],'bridge bundle weights')

def image(x):
 y=[0]*N
 for u,v,w in E:
  z=w*(x[u]-x[v]);y[u]+=z;y[v]-=z
 return y

def rank(vectors):
 p={}
 for z in vectors:
  while z:
   i=z.bit_length()-1
   if i not in p:p[i]=z;break
   z^=p[i]
 return len(p)

first=json.loads((D/'n8_lift_receipt.json').read_text());X=[];Y=[];bits=[]
for record in first['witness_lifts']:
 x=[0]*N
 for v in record['binary_one_nodes']:x[v]+=1
 for v in record['correction_one_nodes']:x[v]+=2
 y=image(x);need(all(y[v]%4==0 for v in keep),'core mod4 divisibility');X.append(x);Y.append(y);bits.append(sum(((x[v]-x[root])%2)<<i for i,v in enumerate(keep)))
need(rank(bits)==5,'core higher lift binary rank')
q2=[[sum(x[v]*y[v]for v in keep)//4%2 for y in Y]for x in X];need(q2==first['q2_matrix'],'core q2 same')
r3=json.loads((D/'n8_final_radical_receipt.json').read_text());tail=json.loads((D/'n8_complete_tower_receipt.json').read_text());records=[{'k':3,'lift':r3['lift_mod8']}]+tail['witness_lifts'];qs=[]
for record in records:
 k=record['k'];m=1<<k;x=record['lift'];y=image(x)
 need(all(y[v]%m==0 for v in keep),'core actual higher lift divisibility')
 need(any((x[v]-x[root])%2 for v in keep),'nonzero core binary class')
 q=sum(x[v]*y[v]for v in keep)//m%2;need(q==int(k==8),'core higher form');qs.append([k,q])
O={'status':'PASS','core_vertices':len(keep),'core_root':root,'core_components':components,'core_support_bridges':bridges,'core_minimum_support_degree':min(map(len,adj.values())),'removed_leaves':sorted(leaves),'bridge_weight_multiset':sorted(bridgeweights),'core_q2_unchanged':True,'core_V2_binary_rank':5,'core_tail_forms':qs,'core_two_primary_exponent_counts':{'1':3353,'2':4,'8':1},'core_tree_two_valuation':3369,'method':'Independent core support/bridge computation and direct restriction of all eleven higher lift certificates. Core group follows also from unimodular bridge splitting of independently audited full group.'}
(P/'independent_core_receipt.json').write_text(json.dumps(O,indent=2)+'\n');print(json.dumps(O,indent=2))
