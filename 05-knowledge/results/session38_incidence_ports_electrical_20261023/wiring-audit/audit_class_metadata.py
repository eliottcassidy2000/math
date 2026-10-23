"""Independently verify every additional class/node metadata field and current."""
import csv,gzip,itertools,json,pathlib,collections
P=pathlib.Path(__file__).resolve().parent;W=P.parent/'wiring'
if not W.exists():W=P.parent/'session38-wiring'
def need(q,m):
 if not q:raise ValueError(m)
C=[{k:int(v)for k,v in r.items()}for r in csv.DictReader((P/'independent_classes.tsv').open(),delimiter='\t')]
I={int(r['canon']):int(r['class'])for r in csv.DictReader((P/'classes_input_n8.tsv').open(),delimiter='\t')}
meta=[]
for r in C:
 w=r['canon'];A=[[0]*8 for _ in range(8)]
 for k,(i,j)in enumerate(itertools.combinations(range(8),2)):A[i][j]=w>>(27-k)&1;A[j][i]=1-A[i][j]
 cyc=sum(A[i][j]==A[j][k]==A[k][i]for i,j,k in itertools.combinations(range(8),3));scores=','.join(map(str,sorted(map(sum,A))))
 count=r['closed_fiber']*r['Aut'];need(count%8==0,'Hamiltonian cycles integrality')
 meta.append({'c3':cyc,'scores':scores,'HC':count//8,'input_class_id':I[w]})
with gzip.open(W/'classes.tsv.gz','rt')as f:
 for p,c,m in zip(csv.DictReader(f,delimiter='\t'),C,meta):
  for key in m:need(str(m[key])==p[key],'class additional field '+key)
nodes=collections.defaultdict(list)
for c in C:nodes[c['node_id']].append(c)
c3={v:meta[g[0]['cid']]['c3']for v,g in nodes.items()};current=collections.Counter()
for r in csv.DictReader((P/'independent_undirected_edges.tsv').open(),delimiter='\t'):
 u,v,w=map(int,(r['source_node'],r['target_node'],r['literal_count']));current[u]+=w*(c3[v]-c3[u]);current[v]+=w*(c3[u]-c3[v])
with gzip.open(W/'nodes.tsv.gz','rt')as f:
 for row in csv.DictReader(f,delimiter='\t'):
  n=int(row['node_id']);G=nodes[n];c=min(G,key=lambda r:r['cid']);expected={'node_id':n,'canon':c['canon'],'class_id':c['cid'],'converse_class':c['converse_cid'],'SC':int(c['cid']==c['converse_cid']),'H':c['H'],'Aut':c['Aut'],'tilings':sum(x['ordinary_fiber']for x in G),'blue_tilings':sum(x['blue_fiber']for x in G),'closed_tilings':sum(x['closed_fiber']for x in G),'c3':c3[n],'current_c3':current[n]}
  need({k:int(v)for k,v in row.items()}==expected,'every node metadata field')
need(sum(current.values())==0,'current conservation')
O={'status':'PASS','additional_class_fields':4*len(C),'all_node_fields':12*len(nodes),'method':'Direct cyclic triples and score decode; Hamiltonian-cycle counts from independently counted closed path-orbits and audited automorphism orders; full integer currents from independently reconstructed undirected edges.'}
(P/'independent_metadata_receipt.json').write_text(json.dumps(O,indent=2)+'\n');print(json.dumps(O,indent=2))
