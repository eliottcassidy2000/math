"""Independent verification of mirrored finite transport counterfeits and E33 minimum.
Uses only independently generated class/port ledgers and explicit witness data.
"""
import collections,csv,itertools,json,pathlib,hashlib
P=pathlib.Path(__file__).resolve().parent;R=P.parent/'theory'
if not R.exists():R=P.parent/'session38-root'
witness_path=R/'residual_counterfeits.json';r=json.loads(witness_path.read_text())
C=[{k:int(v)for k,v in row.items()}for row in csv.DictReader((P/'independent_classes.tsv').open(),delimiter='\t')]
cache={};checks=0
def need(q,m):
 global checks
 checks+=1
 if not q:raise ValueError(m)
def sig(cid,menu):
 key=(cid,menu)
 if key in cache:return cache[key]
 out=[set()for _ in range(8)];w=C[cid]['canon']
 for k,(u,v) in enumerate(itertools.combinations(range(8),2)):
  if w>>(27-k)&1:out[u].add(v)
  else:out[v].add(u)
 inc=[set(range(8))-{u}-out[u]for u in range(8)];e=[[0]*8 for _ in range(8)]
 for u in range(8):
  for v in out[u]:e[u][v]=e[v][u]=len(out[v]&inc[u])
 c=[sum(e[u][v]for v in out[u])for u in range(8)]
 cache[cid,'c3']=sum(c)//3
 if menu=='J4':
  v4=[0]*8
  for B in itertools.combinations(range(8),4):
   if sorted(len(out[u]&set(B))for u in B)==[1,1,2,2]:
    for u in B:v4[u]+=1
  tag=[(len(out[u]),c[u],v4[u])for u in range(8)]
  z=tuple(sorted((tag[u],tuple(sorted(tag[v]for v in out[u])))for u in range(8)))
 else:
  pm=2 if menu=='E_neighbor2' else 3 if menu=='E_neighbor3' else None;qm=2 if menu=='E_edge2' else None
  z=tuple(sorted(tuple(sorted((c[v]if pm is None else c[v]%pm,e[u][v]if qm is None else e[u][v]%qm)for v in out[u]))for u in range(8)))
 cache[key]=z;return z
def mirror(k):
 i,j,a,b,col=k;return C[i]['converse_cid'],C[j]['converse_cid'],7-b,7-a,col
allw=[];needed=set()
for menu,record in r['results'].items():
 for mode in ['ordinary','symmetric']:
  if mode not in record:continue
  w=record[mode];edits={tuple(e[x]for x in ['source_class','target_class','start_score','end_score','blue']):e for e in w['edits']}
  need(len(edits)==len(w['edits']),'unique edit cells');needed.update(edits);allw.append((menu,mode,w,edits))
native={}
for row in csv.DictReader((P/'independent_ports.tsv').open(),delimiter='\t'):
 k=tuple(int(row[x])for x in ['source_cid','target_cid','start_score','end_score','blue'])
 if k in needed:native[k]=int(row['count'])
raw=(P/'independent_class_of_mask.u16le').read_bytes()
def cid(m):return raw[2*m]+256*raw[2*m+1]
tiles=[(b,a)for b in range(6)for a in range(b+2,8)];idx={x:i for i,x in enumerate(tiles)}
def key_of_mask(m):
 adj=[1<<(i-1)if i else 0 for i in range(8)]
 for k,(b,a)in enumerate(tiles):
  if m>>k&1:adj[b]|=1<<a
  else:adj[a]|=1<<b
 ref=sum(((m>>k)&1)<<idx[7-a,7-b]for k,(b,a)in enumerate(tiles))
 return (cid(m),cid(m^((1<<21)-1)),adj[7].bit_count(),adj[0].bit_count(),int(ref==m))
O=[]
for menu,mode,w,edits in allw:
 rowdelta=collections.Counter();coldelta=collections.Counter();merged=collections.Counter();undirected=collections.Counter();delta={k:e['delta']for k,e in edits.items()}
 for k,e in edits.items():
  i,j,a,b,col=k;change=e['delta'];before=native.get(k,0)
  need(before==e['before']and before+change==e['after']and e['after']>=0,'actual independent native cell and bound')
  need(col==0,'blue unchanged');rowdelta[i,a,b,col,sig(j,menu)]+=change;coldelta[j,8-a,6-b,col]+=change
  u,v=C[i]['node_id'],C[j]['node_id'];merged[u,v,col]+=change;undirected[min(u,v),max(u,v),col]+=change
  need(C[i]['canon']==e['source_canon']and C[j]['canon']==e['target_canon'],'canonical word pins')
  need(all(key_of_mask(m)==k for m in e['native_masks']),'independent literal masks')
  if change<0:need(len(set(e['native_masks']))>=-change,'distinct removed native ports')
  sig(i,menu);need(cache[j,'c3']-cache[i,'c3']==a-b-1,'triangle law')
 need(not any(rowdelta.values())and not any(coldelta.values()),'all retained source and target capacities')
 merged={k:v for k,v in merged.items()if v};undirected={k:v for k,v in undirected.items()if v}
 need(merged and undirected,'both projections changed')
 need([[*k,v]for k,v in sorted(merged.items())]==w['directed_merged_delta'],'exact directed projected changes')
 need([[*k,v]for k,v in sorted(undirected.items())]==w['undirected_merged_delta'],'exact undirected projected changes')
 if mode=='symmetric':
  need(all(delta.get(mirror(k),0)==v for k,v in delta.items()),'all reflection equalities')
  need(all(v%2==0 for v in merged.values()),'all literal black parity changes even')
  need(all(k!=mirror(k)or e['after']%2==0 for k,e in edits.items()),'fixed black cells even')
 O.append({'menu':menu,'mode':mode,'changed_cells':len(edits),'projection_entries':len(undirected),'status':'PASS'})
scan=json.loads((P/'converse_menu_scan_receipt.json').read_text());s33=scan['menus']['3,3'];need(s33['singleton_or_converse']and s33['direct_signatures']==6876,'E33 symmetric uniqueness certificate')
need(any(v['menu']=='E_edge2'and v['mode']=='symmetric'for v in O),'q<=2 obstruction')
need(any(v['menu']=='E_neighbor2'and v['mode']=='symmetric'for v in O),'p<=2 obstruction')
Rout={'status':'PASS','witnesses':O,'E33_direct_class_signatures':6876,'E33_reflection_compatible_table_unique':True,'minimal_positive_modular_wiring_alphabet':9,'minimum_scope':'The complete n8 endpoint/color/transported-signature model with full converse-reflection equalities and literal black parity. Not arbitrary encodings and not ordinary class injectivity.','lower_bound':'Reflected head-parity/exact-edge witness survives every q and p=1,2; reflected exact-head/edge-parity witness survives every p and q=1,2. Every p*q<9 has p<=2 or q<=2.','witness_input_sha256':hashlib.sha256(witness_path.read_bytes()).hexdigest(),'explicit_checks':checks}
(P/'independent_reflected_receipt.json').write_text(json.dumps(Rout,indent=2)+'\n');print(json.dumps(Rout,indent=2))
