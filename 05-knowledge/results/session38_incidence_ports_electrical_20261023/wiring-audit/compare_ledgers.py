"""Compare independent exact audit against producer artifacts, preserving loops/colors."""
import array,collections,csv,gzip,hashlib,json,pathlib
HERE=pathlib.Path(__file__).resolve().parent
PROD=HERE.parent/'wiring'
if not PROD.exists():PROD=HERE.parent/'session38-wiring'
checks=0
def need(x,m):
 global checks
 checks+=1
 if not x: raise ValueError(m)
def read(path):
 with (gzip.open(path,'rt') if path.suffix=='.gz' else path.open()) as f:
  yield from csv.DictReader(f,delimiter='\t')
a=(HERE/'independent_class_of_mask.u16le').read_bytes();b=gzip.decompress((PROD/'class_of_mask.u16le.gz').read_bytes())
need(a==b and len(a)==4194304,'all-mask assignment mismatch')
C=list(read(HERE/'independent_classes.tsv'));P=list(read(PROD/'classes.tsv.gz'))
need(len(C)==len(P)==6880,'class cardinality')
fields={'cid':'class_id','canon':'canon','converse_cid':'converse_class','node_id':'node_id','H':'H','Aut':'Aut','ordinary_fiber':'fiber','blue_fiber':'blue_fiber','closed_fiber':'closed_fiber'}
for c,p in zip(C,P):
 for k,v in fields.items(): need(c[k]==p[v],f'class field {k}: {c["cid"]}')
E={}
for r in read(HERE/'independent_directed_edges.tsv'):
 k=(int(r['source_node']),int(r['target_node']));col=int(r['blue']);E.setdefault(k,[0,0])[1-col]=int(r['literal_count'])
need(len(E)==460699,'directed endpoint pairs')
pe={}
for r in read(PROD/'edges.tsv.gz'):
 k=(int(r['source']),int(r['target']));need(k not in pe,'producer duplicate edge row');pe[k]=[int(r['blue']),int(r['black'])]
need(E==pe,'all directed colored multiplicities mismatch')
ports={}
for r in read(HERE/'independent_ports.tsv'):
 u,v,aa,bb,col=map(int,(r['source_cid'],r['target_cid'],r['start_score'],r['end_score'],r['blue']))
 need(int(r['target_start_score'])==8-aa and int(r['target_end_score'])==6-bb,'target endpoint scores')
 k=(u,v,aa,bb,col);need(k not in ports,'independent duplicate port');ports[k]=int(r['count'])
pports={}
for r in read(PROD/'ports.tsv.gz'):
 k=tuple(int(r[x]) for x in ('source_class','target_class','start_score','end_score','blue'))
 need(k not in pports,'producer duplicate port');pports[k]=int(r['multiplicity'])
need(ports==pports,'every port entry mismatch')
# An independent incidence graph connected-components calculation, loops retained but not bridges.
adj=[set() for _ in range(3528)];badj=[set() for _ in range(3528)];blueadj=[set() for _ in range(3528)]
for (u,v),(bl,bk) in E.items():
 for ar,include in ((adj,bl+bk),(badj,bk),(blueadj,bl)):
  if include: ar[u].add(v);ar[v].add(u)
def components(ar):
 unseen=set(range(len(ar)));out=[]
 while unseen:
  s=unseen.pop();todo=[s];cc=[s]
  while todo:
   v=todo.pop();fresh=ar[v]&unseen;unseen.difference_update(fresh);todo.extend(fresh);cc.extend(fresh)
  out.append(tuple(sorted(cc)))
 return sorted(out,key=lambda x:(-len(x),x))
comp={k:components(v) for k,v in [('all',adj),('black',badj),('blue',blueadj)]}
need(len(comp['all'])==1,'merged graph connectedness')
O={'status':'PASS','mask_assignment_sha256':hashlib.sha256(a).hexdigest(),'mask_assignment_compared':len(a)//2,'class_fields_compared':len(C)*len(fields),'port_entries_compared':len(ports),'directed_endpoint_pairs_compared':len(E),'colored_directed_entries_compared':sum(bool(v) for vals in E.values() for v in vals),'components':{k:{'count':len(v),'size_histogram':dict(collections.Counter(map(len,v))),'largest':len(v[0])} for k,v in comp.items()},'explicit_checks':checks,'producer_files':{n:hashlib.sha256((PROD/n).read_bytes()).hexdigest() for n in ['classes.tsv.gz','edges.tsv.gz','ports.tsv.gz']}}
(HERE/'comparison_receipt.json').write_text(json.dumps(O,indent=2,sort_keys=True)+'\n');print(json.dumps(O,indent=2,sort_keys=True))
