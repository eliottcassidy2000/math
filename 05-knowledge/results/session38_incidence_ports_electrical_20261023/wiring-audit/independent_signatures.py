"""Independent incidence audit by directed two-step intersections, not cyclic-triple enumeration.
Class coverage and canonicality are separately certified by independent_wiring.cpp.
"""
import collections,csv,hashlib,itertools,json,pathlib
P=pathlib.Path(__file__).resolve().parent
checks=0
def need(q,msg):
 global checks
 checks+=1
 if not q:raise ValueError(msg)
def unpack(w):
 out=[set() for _ in range(8)]
 for k,(u,v) in enumerate(itertools.combinations(range(8),2)):
  if w>>(27-k)&1:out[u].add(v)
  else:out[v].add(u)
 return out
def data(w):
 out=unpack(w);inc=[set(range(8))-{u}-out[u] for u in range(8)]
 e=[[0]*8 for _ in range(8)]
 for u in range(8):
  for v in out[u]:e[u][v]=e[v][u]=len(out[v]&inc[u])
 c=[sum(e[u][v] for v in out[u]) for u in range(8)]
 need(all(sum(e[u][v] for v in inc[u])==c[u] for u in range(8)),'incoming cyclic incidence identity')
 # Cross-check S4 edge-overlap/vertex identity against direct four-set score test.
 from_edges=[0]*8
 for u in range(8):
  for v in out[u]:
   W=out[v]&inc[u];k=len(W);from_edges[u]+=k*(k-1)//2;from_edges[v]+=k*(k-1)//2
   for z in W:from_edges[z]+=k-1
 from_sets=[0]*8
 for B in itertools.combinations(range(8),4):
  if sorted(len(out[u]&set(B)) for u in B)==[1,1,2,2]:
   for u in B:from_sets[u]+=1
 need(from_edges==from_sets,'vertex S4 incidence by two independent paths')
 tags=[(len(out[u]),c[u],from_sets[u]) for u in range(8)]
 rho=[(tags[u],tuple(sorted(tags[v] for v in out[u]))) for u in range(8)]
 J4=tuple(sorted(rho));J42=tuple(sorted((rho[u],tuple(sorted(rho[v] for v in out[u]))) for u in range(8)))
 return out,c,e,J4,J42
rows=sorted(csv.DictReader((P/'classes_input_n8.tsv').open(),delimiter='\t'),key=lambda r:int(r['canon']))
D=[(int(r['canon']),data(int(r['canon']))) for r in rows]
def E(dat,p=None,q=None):
 out,c,e,*_=dat
 return tuple(sorted(tuple(sorted(((c[v] if p is None else c[v]%p),(e[u][v] if q is None else e[u][v]%q)) for v in out[u])) for u in range(8)))
def summary(f):
 g=collections.defaultdict(list)
 for w,d in D:g[f(d)].append(w)
 coll=[v for v in g.values() if len(v)>1]
 return {'distinct':len(g),'collision_groups':len(coll),'largest_fiber':max(map(len,g.values())),'first_collision':coll[0] if coll else None}
O={}
O['E_exact']=summary(lambda d:E(d))
O['modular_grid']={f'{p},{q}':summary(lambda d:E(d,p,q)) for p in range(1,14) for q in range(1,8)}
O['neighbor_mod_exact_edge']={str(p):summary(lambda d:E(d,p,None)) for p in range(1,5)}
O['exact_neighbor_edge_mod']={str(q):summary(lambda d:E(d,None,q)) for q in range(1,4)}
O['drop_neighbor']=summary(lambda d:tuple(sorted(tuple(sorted(d[2][u][v] for v in d[0][u])) for u in range(8))))
O['drop_edge']=summary(lambda d:tuple(sorted(tuple(sorted(d[1][v] for v in d[0][u])) for u in range(8))))
O['unpaired']=summary(lambda d:tuple(sorted((tuple(sorted(d[1][v] for v in d[0][u])),tuple(sorted(d[2][u][v] for v in d[0][u]))) for u in range(8))))
O['global_edge_joint']=summary(lambda d:tuple(sorted((d[1][u],d[1][v],d[2][u][v]) for u in range(8) for v in d[0][u])))
O['J4']=summary(lambda d:d[3]);O['J4_second_round']=summary(lambda d:d[4])
need(O['E_exact']['distinct']==6880,'E exact injectivity')
need(O['modular_grid']['4,3']['distinct']==6880,'E43 injectivity')
need(all(O['neighbor_mod_exact_edge'][str(p)]['distinct']<6880 for p in range(1,4)),'p lower bound')
need(all(O['exact_neighbor_edge_mod'][str(q)]['distinct']<6880 for q in range(1,3)),'q lower bound')
need(O['drop_neighbor']['distinct']==6749,'drop neighbor')
need(O['drop_edge']['distinct']==6535,'drop edge')
need(O['unpaired']['distinct']==6871,'unpaired')
need(O['global_edge_joint']['distinct']==6832,'global edge joint')
need(O['J4']['distinct']==6878 and O['J4_second_round']['distinct']==6880,'J4 boundary')
# Structural non-monotonicity warning: numerical modulus increase is not in general refinement.
O['scope']='FINITE-EXACT at n=8. (4,3) is the unique componentwise minimum successful pair among all positive-integer modular pairs, by lower-bound failures with the other coordinate exact; no arbitrary-signature minimality claimed.'
O['positive_integer_modular_family_minimum']=[4,3]
O['explicit_checks']=checks
O['input_sha256']=hashlib.sha256((P/'classes_input_n8.tsv').read_bytes()).hexdigest()
(P/'independent_signature_receipt.json').write_text(json.dumps(O,indent=2,sort_keys=True)+'\n')
print(json.dumps({k:v for k,v in O.items() if k not in ('modular_grid',)},indent=2,sort_keys=True))
