"""Independent incidence transport, including every fixed-path mask through n7."""
import itertools,json,pathlib
P=pathlib.Path(__file__).resolve().parent
checks=0
def need(q,m):
 global checks
 checks+=1
 if not q:raise ValueError(m)
def counts(out):
 n=len(out);incoming=[set(range(n))-{v}-out[v] for v in range(n)];e=[[0]*n for _ in range(n)]
 for u in range(n):
  for v in out[u]:e[u][v]=e[v][u]=len(out[v]&incoming[u])
 return [sum(e[u][v] for v in out[u]) for u in range(n)],e
def E(out):
 c,e=counts(out);return tuple(sorted(tuple(sorted((c[v],e[u][v]) for v in out[u])) for u in range(len(out))))
rows=[];host=[]
for n in range(3,8):
 free=[(i,j) for i,j in itertools.combinations(range(n),2) if j>i+1];triples=list(itertools.combinations(range(n),3))
 for m in range(1<<len(free)):
  out=[set() for _ in range(n)];target=[set() for _ in range(n)]
  for i in range(n-1):out[i].add(i+1);target[i].add(i+1)
  for k,(i,j) in enumerate(free):
   if m>>k&1:out[i].add(j);target[j].add(i)
   else:out[j].add(i);target[i].add(j)
  c,e=counts(out);actualc,actuale=counts(target);delta_c=[0]*n;delta_e=[[0]*n for _ in range(n)];affected=0
  for tri in triples:
   pairs=list(itertools.combinations(tri,2));path_edges=[(u,v) for u,v in pairs if v==u+1]
   if not path_edges:continue
   affected+=1
   if len(path_edges)==2:
    i,j,k=tri;cyclic=(j in out[i])==(k in out[j])==(i in out[k]);delta=1-2*cyclic
   else:
    u,v=path_edges[0];w=next(z for z in tri if z!=u and z!=v);delta=int(w in out[u])-int(w in out[v])
   for z in tri:delta_c[z]+=delta
   for u,v in pairs:delta_e[u][v]+=delta;delta_e[v][u]+=delta
  need(affected==(n-2)**2,'affected triple count')
  need([x+y for x,y in zip(c,delta_c)]==actualc,'vertex transport')
  need([[e[u][v]+delta_e[u][v] for v in range(n)]for u in range(n)]==actuale,'edge transport')
  if n==4 and m in [0,1]:host.append((E(out),E(target),[sum(1<<v for v in row)for row in out]))
 rows.append({'n':n,'all_fixed_path_masks':1<<len(free)})
need(host[0][0]==host[1][0] and host[0][1]!=host[1][1],'source E insufficient for actual-path target')
O={'status':'PASS','scope':'All fixed-path masks n3..7; independent edge-intersection incidence counts','universes':rows,'total_masks':sum(r['all_fixed_path_masks']for r in rows),'explicit_checks':checks,'path_dependence_hostile_n4_source_rows':[r[2] for r in host]}
(P/'independent_transport_receipt.json').write_text(json.dumps(O,indent=2)+'\n');print(json.dumps(O,indent=2))
