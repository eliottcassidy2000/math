"""Independent dyadic checks from independently reconstructed undirected edges.
Ground vertex0, rather than producer vertex3527; do not import producer functions.
"""
import csv,json,pathlib,hashlib
P=pathlib.Path(__file__).resolve().parent;D=P.parent/'electrical';N=3528
if not D.exists():D=P.parent/'session38-electrical'
checks=0
def need(q,m):
 global checks
 checks+=1
 if not q:raise ValueError(m)
def popbits(z):
 while z:
  b=z&-z;yield b.bit_length()-1;z-=b
def echelon(rows):
 piv={}
 for z in rows:
  while z:
   p=z.bit_length()-1
   if p not in piv:piv[p]=z;break
   z^=piv[p]
 return piv
edges=[];blue_rows=[0]*(N-1)
with (P/'independent_undirected_edges.tsv').open() as f:
 for r in csv.DictReader(f,delimiter='\t'):
  u,v,col,w=map(int,(r['source_node'],r['target_node'],r['blue'],r['literal_count']))
  if u==v:continue
  edges.append((u,v,col,w))
  if col and w%2:
   z=(1<<(u-1) if u else 0)^(1<<(v-1) if v else 0)
   if u:blue_rows[u-1]^=z
   if v:blue_rows[v-1]^=z
piv=echelon(blue_rows);free=[i for i in range(N-1) if i not in piv];basis=[]
for f in free:
 z=1<<f
 for p in sorted(piv):
  if (piv[p]&z).bit_count()%2:z^=1<<p
 need(all((r&z).bit_count()%2==0 for r in blue_rows),'direct blue kernel basis')
 basis.append(z)
K=len(basis);membership=[0]*N
for k,z in enumerate(basis):
 for j in popbits(z):membership[j+1]|=1<<k
q=[0]*K;blue_gram={}
for u,v,col,w in edges:
 d=membership[u]^membership[v]
 if not col:
  need(w%2==0,'literal black twins')
  if (w//2)%2:
   for i in popbits(d):q[i]^=d
 else:
  di=[(i,1 if membership[u]>>i&1 else -1) for i in popbits(d)]
  for i,a in di:
   for j,b in di:blue_gram[i,j]=blue_gram.get((i,j),0)+w*a*b
qb=[0]*K
for (i,j),val in blue_gram.items():
 need(val%2==0,'blue Gram integrality')
 if val//2%2:q[i]^=1<<j;qb[i]^=1<<j
r1=len(echelon(q));need(K==3363 and r1==3358,'independent q1 rank')
O={'status':'PASS','ground_node':0,'blue_reduced_rank':len(piv),'V1_dimension':K,'q1_rank':r1,'V2_dimension':K-r1,'blue_half_form_rank':len(echelon(qb))}
# Audit external explicit lifts by raw integer Laplacian applications.
p=D/'n8_lift_receipt.json';receipt=json.loads(p.read_text());lifts=[];images=[];binary=[]
for record in receipt['witness_lifts']:
 x=[0]*N
 for i in record['binary_one_nodes']:x[i]+=1
 for i in record['correction_one_nodes']:x[i]+=2
 Lx=[0]*N
 for u,v,col,w in edges:
  a=w*(x[u]-x[v]);Lx[u]+=a;Lx[v]-=a
 need(all(v%4==0 for v in Lx),'explicit full mod4 lift')
 need([v//4 for v in Lx]==record['A_lift_divided_by_four'],'full lift-image witness')
 b=sum(((x[i]-x[0])%2)<<(i-1) for i in range(1,N));binary.append(b)
 lifts.append(x);images.append(Lx)
need(len(echelon(binary))==5==K-r1,'five lifts form full V2 basis')
q2=[[sum(x[i]*y[i] for i in range(N))//4%2 for y in images] for x in lifts]
need(q2==receipt['q2_matrix'],'q2 all entries')
r2=len(echelon([sum(v<<j for j,v in enumerate(row)) for row in q2]));need(r2==4,'q2 rank')
O.update({'mod4_lifts_audited':5,'q2_matrix':q2,'q2_rank':r2,'V3_dimension':5-r2,'complete_two_primary':False,'producer_lift_receipt_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'independent_edge_sha256':hashlib.sha256((P/'independent_undirected_edges.tsv').read_bytes()).hexdigest(),'explicit_checks':checks})
(P/'independent_electrical_receipt.json').write_text(json.dumps(O,indent=2)+'\n');print(json.dumps(O,indent=2))
