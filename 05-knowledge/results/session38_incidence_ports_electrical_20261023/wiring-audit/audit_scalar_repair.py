"""Check the four exact J4-hostile classes' associated overlap scalar."""
import itertools,json,pathlib
P=pathlib.Path(__file__).resolve().parent
def need(q,m):
 if not q:raise ValueError(m)
def get(w):
 out=[set() for _ in range(8)]
 for k,(u,v) in enumerate(itertools.combinations(range(8),2)):
  if w>>(27-k)&1:out[u].add(v)
  else:out[v].add(u)
 incoming=[set(range(8))-{u}-out[u] for u in range(8)];e=[[0]*8 for _ in range(8)]
 for u in range(8):
  for v in out[u]:e[u][v]=e[v][u]=len(out[v]&incoming[u])
 c=[sum(e[u][v] for v in out[u])for u in range(8)];v4=[0]*8
 for B in itertools.combinations(range(8),4):
  if sorted(len(out[u]&set(B))for u in B)==[1,1,2,2]:
   for u in B:v4[u]+=1
 tags=[(len(out[u]),c[u],v4[u])for u in range(8)]
 J4=sorted((tags[u],sorted(tags[v] for v in out[u]))for u in range(8))
 omega=[sum(e[u][v]*(e[u][v]-1)//2 for v in out[u])for u in range(8)]
 Q=sum(omega[u]*omega[v]**2 for u in range(8) for v in out[u])
 return J4,sorted(e[u][v]for u,v in itertools.combinations(range(8),2)),Q
pairs=[]
for (a,b),values in [((6427778,6427781),(2881,2887)),((27467784,27467793),(3008,3000))]:
 x,y=get(a),get(b);need(x[:2]==y[:2],'same J4 and entire edge histogram');need((x[2],y[2])==values,'Q values')
 pairs.append({'canons':[a,b],'Q':values,'J4_and_edge_histogram_equal':True})
r=json.loads((P/'independent_signature_receipt.json').read_text());need(r['J4']['collision_groups']==2 and r['J4']['largest_fiber']==2,'only two pairs in complete census')
O={'status':'PASS','pairs':pairs,'conclusion':'Together with the independent complete J4 census, the scalar Q separates every remaining order-eight J4 collision.'}
(P/'independent_scalar_receipt.json').write_text(json.dumps(O,indent=2)+'\n');print(json.dumps(O,indent=2))
