"""Verify all higher-lift certificates directly against independent integer edges."""
import csv,json,pathlib,hashlib
P=pathlib.Path(__file__).resolve().parent;D=P.parent/'electrical'
if not D.exists():D=P.parent/'session38-electrical'
N=3528
E=[]
for r in csv.DictReader((P/'independent_undirected_edges.tsv').open(),delimiter='\t'):
 u,v,w=map(int,(r['source_node'],r['target_node'],r['literal_count']))
 if u!=v:E.append((u,v,w))
def need(q,m):
 if not q:raise ValueError(m)
def audit(x,image,modulus):
 need(len(x)==len(image)==N,'vector length')
 y=[0]*N
 for u,v,w in E:
  z=w*(x[u]-x[v]);y[u]+=z;y[v]-=z
 need(all(v%modulus==0 for v in y),'integer lift divisibility')
 need([v//modulus for v in y]==image,'entire lift image')
 binary=[(v-x[0])%2 for v in x]
 need(any(binary),'nonzero binary quotient class')
 need(all(0<=v<modulus for v in x),'canonical residue entries')
 q=sum(x[i]*y[i] for i in range(N))//modulus%2
 return q,binary
r3=json.loads((D/'n8_final_radical_receipt.json').read_text())
q,binary=audit(r3['lift_mod8'],r3['A_lift_divided_by_eight'],8)
need(q==r3['q3_matrix'][0][0]==0,'q3 zero')
rr=json.loads((D/'n8_complete_tower_receipt.json').read_text())
print(rr.keys())
layers=[{'k':3,'q':q,'modulus':8}]
need(rr['tail_layers']==[{'k':k,'q':int(k==8)} for k in range(3,9)],'reported tail pattern')
for record in rr['witness_lifts']:
 k=record['k'];qq,bb=audit(record['lift'],record['A_lift_divided_by_modulus'],1<<k)
 need(bb==binary,'same nonzero one-dimensional binary class at each stage')
 need(qq==int(k==8),'direct higher form q value')
 layers.append({'k':k,'q':qq,'modulus':1<<k})
need([v['k'] for v in layers]==list(range(3,9)),'contiguous complete tower')
prior=json.loads((P/'independent_electrical_receipt.json').read_text())
need(prior['V1_dimension']==3363 and prior['V2_dimension']==5 and prior['V3_dimension']==1,'independent initial ranks')
exponents={1:3358,2:4,8:1};valuation=sum(k*v for k,v in exponents.items())
need(rr['two_primary_exponent_counts']=={str(k):v for k,v in exponents.items()},'exact2primary')
need(valuation==3374==rr['spanning_tree_valuation_two'],'exact tree valuation')
O={'status':'PASS','method':'Direct integer Laplacian application to all six frozen lift vectors, independent undirected ledger; no producer elimination imported. Prior independent q1 and q2 rank establish one-dimensional remaining space.','layers':layers,'complete_two_primary':True,'two_primary_exponent_counts':exponents,'spanning_tree_valuation_two':valuation,'witness_coordinates_verified':6*N,'source_sha256':{name:hashlib.sha256((D/name).read_bytes()).hexdigest() for name in ['n8_final_radical_receipt.json','n8_complete_tower_receipt.json']}}
(P/'independent_complete_tower_receipt.json').write_text(json.dumps(O,indent=2)+'\n');print(json.dumps(O,indent=2))
