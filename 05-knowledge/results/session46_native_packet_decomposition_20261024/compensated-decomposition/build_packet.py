from pathlib import Path
from itertools import combinations
import json,csv,hashlib,argparse
H=Path(__file__).resolve().parent;P=H.parent;ap=argparse.ArgumentParser();ap.add_argument('--repo',type=Path,required=True);args=ap.parse_args();R=args.repo;rows=json.loads((H/'inventory.json').read_text())['records'];nodes=list(csv.DictReader((R/'05-knowledge/results/metagraph_order8_20261023/merged_nodes.tsv').open(),delimiter='\t'))
c=json.loads((H/'circuit16.json').read_text());base=c['plus']+c['minus'];search=json.loads((H/'search-results.json').read_text());wanted=[([3015],0),([1394],0),([1394,2905],0)];routes=[]
for buffers,initial in wanted:
 row=next(r for r in search if r['buffers']==buffers and r['initial']==initial);routes.append(dict(name='buffer_'+'_'.join(map(str,buffers)),sources=base+buffers,**row))
def partner(A,p):
 E={frozenset(e) for e in zip(p,p[1:])}
 return [[0 if i==j else A[i][j] if frozenset((i,j)) in E else A[j][i] for j in range(8)] for i in range(8)]
def transpose(A):return list(map(list,zip(*A)))
slots=[]
for row in rows:
 s=row['source'];code=int(nodes[s]['canon_a']);A=[[0]*8 for _ in range(8)]
 for bit,(i,j) in enumerate(combinations(range(8),2)):A[i][j]=code>>(27-bit)&1;A[j][i]=1-A[i][j]
 first,last,k,a,b=row['pins'];mp={first:0,last:7,a:3,b:4};mp.update(zip(sorted(set(range(8))-mp.keys()),sorted(set(range(8))-set(mp.values()))));f=[mp[i] for i in range(8)];C=[[0]*8 for _ in range(8)]
 for i in range(8):
  for j in range(8):C[f[i]][f[j]]=A[i][j]
 fixed=[[j,u] for j,u in enumerate(row['before']['path']) if row['after']['path'][j]==u];options=[]
 for name in ['before','after']:
  r=row[name];p=[f[i] for i in r['path']];K=partner(C,p);rev=p[::-1]
  options.append(dict(name=name,**r,physical_path=p,literal_closure=[dict(action=action,matrix=B,path=pp) for action,B,pp in [('I',C,p),('F',K,p),('R',transpose(C),rev),('FR',transpose(K),rev)]]))
 slots.append(dict(source=s,physical_owner_symbols=[{0:'shared:first',7:'shared:last',3:'shared:arc0',4:'shared:arc1'}.get(i,f'private:{s}:{i}') for i in range(8)],canonical_matrix=A,physical_matrix=C,local_to_physical=f,shared_pins=row['pins'],shared_physical_clock_owners=[[0,0],[7,7],[k,3],[k+1,4]],fixed_local_clock_owners=fixed,fixed_physical_clock_owners=[[j,f[u]] for j,u in fixed],cost=row['cost'],source_label=row['source_label'],options=options))
packet=dict(targets=[1965,2234],base_source_order=base,source_options=slots,routes=routes,bit_convention='bit 0 selects before at1965; bit1 selects after at2234',kernel_lift_scope='All supplied source automorphism groups are trivial. Literal matrices, physical owner bijections, paths, and all four F/R closure representatives are retained. No hidden lift is counted as capacity.',fixed_remainder='A common unchanged remainder must avoid every normalized resource used at any visited state and must have the required shared-coordinate compatibility. The enlarged initial packet current must be balanced by this same remainder; no global remainder is constructed and the full canonical atlas is not rewired.')
(H/'packet.json').write_text(json.dumps(packet,indent=2,sort_keys=True)+'\n')
charges=[dict(name='common_target',charge={1965:2,2234:-2})]
for route in routes:
 for j,(a,b) in enumerate(zip(route['path'],route['path'][1:])):
  q={s:2*((b>>i&1)-(a>>i&1)) for i,s in enumerate(route['sources']) if ((a^b)>>i)&1};charges.append(dict(name=f"{route['name']}_step_{j+1}",charge=q))
(H/'charges.json').write_text(json.dumps(charges,indent=2,sort_keys=True)+'\n')
files=['05-knowledge/results/metagraph_order8_20261023/'+p for p in ['merged_nodes.tsv','owner_map.u16.gz','merged_edges.tsv.gz']]+['05-knowledge/results/metagraph_order8_electrical_20261023/receipt.json']
pins=dict(originals={p:hashlib.sha256((R/p).read_bytes()).hexdigest() for p in files},inherited_artifacts={str(p.relative_to(P)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [H/'circuit16.json',H/'inventory.json']})
(H/'SOURCE_PINS.json').write_text(json.dumps(pins,indent=2,sort_keys=True)+'\n')
print('packet',len(slots),'routes',len(routes),'charges',len(charges))
