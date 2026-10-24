#!/usr/bin/env python3
"""Independent native authentication. No producer imports; explicit checks survive -O."""
from pathlib import Path
from itertools import permutations,combinations,product
from collections import Counter,defaultdict
from hashlib import sha256
from array import array
import argparse,csv,gzip,json,sys
P=Path(__file__).resolve().parent;I=P/'inputs';N=0

def need(v,m):
 global N
 N+=1
 if not v: raise RuntimeError(m)
def read(n):return json.loads((I/n).read_text())
def digest(p):
 b=p.read_bytes()
 if p.suffix in ('.json','.py','.md','.tsv','.txt'):b=b.replace(b'\r\n',b'\n')
 return sha256(b).hexdigest()
PAIRS=list(combinations(range(8),2));CH=[(i,j) for i in range(8) for j in range(i+2,8)];FULL=(1<<21)-1

def mat(c):
 A=[[0]*8 for _ in range(8)]
 for k,(i,j) in enumerate(PAIRS):A[i][j]=(c>>(27-k))&1;A[j][i]=1-A[i][j]
 return A
def word(A,p):return sum(A[p[i]][p[j]]<<k for k,(i,j) in enumerate(CH))
def reflect(m):return sum(((m>>CH.index((7-j,7-i)))&1)<<k for k,(i,j) in enumerate(CH))
def orbit(m):return {m,reflect(m),FULL^m,FULL^reflect(m)}
def paths(A):
 for p in permutations(range(8)):
  if all(A[a][b] for a,b in zip(p,p[1:])):yield p,word(A,p)
def partner(A,p):
 ix={a:i for i,a in enumerate(p)}
 return [[0 if a==b else A[a][b] if abs(ix[a]-ix[b])==1 else A[b][a] for b in range(8)] for a in range(8)]
def add(a,b):return ((a[0]+b[0])%128,a[1]^b[1])

def main():
 old=read('native_packet_decomposition_20261023/native/domains.json')['sources'];supply=read('compatible_charged_helpers_20261023/supply/selected_supply.json');repair=read('guarded_component_completion_20261024/sections/native_2125_receipt.json');host=read('guarded_component_completion_20261024/native/xor_tree_receipt.json')
 e=read('metagraph_order8_electrical_20261023/receipt.json')['order8'];X=e['higher_layers'][-1]['lifts'][0];Y=e['radical_lifts_mod4'];Z=e['eighth_radical_lifts_mod8'][0]
 for name,h in host['input_hashes'].items():need(digest(I/name)==h,'recursive source hash '+name)
 need(digest(I/'compatible_charged_helpers_20261023/supply/selected_supply.json')==repair['supply_hash'],'native2125 source hash')
 nodes=list(csv.DictReader((I/'metagraph_order8_20261023/merged_nodes.tsv').open(),delimiter='\t'))
 merged={int(r[k]):int(r['merged_id']) for r in nodes for k in ('class_a','class_b')}
 owners=array('H');owners.frombytes(gzip.decompress((I/'metagraph_order8_20261023/owner_map.u16.gz').read_bytes()))
 if sys.byteorder!='little':owners.byteswap()
 labels={s:(X[s]%128,sum((Y[j][s]%2)<<(j-1) for j in range(1,5))) for s in range(len(nodes))}
 # Target classes authenticated from literal canonical tournaments, independently of owner map.
 target_words={t:set() for t in (1965,2234)}
 for t in target_words:
  for c in set(int(nodes[t][k]) for k in ('canon_a','canon_b')):
   target_words[t].update(m for p,m in paths(mat(c)))
 need(not target_words[1965]&target_words[2234],'two actual target classes distinct')
 data={r['source']:{'records':sorted(r['exact_all_shared_clocks_target_domain'],key=lambda x:x['target']),'fixed':r['fixed_clocks'],'blocks':[]} for r in old}
 cores=sorted(data)
 for r in supply['helpers']:
  s=r['selected'];data[r['source']]={'records':[s['before'],s['after']],'fixed':s['common_clocks'],'blocks':[]}
 data[2125]={'records':repair['records'],'fixed':[],'blocks':repair['retained_blocks']}
 frames={368:list(range(8)),776:[8,9,4,10,11,12,13,3],3139:[14,15,16,17,18,3,4,19],2125:list(range(20,28))};next_owner=28
 for s in sorted(set(data)-set(frames)):frames[s]=list(range(next_owner,next_owner+8));next_owner+=8
 need(next_owner==652 and len(data)==82,'actual82 in652-owner source frame')
 need(frames=={int(k):v for k,v in host['source_frames'].items()},'host frames rebuilt without adopting certificate')
 for a,b in combinations(data,2):need(len(set(frames[a])&set(frames[b]))==(2 if {a,b}<={368,776,3139} else 0),'private core and exact tree intersections')
 mats={s:mat(int(nodes[s]['canon_a'])) for s in data};occupied=set();cert={};costs={};source_arcs={};edge_units={}
 for s,d in sorted(data.items()):
  A=mats[s];em=frames[s]
  for i,j in PAIRS:
   a,b=em[i],em[j];edge=tuple(sorted((a,b)));v=A[i][j] if a<b else A[j][i]
   need(edge not in source_arcs or source_arcs[edge]==v,'coherent source arc');source_arcs[edge]=v
  allpaths=list(paths(A)); menu=[]; kernels=[[],[]];orbits=[orbit(r['mask']) for r in d['records']]
  for p,m in allpaths:
   need(merged[owners[m]]==s,'source map agrees with canonical relabeling')
   for bit,t in enumerate((1965,2234)):
    need((FULL^m in target_words[t])==(merged[owners[FULL^m]]==t),'target classification independently canonical')
    if m in orbits[bit]:kernels[bit].append(list(p))
   guard=all(p[t]==v for t,v in d['fixed']) and all(any(list(p[i:i+len(b)])==b for i in range(9-len(b))) for b in d['blocks'])
   if guard:
    for t in (1965,2234):
     if FULL^m in target_words[t]:menu.append((p,m,t))
  need(set(menu)=={(tuple(r['path']),r['mask'],r['target']) for r in d['records']},'full guarded menu exactly declared two records')
  for bit,r in enumerate(d['records']):
   p=r['path'];m=r['mask'];t=(1965,2234)[bit]
   need(r['target']==t and word(A,p)==m and all(A[a][b] for a,b in zip(p,p[1:])),'literal chosen full path')
   need(word(partner(A,p),p)==FULL^m,'actual inverse partner')
   op=list(map(list,zip(*A)));need(word(op,p[::-1])==reflect(m),'actual converse word')
   need(len(orbits[bit])==4 and not occupied&orbits[bit],'global reserved resource unit capacity');occupied.update(orbits[bit])
   need((X[s]-X[t])%8==r['phase'],'chosen scalar phase')
   need(list(p) in kernels[bit],'selected lift retained, never substituted')
  for t in (1965,2234):edge_units[tuple(sorted((s,t)))]=len({min(orbit(m)) for p,m in allpaths if FULL^m in target_words[t]})
  costs[s]=sum(a!=b for a,b in zip(*(r['path'] for r in d['records'])))
  cert[str(s)]={'hamiltonian_paths':len(allpaths),'guarded_target_records':len(menu),'full_resource_kernel_lifts':kernels,'cost':costs[s],'normalized_resources':[min(z) for z in orbits]}
 need(len(occupied)==656,'164 reserved units,82 occupied at every choice')
 for s in (368,776,3139,2125):need(cert[str(s)]['full_resource_kernel_lifts']==host['full_four_cell_kernels'][str(s)]['full_canonical_resource_lifts'],'independent full kernel agrees with source certificate')
 truth=[]
 for bits in product((0,1),repeat=3):
  owner_values=[];partner_arcs={}
  for s,bit,clock in zip((368,776,3139),bits,(1,4,7)):
   p=data[s]['records'][bit]['path'];gp=[frames[s][v] for v in p];owner_values.append(gp[clock]);B=partner(mats[s],p)
   for i,j in PAIRS:
    a,b=frames[s][i],frames[s][j];edge=tuple(sorted((a,b)));v=B[i][j] if a<b else B[j][i]
    need(edge not in partner_arcs or partner_arcs[edge]==v,'all tree partner states glue');partner_arcs[edge]=v
  allowed=len(set(owner_values))==1;need(allowed==(bits in ((0,1,1),(1,0,0))),'actual owner relation, not an abstract bit substitute')
  truth.append({'bits':bits,'owner_values':owner_values,'allowed':allowed})
 need([r['bits'] for r in truth if r['allowed'] and r['owner_values'][0]==4]==[(0,1,1)],'external owner4 excludes owner3 complete phase')
 free=sorted(set(data)-set(cores)-{368,776,3139});pair_charges=Counter(add(labels[a],labels[b]) for a,b in combinations(free,2))
 need(len(free)==51 and sum(pair_charges.values())==1275 and len(pair_charges)==972,'actual owner4 helper pair holes')
 need(sum(sorted((costs[s] for s in free),reverse=True)[:2])==15,'helper edge cost15')
 need(sum(sorted((costs[s] for s in free),reverse=True)[:4])==29,'direct helper endpoint cost29')
 # Literal original potential equations, retaining black multiplicity two.
 F=[X]+Y;mods=[256]+[4]*len(Y);L=[[0]*len(nodes) for _ in F];edges=0
 with gzip.open(I/'metagraph_order8_20261023/merged_edges.tsv.gz','rt') as f:
  for row in csv.DictReader(f,delimiter='\t'):
   a,b=int(row['source']),int(row['target']);w=int(row['blue'])+int(row['black']);edges+=1
   if (a,b) in edge_units:
    need(int(row['blue'])==0 and int(row['black'])==2*edge_units.pop((a,b)),'literal black weight equals twice independently enumerated full resource units')
   for j,v in enumerate(F):delta=w*(v[a]-v[b]);L[j][a]+=delta;L[j][b]-=delta
 need(not edge_units,'all164 native source/target capacity rows located')
 for j,mod in enumerate(mods):need(all(v%mod==0 for v in L[j]),'every physical row harmonic '+str(j))
 need(all((X[s]%2)==((Y[0][s]+Y[4][s])%2) for s in range(len(nodes))),'Y0 redundancy needed by paired Q adapter')
 need(all(2*(v[1965]-v[2234])%m==0 for v,m in zip(F,mods)),'all source rows unchanged by U/V at paired weight2')
 result={'status':'PASS independent native full records and original harmonic rows; saturation theorem inherited explicitly','source_commit':'c09b642a629b7d12b529b88c2697abaee78e73e7','slots':82,'host_owners':652,'reserved_units':164,'reserved_closed_words':656,'arrival_units':82,'arrival_closed_words':328,'cores':cores,'free_helpers':free,'costs':costs,'labels':{s:16*labels[s][0]+labels[s][1] for s in data},'full_native_records':cert,'tree_truth':truth,'helper_pair_charge_count':972,'helper_edge_clock_bound':15,'helper_direct_clock_bound':29,'physical_edges':edges,'checks':N}
 return json.loads(json.dumps(result))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');a=ap.parse_args();result=main();out=P/'native_certificate.json'
 if a.write:out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
 else:need(result==json.loads(out.read_text()),'frozen native certificate exact')
 print(json.dumps({'status':'PASS','checks':N,'sources':82,'all_full_menus':True,'all_full_kernel_lifts':True,'owner4_helper_image':972},sort_keys=True))
