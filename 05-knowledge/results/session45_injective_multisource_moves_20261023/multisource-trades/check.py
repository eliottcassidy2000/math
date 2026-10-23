#!/usr/bin/env python3
"""Independent exact checker: actual paths, physical pins, F/R resources, and full image witnesses."""
import argparse,csv,gzip,json,sys
from array import array
from collections import Counter,defaultdict
from hashlib import sha256
from itertools import combinations,permutations,product
from pathlib import Path
HERE=Path(__file__).resolve().parent;CHECKS=0

def need(test,detail):
 global CHECKS
 CHECKS+=1
 if not test:raise RuntimeError(detail)
def transpose(A):return [list(r) for r in zip(*A)]
def key(A,p):return tuple(tuple(r) for r in A),tuple(p)
def complement(A,p):
 B=[r[:] for r in A];adj={frozenset((u,v)) for u,v in zip(p,p[1:])}
 for u,v in combinations(range(8),2):
  if frozenset((u,v)) not in adj:B[u][v],B[v][u]=B[v][u],B[u][v]
 return B

def main():
 p=argparse.ArgumentParser();p.add_argument('--repo',type=Path);p.add_argument('--out',type=Path);a=p.parse_args()
 repo=a.repo or next((q for q in HERE.parents if (q/'05-knowledge/results/metagraph_order8_20261023/merged_nodes.tsv').is_file()),None)
 if repo is None:raise RuntimeError('Supply --repo containing original native atlas.')
 pins=json.loads((HERE/'SOURCE_PINS.json').read_text());blobs={}
 for f in pins['files']:
  b=(repo/f['path']).read_bytes();need(sha256(b).hexdigest()==f['sha256'],'exact original source hash '+f['path']);blobs[Path(f['path']).name]=b
 nodes=list(csv.DictReader(blobs['merged_nodes.tsv'].decode().splitlines(),delimiter='\t'));owner=array('H');owner.frombytes(gzip.decompress(blobs['owner_map.u16.gz']))
 if sys.byteorder!='little':owner.byteswap()
 merged={int(r[k]):int(r['merged_id']) for r in nodes for k in ('class_a','class_b')};N=len(nodes);need(N==3528,'native vertex universe')
 edges=[tuple(int(r[k]) for k in ('source','target','blue','black')) for r in csv.DictReader(gzip.decompress(blobs['merged_edges.tsv.gz']).decode().splitlines(),delimiter='\t')];colored={(i,j):(b,k) for i,j,b,k in edges}
 D=json.loads(blobs['receipt.json'])['order8'];X=D['higher_layers'][-1]['lifts'][0];Y=D['radical_lifts_mod4'];Z=D['eighth_radical_lifts_mod8'][0]
 def apply(x):
  y=[0]*N
  for i,j,b,k in edges:d=(b+k)*(x[i]-x[j]);y[i]+=d;y[j]-=d
  return y
 need(all(z%256==0 for z in apply(X)),'primitive original potential harmonicmod256')
 sparse=[0]*N;sparse[1965]=4;need(all(z%8==0 for z in apply(sparse)),'actual sparse4e1965 harmonic character')
 need((sparse[1965]-sparse[2234])%8==4,'universal individual port phases are NOT retained')
 separator=[2*z%8 for z in Y[2]];need(all(z%8==0 for z in apply(separator)),'actual original separating character2Y2')
 W=json.loads((HERE/'image_witnesses.json').read_text());need(W['modulus']==256 and W['input_sha256']==sha256(blobs['merged_edges.tsv.gz']).hexdigest(),'image witness exact matrix and modulus')
 need(len(W['solutions'])==2,'exactly two image witnesses')
 wanted=[{1425:2,1921:2,1477:-2,2837:-2},{1965:2,2234:-2}];image_reports=[]
 for row,want in zip(W['solutions'],wanted):
  need(dict(row['charge'])==want,'exact signed charge');z=row['potential'];need(len(z)==N and z[-1]==0 and all(type(t) is int and 0<=t<256 for t in z),'full grounded integral image witness');y=apply(z)
  need(all((v-want.get(i,0))%256==0 for i,v in enumerate(y)),'entire original Laplacian image equation')
  need(row['integer_correction_nonzero']==[[i,(v-want.get(i,0))//256] for i,v in enumerate(y) if v!=want.get(i,0)],'exact integer correction ledger')
  image_reports.append(dict(charge=row['charge'],modulus=256,vector_sha256=sha256(json.dumps(z,separators=(',',':')).encode()).hexdigest()))
 packet=json.loads((HERE/'packet.json').read_text());need(packet['targets']==[1965,2234] and len(packet['slots'])==4,'exact declared packet')
 chords=[(i,j) for i in range(8) for j in range(i+2,8)];pairs=list(combinations(range(8),2));full=(1<<21)-1;indices={e:i for i,e in enumerate(chords)}
 def word(A,p):return sum(A[p[i]][p[j]]<<bit for bit,(i,j) in enumerate(chords))
 def reflect(q):return sum(((q>>indices[(7-j,7-i)])&1)<<bit for bit,(i,j) in enumerate(chords))
 raw_ports={};tile_resources={};slot_reports=[];fullpath_counts={};states={b:[] for b in ('before','after')}
 for slot,r in enumerate(packet['slots']):
  s=r['source'];canon=int(nodes[s]['canon_a']);A=[[0]*8 for _ in range(8)]
  for bit,(i,j) in enumerate(pairs):A[i][j]=canon>>(27-bit)&1;A[j][i]=1-A[i][j]
  need(A==r['canonical_matrix'] and canon==r['canon'],'canonical source matrix retained')
  f=r['local_to_physical'];need(sorted(f)==list(range(8)),'explicit bijective common physical-owner map');C=[[0]*8 for _ in range(8)]
  for i in range(8):
   for j in range(8):C[f[i]][f[j]]=A[i][j]
  need(C==r['physical_matrix'],'actual relabelled physical source matrix retained')
  paths=[]
  def dfs(p,left):
   if not left:paths.append(p);return
   for v in sorted(left):
    if not p or A[p[-1]][v]:dfs(p+(v,),left-{v})
  dfs((),set(range(8)));brute=[p for p in permutations(range(8)) if all(A[u][v] for u,v in zip(p,p[1:]))];need(paths==brute,'independent DFS and all40320 permutations agree')
  autos=[g for g in permutations(range(8)) if all(A[i][j]==A[g[i]][g[j]] for i,j in pairs)];need(autos==[tuple(range(8))],'actual source automorphism group trivial')
  fullpath_counts[s]=len(paths);need(2*len(paths)==int(nodes[s]['fibre']),'native nonselfconverse paired path inventory');need(nodes[s]['class_a']!=nodes[s]['class_b'],'ordinary and converse source owners distinct')
  fp=r['all_fixed_local_clock_owners'];need(len(fp)==6,'all six unchanged literal owner clocks pinned')
  eligible=[list(p) for p in paths if all(p[k]==v for k,v in fp)];need(sorted(eligible)==sorted([r['before']['path'],r['after']['path']]),'complete six-pin family exactly two actual paths')
  targetcounts=Counter()
  for pp in paths:
   qq=word(A,pp);need(owner[qq]==int(nodes[s]['class_a']),'every complete native path retains the actual ordinary source owner');targetcounts[merged[owner[full^qq]]]+=1
  pins_report=[];rows=[]
  for state in ('before','after'):
   row=r[state];p=row['path'];physical=row['physical_path'];need(p in eligible and physical==[f[v] for v in p],'actual complete eligible local path and mapped owner path')
   need(all(C[u][v] for u,v in zip(physical,physical[1:])),'all seven actual physical arcs')
   need(all(physical[k]==v for k,v in r['shared_physical_pins']),'every declared shared indexed owner retained')
   need(all(physical[k]==v for k,v in r['all_fixed_physical_clock_owners']),'every remaining actual fixed clock retained')
   need(C[3][4]==1,'common actual shared arc3to4 retained')
   q=word(C,physical);t=row['target'];need(q==row['mask'] and owner[q]==int(nodes[s]['class_a']) and merged[owner[full^q]]==t,'original forward gauge and target owner authenticated')
   color=colored[tuple(sorted((s,t)))];need(color==(0,2*targetcounts[t]),'pure black factor2 equals actual ordinary path inventory')
   K=complement(C,physical);R=transpose(C);KR=transpose(K);rev=physical[::-1]
   need(complement(K,physical)==C and complement(R,rev)==KR,'literal K involution and actual converse square')
   quartet=[(C,physical),(K,physical),(R,rev),(KR,rev)];words=[q,full^q,reflect(q),full^reflect(q)];need(len(set(words))==4,'native pure black size-four F/R resource tile')
   resource=min(words);need(resource not in tile_resources,'no selected before/after resource tile reused across sources');tile_resources[resource]=(slot,state)
   closure=[]
   for (B,w),mask in zip(quartet,words):
    need(all(B[u][v] for u,v in zip(w,w[1:])),'every actual K/converse closure path Hamiltonian');need(word(B,w)==mask,'literal closure equals F/R word action')
    literal=key(B,w);need(literal not in raw_ports,'no literal matrix-and-owner-path raw port reused');raw_ports[literal]=(slot,state)
    closure.append(dict(matrix=B,path=w,forward_word=mask,ordinary_source=int(owner[mask]),merged_source=merged[owner[mask]],merged_target=merged[owner[full^mask]]))
   states[state].append((s,t,resource));rows.append(dict(state=state,target=t,phase_mod8=(X[s]-X[t])%8,resource_tile=sorted(words),resource_id=resource,original_black_weight=color[1],complete_ordinary_target_capacity=targetcounts[t],closure=closure))
  need((X[r['before']['target']]-X[r['after']['target']])%8==0,'individual original phase retained')
  cost=sum(x!=y for x,y in zip(r['before']['path'],r['after']['path']));need(cost==2,'actual costtwo primitive reorder')
  slot_reports.append(dict(source=s,ordinary_source=int(nodes[s]['class_a']),all_hamiltonian_paths=len(paths),exact_six_pin_family=eligible,cost=cost,states=rows))
 need(len(raw_ports)==32 and len(tile_resources)==8,'two endpoints use32 distinct raw closure ports and8 distinct black tiles')
 need(Counter(t for s,t,q in states['before'])==Counter(t for s,t,q in states['after'])=={1965:2,2234:2},'actual target margins retained')
 signed_edges=defaultdict(int)
 for state,sgn in [('before',-2),('after',2)]:
  for s,t,q in states[state]:signed_edges[tuple(sorted((s,t)))]+=sgn
 need(len(signed_edges)==8 and all(abs(v)==2 for v in signed_edges.values()),'nonzero incidence change on8 distinct original edges')
 for (s,t),c in signed_edges.items():need(colored[s,t][1]+c>=0,'actual edge nonnegativity; not used as distinct-resource proof')
 def delta(x,changed):
  y=[0]*N
  for i in changed:
   r=packet['slots'][i];s=r['source']
   for state,c in [('before',-2),('after',2)]:
    t=r[state]['target'];z=c*(x[s]-x[t]);y[s]+=z;y[t]-=z
  return y
 need(not any(delta(X,range(4))),'chosen frozen modulus256 lift has EXACT integer debtzero')
 q=[0]*N;v=[0]*N
 for r in packet['slots']:q[r['source']]=r['sign']
 v[1965]=1;v[2234]=-1
 for i in range(N):
  if not q[i] and not v[i]:continue
  for j in range(N):
   if not q[j] and not v[j]:continue
   expected=2*(q[i]*v[j]+v[i]*q[j]);actual=-signed_edges.get(tuple(sorted((i,j))),0) if i!=j else sum(c for (s,t),c in signed_edges.items() if i in (s,t));need(actual==expected,'exact signed-matrix identity DeltaL=2(qvT+vqT)')
 accepted=[];partial=[]
 for bits in product((0,1),repeat=4):
  chosen=[i for i,b in enumerate(bits) if b];targets=[r['after' if bits[i] else 'before']['target'] for i,r in enumerate(packet['slots'])];balanced=Counter(targets)=={1965:2,2234:2};debt=delta(separator,chosen);passes=balanced and all(z%8==0 for z in debt)
  if passes:accepted.append(list(bits))
  partial.append(dict(changed=list(bits),balanced=balanced,separator_debt_nonzero=[[i,z%8] for i,z in enumerate(debt) if z%8]))
 need(accepted==[[0,0,0,0],[1,1,1,1]],'complete16-state pinned packet has exactly two universallymod8-safe balanced endpoints')
 rectangles=[]
 for ii in ([0,2],[1,3]):
  dx=delta(X,ii);dy=delta(separator,ii);need(any(z%8 for z in dy) and any(z%256 for z in dx),'each rectangle separately fails universalmod8 and chosenmod256')
  rectangles.append(dict(slots=ii,frozen_integer_debt=[[i,z] for i,z in enumerate(dx) if z],harmonic_separator_mod8_debt=[[i,z%8] for i,z in enumerate(dy) if z%8]))
 out=dict(status='FINITE-EXACT INJECTIVE FOUR-SOURCE INCIDENCE-CHANGING CIRCUIT',checks=CHECKS,original_edges=len(edges),source_hashes=pins['files'],slot_reports=slot_reports,distinct_resource_tiles_both_endpoints=len(tile_resources),distinct_raw_closure_ports_both_endpoints=len(raw_ports),per_endpoint_ordinary_slots=4,per_endpoint_converse_paired_observers=8,per_endpoint_full_K_converse_raw_ports=16,total_actual_changed_owner_clocks=8,incidence_delta_edges=[[i,j,c] for (i,j),c in sorted(signed_edges.items())],image_certificates=image_reports,rectangles=rectangles,primitive_16_state_test=partial,accepted_states=accepted,frozen_values={str(i):X[i] for i in [1425,1921,1477,2837,1965,2234]},frozen_integer_debt=[],per_slot_universal_guard_hostile=dict(character='4e1965 mod8',target_phase_difference=4,conclusion='Original per-slot X phase is retained; universal per-slot phase is not. Full harmonic-bank claims refer to aggregate node currents.'),scope='Selected actual finite packets and their literal K/converse closures. Full original H256 and H8 are preserved by exact image witnesses; old/new kernels need not be equal. An unchanged global remainder must avoid all8 before/after resource tiles; the full canonical all-path atlas is not rewired.')
 text=json.dumps(out,indent=2,sort_keys=True)+'\n'
 if a.out:a.out.write_text(text)
 else:print(text,end='')
if __name__=='__main__':main()
