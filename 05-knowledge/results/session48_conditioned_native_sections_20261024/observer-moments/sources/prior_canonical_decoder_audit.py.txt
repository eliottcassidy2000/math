#!/usr/bin/env python3
"""Independent literal-port, harmonic-debt, and completion-boundary audit.
No imports from the production checkers; no validation-critical assertions.
"""
import argparse, csv, gzip, json, sys
from array import array
from collections import Counter
from hashlib import sha256
from itertools import combinations, permutations, product
from pathlib import Path
COUNT=0

def require(ok, message):
 global COUNT
 COUNT+=1
 if not ok: raise RuntimeError(message)

def matrix(code,n=8):
 bits=f'{code:0{n*(n-1)//2}b}'; a=[[0]*n for _ in range(n)]
 for (i,j),bit in zip(combinations(range(n),2),bits): a[i][j]=int(bit);a[j][i]=1-int(bit)
 return tuple(map(tuple,a))

def transpose(a):return tuple(zip(*a))
def legal(a,p):return len(p)==len(a) and set(p)==set(range(len(a))) and all(a[u][v] for u,v in zip(p,p[1:]))
def partner(a,p):
 b=[list(r) for r in transpose(a)]
 for u,v in zip(p,p[1:]): b[u][v]=a[u][v];b[v][u]=a[v][u]
 return tuple(map(tuple,b))
def closure(a,p):
 b=partner(a,p); c=transpose(a);d=transpose(b);r=p[::-1]
 require(partner(b,p)==a and partner(c,r)==d,'literal commuting involutions')
 ports=frozenset(((a,p),(b,p),(c,r),(d,r)))
 require(len(ports)==4 and all(legal(x,y) for x,y in ports),'literal four-port free orbit')
 return ports

def all_paths(a):return [p for p in permutations(range(len(a))) if legal(a,p)]
def port_hash(port):return sha256(json.dumps(port,separators=(',',':')).encode()).hexdigest()
def hall(ns):
 return all(len(set().union(*(ns[i] for i in range(len(ns)) if mask>>i&1)))>=mask.bit_count() for mask in range(1,1<<len(ns)))
def assignments(ns):return [a for a in product(*[sorted(s) for s in ns]) if len(a)==len(set(a))]

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--repo',type=Path,required=True);ap.add_argument('--workspace',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
 base=args.workspace
 def workspace_file(logical_path):
  """Resolve original session paths or their stripped publication siblings."""
  logical=Path(logical_path);original=base/logical
  if original.is_file():return original
  if logical.parts and logical.parts[0].startswith('session45-'):
   published=base/Path(logical.parts[0][len('session45-'):],*logical.parts[1:])
   if published.is_file():return published
   raise FileNotFoundError(f'Neither source nor publication layout contains {logical_path}: {original}; {published}')
  raise FileNotFoundError(f'Workspace input does not exist: {original}')
 sourcepins=json.loads(workspace_file('session45-injective-completion/SOURCE_PINS.json').read_text());files={}
 for row in sourcepins['files']+sourcepins['inherited_results']:
  raw=(args.repo/row['path']).read_bytes();require(len(raw)==row['bytes'] and sha256(raw).hexdigest()==row['sha256'],'source-byte pin '+row['path']);files[Path(row['path']).name]=raw
 rows=list(csv.DictReader(files['merged_nodes.tsv'].decode().splitlines(),delimiter='\t'))
 require(len(rows)==3528 and all(int(r['merged_id'])==i for i,r in enumerate(rows)),'literal merged-node indexing')
 merged={int(r[k]):int(r['merged_id']) for r in rows for k in ('class_a','class_b')}
 owner=array('H');owner.frombytes(gzip.decompress(files['owner_map.u16.gz']))
 if sys.byteorder!='little':owner.byteswap()
 require(len(owner)==1<<21,'complete original word owner function')
 chords=[(i,j) for i in range(8) for j in range(i+2,8)];full=(1<<21)-1
 def word(a,p):return sum((1<<b) for b,(i,j) in enumerate(chords) if a[p[i]][p[j]])
 def annotation(a,p):
  require(legal(a,p),'actual Hamiltonian path')
  q=word(a,p);return q,merged[owner[q]],merged[owner[full^q]]
 edges=[tuple(map(int,(r['source'],r['target'],r['blue'],r['black']))) for r in csv.DictReader(gzip.decompress(files['merged_edges.tsv.gz']).decode().splitlines(),delimiter='\t')]
 require(len(edges)==438216,'complete original edge set')
 colors={(u,v):(b,k) for u,v,b,k in edges}
 d=json.loads(files['receipt.json'])['order8'];X=d['higher_layers'][-1]['lifts'][0];Y=d['radical_lifts_mod4'];Z=d['eighth_radical_lifts_mod8'][0]
 def action(x):
  out=[0]*len(rows)
  for u,v,b,k in edges: t=(b+k)*(x[u]-x[v]);out[u]+=t;out[v]-=t
  return out
 lx=action(X);require(all(t%256==0 for t in lx) and any(t%512 for t in lx),'frozen X harmonic modulo256 but not512')
 require(all(t%512==0 for t in action([2*v for v in X])),'scaled 2X genuinely harmonic modulo512')
 for j,y in enumerate(Y):require(all(t%4==0 for t in action(y)),'actual original Y harmonic '+str(j))
 require(all(t%8==0 for t in action(Z)),'actual original Z harmonic')
 # Raw labelled ports and original unit-capacity normalized tilings are different carriers.
 quotient_a=matrix(int(rows[43]['canon_a']));qp=(7,6,5,3,2,1,0,4);qr=(7,6,5,2,1,3,0,4)
 qo1=closure(quotient_a,qp);qo2=closure(quotient_a,qr);qw1={word(A,P) for A,P in qo1};qw2={word(A,P) for A,P in qo2}
 require(annotation(quotient_a,qp)[0]==annotation(quotient_a,qr)[0]==344062,'two hidden owner lifts share one normalized forward word')
 require(qw1==qw2 and len(qw1)==4 and qo1.isdisjoint(qo2),'one normalized black tile resource has two distinct raw labelled K/C lifts')
 quotient_loss={'source':43,'canon':int(rows[43]['canon_a']),'paths':[list(qp),list(qr)],'normalized_word':344062,'normalized_F_R_words':sorted(qw1),'raw_union_ports':8,'normalized_tile_capacity':1,'meaning':'Original metagraph tilings have unit capacity per normalized F/R tile. Raw owner-labelled ports define a different carrier; hidden owner lifts remain sidecars, not additional tiling capacity.'}
 # Entire native hostile fibre reconstructed independently by 8! permutations.
 source=2242;a=matrix(int(rows[source]['canon_a']));paths=all_paths(a);require(len(paths)==407,'source2242 complete path count')
 fibre=[p for p in paths if annotation(a,p)[2]==3442]
 require(fibre==[(1,0,5,3,2,4,6,7),(1,5,4,3,2,0,6,7),(4,1,5,3,2,0,6,7)],'complete hostile target fibre')
 eligible=[p for p in fibre if p[0]==1 and p[3:5]==(3,2) and p[6:8]==(6,7)]
 require(eligible==fibre[:2] and all((X[source]-X[3442])%8==5 for _ in fibre),'exact native pinned neighborhood and phase')
 resources=[closure(a,p) for p in fibre];require(len(set().union(*resources))==12,'hostile three disjoint raw K/C resources')
 qhostile=[frozenset(word(B,P) for B,P in r) for r in resources]
 require(all(len(q)==4 for q in qhostile) and len(set().union(*qhostile))==12,'native Hall hostile also has three distinct pureblack normalized resources')
 require(colors[tuple(sorted((2242,3442)))][0]==0,'native Hall hostile is on an original pureblack edge')
 ns=[{0,1}]*3;require(all(assignments([ns[i],ns[j]]) for i,j in combinations(range(3),2)) and not hall(ns) and not assignments(ns),'all pairs injective but triple impossible')
 # Even a single occupied port excludes its entire closed demand resource.
 one_occupied={next(iter(resources[0]))};available={i for i,r in enumerate(resources) if not r&one_occupied}
 require(available=={1,2} and hall([available]*2) and not hall([available]*3),'nonclosed occupied remainder removes full-orbit capacity')
 normalized=[]
 for rho in (tuple(range(8)),(1,2,3,4,5,6,7,0),(7,6,5,4,3,2,1,0)):
  inv=tuple(rho.index(x) for x in range(8));local=tuple(tuple(a[rho[i]][rho[j]] for j in range(8)) for i in range(8));ports=[]
  for p in eligible:
   q=tuple(inv[v] for v in p);require(legal(local,q),'local owner transport path')
   physical=tuple(tuple(local[inv[i]][inv[j]] for j in range(8)) for i in range(8)),tuple(rho[v] for v in q)
   require(physical==(a,p),'exact physical normalization equality');ports.append(port_hash(physical))
  normalized.append(ports)
 require(normalized[0]==normalized[1]==normalized[2],'local label changes do not create capacity')
 # Independent full Hall-subset criterion vs assignment enumeration, including empty neighborhoods.
 subsets=[{j for j in range(4) if m>>j&1} for m in range(16)];tests=laminartests=0
 for size in range(1,5):
  for ns in product(subsets,repeat=size):
   tests+=1;exist=bool(assignments(ns));require(hall(ns)==exist,'Hall iff finite assignment control')
   if all(not (s&t) or s<=t or t<=s for s,t in combinations(ns,2)):
    laminartests+=1; compressed=all(sum(t<=s for t in ns)<=len(s) for s in set(map(frozenset,ns)))
    require(compressed==exist,'laminar compressed criterion including empty neighborhoods')
 # Boundary hostiles for factors and reservation assumptions.
 matching_global_guard_hostile={'neighbors':[[0],[1]],'edge_lift_values':[[0],[0]],'joint_guard':'lift_0 != lift_1'}
 require(bool(assignments([{0},{1}])) and not any(x!=y for x,y in product([0],[0])),'resource matching cannot certify an independent joint guard')
 require(bool(assignments([{0}])) and not assignments([{0},{0}]),'individually feasible regions can reuse one resource')
 C9=[(u,v) for u,v in product(range(9),repeat=2) if (u-v)%3==0]
 joined=[(a,b,c) for a,b in C9 for bb,c in C9 if b==bb]
 require(len(C9)==27 and len(joined)==81,'all full C9 lifts retained in actual middle-owner join')
 require(not [(a,b,c) for a,b in C9 for bb,c in C9 if b==0 and bb==3 and b==bb],'matching quotient residues cannot glue unequal full owners')
 cycle=[[(0,0),(1,1)],[(0,0),(1,1)],[(0,1),(1,0)]]
 require(not [(a,b,c) for a,b,c in product(range(2),repeat=3) if (a,b) in cycle[0] and (b,c) in cycle[1] and (a,c) in cycle[2]],'cyclic pairwise compatibility has no global row')
 packet_reports=[]
 for kind,file in [('four','session45-multisource-trades/packet.json'),('three','session45-source-circuits/three_source_packet.json')]:
  raw=workspace_file(file).read_bytes();packet=json.loads(raw);states=[[],[]];qresources=[[],[]];maps={};slot_reports=[];physical_sources=[];currents=[Counter(),Counter()]
  records=packet['slots'];target_counter=[Counter(),Counter()]
  for i,slot in enumerate(records):
   s=slot['source'];A=matrix(int(rows[s]['canon_a']));rho=tuple(slot['local_to_physical'] if kind=='four' else slot['local_to_physical_owners']);require(sorted(rho)==list(range(8)),'bijective physical owner map')
   require(s not in maps or maps[s]==rho,'repeated source uses same whole owner map');maps[s]=rho
   inv=tuple(rho.index(j) for j in range(8));physical=tuple(tuple(A[inv[u]][inv[v]] for v in range(8)) for u in range(8));physical_sources.append(physical)
   require(tuple(map(tuple,slot['canonical_matrix'] if kind=='four' else slot['adjacency']))==A,'source canonical matrix byte interpretation')
   if kind=='four':require(tuple(map(tuple,slot['physical_matrix']))==physical,'whole physical source matrix retained')
   pins=slot['all_fixed_physical_clock_owners'] if kind=='four' else slot['preserved_owner_pins'];pair=[];slotphases=[]
   complete=all_paths(A);counts=Counter(annotation(A,p)[2] for p in complete)
   autos=[p for p in permutations(range(8)) if all(A[u][v]==A[p[u]][p[v]] for u,v in combinations(range(8),2))]
   for state,key in enumerate(('before','after') if kind=='four' else ('old','new')):
    record=slot[key];p=tuple(record['path']);pp=tuple(rho[v] for v in p);pair.append(pp)
    if kind=='four':require(pp==tuple(record['physical_path']),'full physical path transport')
    require(all(pp[j]==v for j,v in pins),'all retained indexed owner pins')
    require(all(pp[j]==v for j,v in slot.get('shared_physical_pins',[])),'all explicit cross-source indexed owners')
    q,ss,t=annotation(physical,pp);require((q,ss,t)==(record['mask'],s,record['target']),'actual original source word and target')
    orbit=closure(physical,pp);words={word(B,P) for B,P in orbit};require(len(words)==4,'pureblack size-four normalized F/R orbit')
    require(colors[tuple(sorted((s,t)))][0]==0,'actual original edge has no blue weight')
    # Converse sends a pin at j to exactly 7-j, not back to the same clock.
    require(all(pp[::-1][7-j]==v for j,v in pins),'converse fixed-pin clock transport')
    require({merged[owner[w]] for w in words}=={s,t},'raw quartet has exact original source/target owners')
    states[state].append(orbit);qresources[state].append(min(words));target_counter[state][t]+=1
    phase=(X[s]-X[t])%8;slotphases.append(phase)
    currents[state][s]+=2*(X[s]-X[t]);currents[state][t]-=2*(X[s]-X[t])
   require(slotphases[0]==slotphases[1],'selected X mod8 phase preserved per slot')
   unchanged=[[j,pair[0][j]] for j in range(8) if pair[0][j]==pair[1][j]]
   require(sorted(pins)==unchanged,'every unchanged indexed physical owner is declared')
   slot_reports.append({'source':s,'source_paths':len(complete),'source_automorphisms':len(autos),'source_target_counts':{str(t):counts[t] for t in packet['targets']},'fixed_physical_pins':pins,'phase_mod8':slotphases[0]})
  require(all(len(set().union(*st))==16 for st in states) and len(set().union(*(states[0]+states[1])))==32,'both packet endpoints and their union literally injective')
  require(len(set(qresources[0]+qresources[1]))==8,'eight stronger normalized resource keys distinct')
  require(target_counter[0]==target_counter[1],'target multiplicity preserved')
  def debt(char,indices):
   out=Counter()
   for j in indices:
    slot=records[j];s=slot['source']
    for key,sign in zip(('before','after') if kind=='four' else ('old','new'),(-2,2)):
     t=slot[key]['target'];value=sign*(char[s]-char[t]);out[s]+=value;out[t]-=value
   return {str(k):v for k,v in sorted(out.items()) if v}
  dx=debt(X,range(4));require(dx==({} if kind=='four' else {'2625':256,'2628':-256}),'exact full integral X debt')
  require(all(t%256==0 for t in dx.values()),'frozen X debt vanishes modulo256')
  require(all(t%512==0 for t in debt([2*x for x in X],range(4)).values()),'scaled 2X joint debt vanishes modulo512')
  yreports=[];ordinary_phase_changes=[]
  for j,y in enumerate(Y):
   dy=debt(y,range(4));require(all(v%4==0 for v in dy.values()),'nodewise joint debt for actual Y generator')
   yreports.append(dy)
   for i,slot in enumerate(records):
    before,after=[slot[k]['target'] for k in (('before','after') if kind=='four' else ('old','new'))];change=(y[before]-y[after])%4
    if change:ordinary_phase_changes.append({'slot':i,'character':'Y'+str(j),'modulus':4,'ordinary_phase_change':change})
  sparse=[0]*len(rows);sparse[packet['targets'][0]]=4
  require(all(v%8==0 for v in action(sparse)),'sparse4eU is an actual original mod8 harmonic')
  require(all(v%8==0 for v in debt(sparse,range(4)).values()),'sparse4eU nodewise joint debt vanishes')
  sparse_phase_changes=[(sparse[slot[k0]['target']]-sparse[slot[k1]['target']])%8 for slot in records for k0,k1 in [(('before','after') if kind=='four' else ('old','new'))]]
  require(sparse_phase_changes==[4]*4,'universal ordinary per-port phase preservation is explicitly false')
  rectangles=[{'slots':ids,'integer_X_debt':debt(X,ids)} for ids in ([0,2],[1,3])]
  for r in rectangles:require(any(v%256 for v in r['integer_X_debt'].values()),'constituent rectangle alone fails higher guard')
  allowner_conflicts=[{'owners':[u,v],'source_arc_values':[A[u][v] for A in physical_sources]} for u,v in combinations(range(8),2) if len({A[u][v] for A in physical_sources})>1]
  require(allowner_conflicts,'shared named owners are not a common host amalgam')
  if kind=='four':require(all(A[3][4]==1 for A in physical_sources),'declared cross-source physical arc 3 to4 agrees')
  packet_reports.append({'kind':kind,'input_sha256':sha256(raw).hexdigest(),'slots':slot_reports,'raw_ports_per_endpoint':16,'raw_union_ports':32,'normalized_union_resources':8,'frozen_integer_X_debt':dx,'Y_nodewise_debts':yreports,'ordinary_per_port_phase_changes':ordinary_phase_changes,'sparse_4eU_mod8_per_port_phase_changes':sparse_phase_changes,'constituent_rectangles':rectangles,'shared_owner_arc_conflicts':allowner_conflicts})
 # Full H_256 and H_8 preservation authenticated from image witnesses, not a list of characters.
 wraw=workspace_file('session45-multisource-trades/image_witnesses.json').read_bytes();w=json.loads(wraw);image_reports=[]
 want=[{1425:2,1921:2,1477:-2,2837:-2},{1965:2,2234:-2}]
 for row,charge in zip(w['solutions'],want):
  require(dict(row['charge'])==charge,'image witness exact source/target charge')
  potential=row['potential'];require(len(potential)==3528 and all(type(t)==int for t in potential),'full integer image vector')
  lp=action(potential);require(all((t-charge.get(j,0))%256==0 for j,t in enumerate(lp)),'direct entire Laplacian image equation')
  require(row['integer_correction_nonzero']==[[j,(t-charge.get(j,0))//256] for j,t in enumerate(lp) if t!=charge.get(j,0)],'full integer image correction ledger')
  image_reports.append({'charge':row['charge'],'image_modulus':256,'potential_sha256':sha256(json.dumps(potential,separators=(',',':')).encode()).hexdigest()})
 out={'status':'PASS AFTER EDGEWISE-HYPOTHESIS REPAIR','checks':COUNT,'source_commit':sourcepins['commit'],'source_pins':sourcepins['files'],'resource_carrier_quotient_hostile':quotient_loss,'native_hostile':{'source':2242,'target':3442,'full_source_paths':407,'complete_target_fibre':[list(p) for p in fibre],'eligible':[list(p) for p in eligible],'pairwise_feasible':True,'triple_feasible':False,'normalized_port_hashes':normalized,'nonclosed_remainder_available_resources':sorted(available),'normalized_F_R_resources':[sorted(q) for q in qhostile]},'finite_matching_tests':{'families_including_empty_neighborhoods':tests,'laminar_families':laminartests},'global_guard_hostile':matching_global_guard_hostile,'join_controls':{'full_C9_local_lifts':27,'full_C9_joint_lifts':81,'resource_disjointness_required':True,'running_intersection_required':True},'packets':packet_reports,'four_source_image_certificates':image_reports,'harmonic_boundary':{'raw_X_harmonic_mod256':True,'raw_X_harmonic_mod512':False,'scaled_2X_harmonic_mod512':True,'full_bank_meaning':'nodewise harmonicity, not every ordinary per-port phase','three_source_all_dyadic_scope':'depends on independently audited complete original electrical filtration'}}
 args.out.write_text(json.dumps(out,sort_keys=True,indent=2)+'\n');print(json.dumps({'checks':COUNT,'sha256':sha256(args.out.read_bytes()).hexdigest()}))
if __name__=='__main__':main()
