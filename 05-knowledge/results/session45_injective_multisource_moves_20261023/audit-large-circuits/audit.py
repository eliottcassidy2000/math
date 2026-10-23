#!/usr/bin/env python3
"""Independent audit of large native circuits. No producer imports; explicit failures survive -O."""
import argparse,csv,gzip,hashlib,itertools,json,struct
from pathlib import Path
from collections import Counter
A=argparse.ArgumentParser();A.add_argument('--repo',type=Path,required=True);A.add_argument('--packet',type=Path,required=True);A.add_argument('--out',type=Path,required=True);a=A.parse_args()
checks=0

def require(test,message):
 global checks
 checks+=1
 if not test:raise RuntimeError(message)

D=a.repo/'05-knowledge/results';B=D/'metagraph_order8_20261023';P=a.packet
paths={'nodes':B/'merged_nodes.tsv','owners':B/'owner_map.u16.gz','edges':B/'merged_edges.tsv.gz','electrical':D/'metagraph_order8_electrical_20261023/receipt.json','pool':P/'pool.json','candidates':P/'candidates.json','images':P/'image_witnesses.json'}
hashes={k:hashlib.sha256(v.read_bytes()).hexdigest()for k,v in paths.items()}
nodes=[{k:int(v)for k,v in row.items()}for row in csv.DictReader(paths['nodes'].open(),delimiter='\t')]
owner=struct.unpack('<'+'H'*(1<<21),gzip.decompress(paths['owners'].read_bytes()))
merge={row['class_a']:row['merged_id']for row in nodes};merge.update({row['class_b']:row['merged_id']for row in nodes})
require(len(nodes)==3528,'full original vertex count')
edges=[tuple(int(row[k])for k in ('source','target','blue','black'))for row in csv.DictReader(gzip.open(paths['edges'],'rt'),delimiter='\t')]
lookup={(i,j):(b,k)for i,j,b,k in edges}
E=json.loads(paths['electrical'].read_text())['order8'];Y=E['radical_lifts_mod4'];X=E['higher_layers'][-1]['lifts'][0]
records=json.loads(paths['pool'].read_text())['options'];options={r['source']:r for r in records};cases=json.loads(paths['candidates'].read_text())['found']
U,V=1965,2234;full=(1<<21)-1;chords=[(i,j)for i in range(8)for j in range(i+2,8)];chord_index={e:k for k,e in enumerate(chords)}
perms=list(itertools.permutations(range(8)))

def matrix(code):
 upper=[(code>>k)&1 for k in range(27,-1,-1)];out=[[0]*8 for _ in range(8)]
 for (i,j),bit in zip(itertools.combinations(range(8),2),upper):out[i][j]=bit;out[j][i]=1-bit
 return out

def transpose(M):return [list(row)for row in zip(*M)]
def forward_word(M,p):return sum(M[p[i]][p[j]]*(1<<k)for k,(i,j)in enumerate(chords))
def reflected(w):
 bits=[(w>>k)&1 for k in range(21)]
 return sum(bits[chord_index[7-j,7-i]]*(1<<k)for k,(i,j)in enumerate(chords))
def tile(w):return {w,full-w,reflected(w),full-reflected(w)}
def switch(M,p):
 keep={frozenset((p[i],p[i+1]))for i in range(7)}
 return [[0 if i==j else M[i][j]if frozenset((i,j))in keep else 1-M[i][j]for j in range(8)]for i in range(8)]
def hamiltonian(M,p):return len(set(p))==8 and all(M[p[i]][p[i+1]]for i in range(7))
def lap(v):
 out=[0]*3528
 for i,j,b,k in edges:
  flow=(b+k)*(v[i]-v[j]);out[i]+=flow;out[j]-=flow
 return out

require(len(options)==len(records)==20,'twenty distinct actual source choices')
for yy in Y:require(all(v%4==0 for v in lap(yy)),'original order-four harmonic equation')
require(all(v%256==0 for v in lap(X)),'original order-256 harmonic equation')
require(X[U]==X[V]==146,'actual chosen target potential')
require(all((yy[U]-yy[V])%2==0 for yy in Y),'all original order-four paired target phases')
source_report=[];all_tiles=[];saved={}
for source,row in options.items():
 M=matrix(nodes[source]['canon_a']);degree=tuple(map(sum,M));actual=[];automorphisms=[]
 for p in perms:
  if hamiltonian(M,p):actual.append(p)
  if all(degree[i]==degree[p[i]]for i in range(8))and all(M[i][j]==M[p[i]][p[j]]for i in range(8)for j in range(i+1,8)):automorphisms.append(p)
 require(automorphisms==[tuple(range(8))],f'full 40320-permutation automorphism census {source}')
 require(nodes[source]['class_a']!=nodes[source]['class_b'],'nonselfconverse source carrier')
 require(2*len(actual)==nodes[source]['fibre'],'complete source word count matches two ordinary/converse path inventories')
 capacities=Counter(merge[owner[full-forward_word(M,p)]]for p in actual)
 for target in (U,V):require(capacities[target]==row['complete_target_capacities'][str(target)],'independent all-permutation native target capacity')
 p=tuple(row['before']['path']);q=tuple(row['after']['path']);fixed=[i for i in range(8)if p[i]==q[i]]
 require(len(fixed)==6,'six fixed owner clocks')
 complete=[z for z in actual if all(z[i]==p[i]for i in fixed)]
 require(set(complete)=={p,q},'all actual Hamiltonian paths satisfying all six literal clocks are exactly the two endpoints')
 start,end,j,c,d=row['pins'];require(j in range(1,6)and len({start,end,c,d})==4,'valid distinct shared-owner pin roles')
 require(p[0]==q[0]==start and p[-1]==q[-1]==end and p[j:j+2]==q[j:j+2]==(c,d),'shared endpoints and indexed oriented internal arc')
 roles={start:'first',end:'last',c:'arc0',d:'arc1'};owner_names={v:('shared',roles[v])if v in roles else('private',source,v)for v in range(8)}
 require(len(set(owner_names.values()))==8,'injective actual owner occurrence map')
 require(M[c][d]==1,'actual shared arc orientation')
 require(all(owner_names[p[i]]==owner_names[q[i]]for i in fixed),'every fixed clock retained in explicit physical owner map')
 resource=[]
 for record in (row['before'],row['after']):
  z=tuple(record['path']);w=forward_word(M,z);K=switch(M,z);C=transpose(M);KC=transpose(K);rz=z[::-1]
  require(z in actual and w==record['mask'],'actual Hamiltonian path and independent forward word')
  require(merge[owner[w]]==source and merge[owner[full-w]]==record['target'],'literal original atlas source and K-target ownership')
  require(hamiltonian(K,z)and hamiltonian(C,rz)and hamiltonian(KC,rz),'all four K/converse paths valid')
  require(switch(K,z)==M and switch(C,rz)==KC,'literal K involution and converse commuting square')
  words=[forward_word(T,pp)for T,pp in ((M,z),(K,z),(C,rz),(KC,rz))]
  require(words==[w,full-w,reflected(w),full-reflected(w)],'all four actual normalized literal records equal F/R orbit')
  require(len(set(words))==4,'genuine black four-port resource tile')
  colors=lookup[tuple(sorted((source,record['target'])))];require(colors[0]==0 and colors[1]>=2,'original colored edge has real paired black resource')
  resource.append(set(words));all_tiles.append(set(words))
  require((X[source]-X[record['target']])%8==(X[source]-X[U])%8,'original chosen per-slot mod-eight phase')
 source_report.append({'source':source,'all_hamiltonian_paths':len(actual),'automorphisms':len(automorphisms),'target_capacities':{str(t):capacities[t]for t in (U,V)},'fixed_clocks':fixed,'resource_ids':[min(s)for s in resource]})
 saved[source]={'matrix':M,'paths':(p,q),'tiles':resource}
require(sum(map(len,all_tiles))==len(set.union(*all_tiles))==160,'all forty pool resources disjoint, stronger than each circuit separately')

images=json.loads(paths['images'].read_text());require(images['modulus']==256 and images['input_sha256']==hashes['edges'],'image witness input and modulus')
verified_charges={}
for solution in images['solutions']:
 values=solution['potential'];charge={i:v for i,v in solution['charge']};correction={i:v for i,v in solution['integer_correction_nonzero']}
 require(len(values)==3528 and values[-1]==0 and all(type(v)is int and 0<=v<256 for v in values),'full grounded integral image witness schema')
 require(sum(charge.values())==0 and sum(correction.values())==0,'zero-sum image equation charges and correction')
 image=lap(values)
 for i,value in enumerate(image):require(value==charge.get(i,0)+256*correction.get(i,0),'independently multiplied exact integer image equation')
 verified_charges[tuple(sorted(charge.items()))]=values
require(tuple([(U,2),(V,-2)])in verified_charges,'target charge 2v has full image antecedent')

# Read only audited compression evidence: direct checks below use actual harmonic
# vectors; the separate dyadic audit authenticates completeness and exponent 256.
characters256=[[64*y for y in yy]for yy in Y[:4]]+[X]
characters8=[[2*y for y in yy]for yy in Y[:4]]+[X]
for m,bank in ((8,characters8),(256,characters256)):
 for character in bank:require(all(v%m==0 for v in lap(character)),'literal separator harmonic equation')

reports=[]
for key,case in sorted(cases.items(),key=lambda item:int(item[0])):
 positives=case['plus'];negatives=case['minus'];sources=positives+negatives;signs=[1]*len(positives)+[-1]*len(negatives);n=len(sources)
 require(n==case['slots']==2*int(key)and len(set(sources))==n and U not in sources and V not in sources,'distinct source slots disjoint from targets')
 require(sum(signs)==0,'balanced positive and negative partition')
 charge=tuple(sorted((s,2*t)for s,t in zip(sources,signs)));require(charge in verified_charges,'exact combined source charge 2q has verified full image antecedent')
 before=[options[s]['before'if t==1 else'after']for s,t in zip(sources,signs)];after=[options[s]['after'if t==1 else'before']for s,t in zip(sources,signs)]
 require(Counter(z['target']for z in before)==Counter(z['target']for z in after),'literal target margins')
 old_tiles=[tile(z['mask'])for z in before];new_tiles=[tile(z['mask'])for z in after]
 require(len(set.union(*old_tiles))==4*n and len(set.union(*new_tiles))==4*n and len(set.union(*(old_tiles+new_tiles)))==8*n,'actual full K/converse packet resource injection at both endpoints and their union')
 # Build changed conductances from the literal before/after targets. Off-diagonal
 # -Delta w and diagonal row sums recover Delta L without assuming its formula.
 delta=Counter()
 for s,old,new in zip(sources,before,after):delta[tuple(sorted((s,old['target'])))]-=2;delta[tuple(sorted((s,new['target'])))]+=2
 require(len(delta)==2*n and all(abs(v)==2 for v in delta.values()),'actual selected incidence change has exactly two edges per changed source')
 def changed_action(chi,selected=None):
  out=Counter()
  for i,(s,old,new)in enumerate(zip(sources,before,after)):
   if selected is not None and i not in selected:continue
   for target,weight in ((old['target'],-2),(new['target'],2)):
    current=weight*(chi[s]-chi[target]);out[s]+=current;out[target]-=current
  return out
 for chi in characters256:
  computed=changed_action(chi);qdot=sum(t*chi[s]for s,t in zip(sources,signs));vdot=chi[U]-chi[V]
  expected={s:2*t*vdot for s,t in zip(sources,signs)};expected[U]=2*qdot;expected[V]=-2*qdot
  require(all(computed[k]==v for k,v in expected.items()),'actual changed edges agree with rank-two row-owned current formula')
  require(all(v%256==0 for v in computed.values()),'accepted endpoint change preserves each actual higher separator')
 accepted8=[];accepted256=[];balanced=0;separators=Counter();minimal_h8=None
 for bitset in range(1<<n):
  chosen=[j for j in range(n)if bitset&(1<<j)]
  margin=Counter(z['target']for z in before)
  for j in chosen:margin[before[j]['target']]-=1;margin[after[j]['target']]+=1
  if margin!=Counter(z['target']for z in before):continue
  balanced+=1
  ok8=all(all(z%8==0 for z in changed_action(chi,chosen).values())for chi in characters8)
  # Actual harmonic vector, actual row and residue are authenticated for every rejection.
  bad=None
  for k,chi in enumerate(characters256):
   action=changed_action(chi,chosen)
   if any(v%256 for v in action.values()):
    bad=(k,next((i,v%256)for i,v in action.items()if v%256));break
  if ok8:accepted8.append(bitset)
  if bad is None:accepted256.append(bitset)
  else:
   k,(vertex,value)=bad;require(vertex in (U,V)and value!=0,'proper balanced subset has actual nonzero owned-target harmonic current')
   separators[k]+=1
  if ok8 and 0<len(chosen)<n:minimal_h8=min(len(chosen),minimal_h8 or n)
 require(accepted256==[0,(1<<n)-1],'complete full-bank guarded fibre has exactly two endpoints')
 require(len(accepted8)=={6:2,8:2,10:8,12:16,14:68,16:196}[n],'complete H8 state count using independently audited full paired-bank coordinates')
 require(sum(sum(x!=y for x,y in zip(saved[s]['paths'][0],saved[s]['paths'][1]))for s in sources)==2*n,'sharp two-owner-clock minimum per changed permutation')
 # Every H8 vertex was authenticated by literal changed currents above. Determine
 # the exact bottleneck by exhausting all smaller even slot-support thresholds.
 lower_reachable=[];h8_path=None;h8_budget=None
 for budget in range(2,n+1,2):
  parents={0:None};queue=[0]
  for state in queue:
   for nxt in accepted8:
    if nxt not in parents and (state^nxt).bit_count()<=budget:parents[nxt]=state;queue.append(nxt)
  if (1<<n)-1 in parents:
   h8_budget=budget;h8_path=[(1<<n)-1]
   while h8_path[-1]:h8_path.append(parents[h8_path[-1]])
   h8_path.reverse();break
  lower_reachable.append({'slot_budget':budget,'reachable_states_from_before':len(parents)})
 require(h8_budget=={6:6,8:8,10:4,12:4,14:4,16:4}[n],'exact H8-only minimax simultaneous slot budget')
 require(h8_path[0]==0 and h8_path[-1]==(1<<n)-1 and all(x in accepted8 for x in h8_path),'literal authenticated H8 connecting path endpoints and vertices')
 require(max((x^y).bit_count()for x,y in zip(h8_path,h8_path[1:]))==h8_budget,'explicit path attains proven H8 bottleneck')
 reports.append({'slots':n,'all_joint_states':1<<n,'balanced_states':balanced,'H8_accepted':len(accepted8),'H256_accepted':len(accepted256),'H8_minimum_nontrivial_subset_slots':minimal_h8,'H8_minimax_slot_budget':h8_budget,'H8_path_masks':h8_path,'H8_smaller_budget_reachability':lower_reachable,'proper_balanced_rejected':balanced-2,'separator_counts':dict(separators),'ordinary_owner_clock_cost':2*n,'closed_ports_per_endpoint':4*n,'closed_ports_union':8*n,'integer_X_target_debt':2*sum(t*X[s]for s,t in zip(sources,signs))})

report={'status':'PASS','checks':checks,'source_sha256':hashes,'sources':source_report,'circuits':reports,'verified_full_image_equations':len(verified_charges),'all_pool_union_ports':160,'scope':'Twenty sources exhaustively tested on all 40320 owner permutations each. Six complete pinned binary menus exhaustively audited. Full-bank completeness and all-dyadic extension depend on separately independently audited original dyadic filtration; no full-atlas or common-host-amalgam existence asserted.'}
a.out.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');print(json.dumps({'status':'PASS','checks':checks,'output':str(a.out),'circuits':reports},sort_keys=True))
