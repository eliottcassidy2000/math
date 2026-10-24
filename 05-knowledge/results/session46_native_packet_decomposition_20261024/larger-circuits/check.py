"""Literal native validation of large primitive paired-current circuits.
The optional image generator is untrusted by this checker: all Lz equations are multiplied anew.
"""
from pathlib import Path
from itertools import combinations,permutations,product
from array import array
from collections import Counter
import argparse,csv,gzip,json,hashlib,sys
P=argparse.ArgumentParser();P.add_argument('--repo',type=Path,required=True);a=P.parse_args();R=a.repo/'05-knowledge/results';HERE=Path(__file__).resolve().parent
checks=0

def need(v,m):
 global checks
 checks+=1
 if not v:raise RuntimeError(m)

pins=json.loads((HERE/'SOURCE_PINS.json').read_text())
for name,digest in pins['originals'].items():need(hashlib.sha256((a.repo/name).read_bytes()).hexdigest()==digest,'source hash '+name)
need(hashlib.sha256((HERE/'source_inventory.json').read_bytes()).hexdigest()==pins['inherited_inventory_sha256'],'inherited full28 inventory hash')
pool=json.loads((HERE/'pool.json').read_text());
need(pool['options']==json.loads((HERE/'source_inventory.json').read_text())['records'],'pool exactly inherits source inventory')
options={r['source']:r for r in pool['options']};cases=json.loads((HERE/'candidates.json').read_text())['found']
f=R/'metagraph_order8_20261023';nodes=list(csv.DictReader((f/'merged_nodes.tsv').open(),delimiter='\t'))
owner=array('H');owner.frombytes(gzip.decompress((f/'owner_map.u16.gz').read_bytes()))
if sys.byteorder!='little':owner.byteswap()
merged={int(n[k]):int(n['merged_id']) for n in nodes for k in ('class_a','class_b')}
pot=json.loads((R/'metagraph_order8_electrical_20261023/receipt.json').read_text())['order8'];Y=pot['radical_lifts_mod4'];X=pot['higher_layers'][-1]['lifts'][0]
pairs=list(combinations(range(8),2));chords=[(i,j) for i,j in pairs if j>=i+2];idx={e:i for i,e in enumerate(chords)};full=(1<<21)-1

def mirror(w):return sum(((w>>i)&1)<<idx[7-j,7-k] for i,(k,j) in enumerate(chords))
def orbit(w):return {w,full^w,mirror(w),full^mirror(w)}
def mat(code):
 B=[[0]*8 for _ in range(8)]
 for bit,(u,v) in enumerate(pairs):B[u][v]=(code>>(27-bit))&1;B[v][u]=1-B[u][v]
 return B

def partner(B,p):
 E={frozenset(e) for e in zip(p,p[1:])}
 return [[0 if u==v else B[u][v] if frozenset((u,v)) in E else B[v][u] for v in range(8)] for u in range(8)]

def autos(B):
 deg=list(map(sum,B));groups=[[j for j,d in enumerate(deg) if d==deg[i]] for i in range(8)];assignment={};out=[]
 def rec():
  if len(assignment)==8:out.append(tuple(assignment[i] for i in range(8)));return
  i=next(i for i in range(8) if i not in assignment)
  for j in groups[i]:
   if j in assignment.values() or any(B[i][u]!=B[j][v] for u,v in assignment.items()):continue
   assignment[i]=j;rec();del assignment[i]
 rec();return out

U,V=1965,2234
need(X[U]==X[V]==146,'equal literal primitive target values')
need(all((Y[j][U]-Y[j][V])%2==0 for j in range(5)),'paired target parity signatures')
labels={};certs={}
for s,r in options.items():
 B=mat(int(nodes[s]['canon_a']));need(autos(B)==[tuple(range(8))],'literal trivial source automorphism')
 p,q=map(tuple,(r['before']['path'],r['after']['path']));fixed=[j for j in range(8) if p[j]==q[j]];free=[j for j in range(8) if j not in fixed]
 need(len(free)==r['cost'] and len(fixed)>=4,'exact complete fixed-common-clock record')
 completions=[]
 for owners in permutations([p[j] for j in free]):
  pp=list(p)
  for j,v in zip(free,owners):pp[j]=v
  if all(B[x][y] for x,y in zip(pp,pp[1:])):completions.append(tuple(pp))
 completions_with_targets=[]
 for pp in completions:
  ww=sum(B[pp[i]][pp[j]]<<b for b,(i,j) in enumerate(chords));completions_with_targets.append({'path':list(pp),'mask':ww,'target':merged[owner[full^ww]]})
 allowed_completions=[tuple(z['path']) for z in completions_with_targets if z['target'] in [U,V]]
 need(set(allowed_completions)=={p,q},'complete target-margin-compatible native pinned fibre')
 need((s==1923 and len(completions)==3) or (s!=1923 and len(completions)==2),'complete raw native pinned domain')
 start,end,j,c,d=r['pins'];need(1<=j<=5 and len({start,end,c,d})==4,'four distinct shared pin roles')
 need(p[0]==q[0]==start and p[7]==q[7]==end and p[j:j+2]==q[j:j+2]==(c,d),'actual indexed pins')
 shared={start:'shared:first',end:'shared:last',c:'shared:arc0',d:'shared:arc1'}
 owner_map={i:shared.get(i,f'private:{s}:{i}') for i in range(8)}
 need(len(set(owner_map.values()))==8,'injective literal owner naming')
 literal_states=[]
 for rec in (r['before'],r['after']):
  pp=tuple(rec['path']);w=sum(B[pp[i]][pp[j]]<<b for b,(i,j) in enumerate(chords));need(w==rec['mask'],'actual normalized forward mask')
  need(merged[owner[w]]==s and merged[owner[full^w]]==rec['target'],'original source and target owners')
  need(len(orbit(w))==4,'genuine black F/R resource')
  K=partner(B,pp);need(partner(K,pp)==B,'literal partner involution')
  need(all(K[x][y] for x,y in zip(pp,pp[1:])),'retained directed partner path')
  C=list(map(list,zip(*B)));rev=pp[::-1];need(all(C[x][y] for x,y in zip(rev,rev[1:])),'actual converse path')
  need(partner(C,rev)==list(map(list,zip(*K))),'converse partner square')
  cw=sum(C[rev[i]][rev[j]]<<b for b,(i,j) in enumerate(chords));need(cw==mirror(w),'actual reflected converse gauge')
  for clock in fixed:need(owner_map[pp[clock]]==owner_map[p[clock]],'every fixed owner coordinate preserved')
  literal_states.append({'path':list(pp),'clock_by_owner':{owner_map[v]:j for j,v in enumerate(pp)},'source_matrix':B,'partner_matrix':K,'converse_source_matrix':C,'converse_partner_matrix':list(map(list,zip(*K))),'forward_mask':w,'reflected_converse_mask':cw,'normalized_resource':min(orbit(w)),'closed_ports':sorted(orbit(w))})
 labels[s]=tuple(Y[j][s]%2 for j in range(4))+(X[s]%128,)
 stored=pool['labels'][str(s)];need(stored==[sum(labels[s][j]<<j for j in range(4)),labels[s][-1]],'search label exactly matches original native electrical coordinates')
 certs[s]={'source':s,'fixed_clocks':fixed,'fixed_owners':[p[j] for j in fixed],'changed_clocks':free,'owner_map':owner_map,'source_canonical_code':int(nodes[s]['canon_a']),'source_automorphisms':[list(range(8))],'pins':r['pins'],'complete_pinned_paths':completions_with_targets,'literal_kernel_lifts':literal_states,'before':r['before'],'after':r['after'],'resource_ids':[min(orbit(r[k]['mask'])) for k in ('before','after')]}
need(len(set(labels.values()))==len(labels),'distinct source charge labels')
all_native_orbits=[orbit(z['mask']) for r in certs.values() for z in r['complete_pinned_paths']]
need(len(all_native_orbits)==57 and len(set.union(*all_native_orbits))==228,'all57 pinned alternatives have disjoint normalized F/R resources')

raw_edges=(f/'merged_edges.tsv.gz').read_bytes();edges=[tuple(int(row[k]) for k in ('source','target','blue','black')) for row in csv.DictReader(gzip.decompress(raw_edges).decode().splitlines(),delimiter='\t')]
def lap(v):
 out=[0]*3528
 for i,j,b,k in edges:
  delta=(b+k)*(v[i]-v[j]);out[i]+=delta;out[j]-=delta
 return out
chars=[[(64*Y[j][i])%256 for i in range(3528)] for j in range(4)]+[X]
for chi in chars:need(all(z%256==0 for z in lap(chi)),'literal higher harmonic separator')
ordinary=[0]*3528;ordinary[U]=4
need(all(z%8==0 for z in lap(ordinary)),'ordinary4eU is authentic mod8 harmonic')
need((ordinary[U]-ordinary[V])%8==4,'ordinary per-port phase is not preserved')
need(pot['complete_two_primary'] and pot['critical_two_rank']==3363 and pot['exact_valuation_one']==3358 and pot['exact_valuation_two']==4 and pot['exact_valuation_three']==0,'authenticated full two-primary invariant counts')
need([(r['level'],r['dimension'],r['rank'],r['next_dimension']) for r in pot['higher_layers']]==[(4,1,0,1),(5,1,0,1),(6,1,0,1),(7,1,0,1),(8,1,1,0)],'authenticated dyadic filtration exponent256')
images=json.loads((HERE/'image_witnesses.json').read_text());need(images['modulus']==256 and images['input_sha256']==hashlib.sha256(raw_edges).hexdigest(),'image antecedent pin')
for sol in images['solutions']:
 expected=dict(sol['charge']);correction=dict(sol['integer_correction_nonzero']);actual=lap(sol['potential'])
 need(all(z==expected.get(i,0)+256*correction.get(i,0) for i,z in enumerate(actual)),'full integer image witness')
image_charges={tuple(map(tuple,s['charge'])) for s in images['solutions']}
need(tuple(sorted([(U,2),(V,-2)])) in image_charges,'common target charge has explicit image')
reports=[]
for key,c in sorted(cases.items(),key=lambda z:int(z[0])):
 plus=c['plus'];minus=c['minus'];ss=plus+minus;n=len(ss);k=len(plus);sign=[1]*k+[-1]*k
 need(set(plus).isdisjoint(minus) and len(set(ss))==n,'no repeated source')
 charge=tuple(sorted([(s,2*sgn) for s,sgn in zip(ss,sign)]));need(charge in image_charges,'exact circuit charge image supplied')
 old=[options[s]['before' if sg==1 else 'after'] for s,sg in zip(ss,sign)];new=[options[s]['after' if sg==1 else 'before'] for s,sg in zip(ss,sign)]
 old_orbits=[orbit(r['mask']) for r in old];new_orbits=[orbit(r['mask']) for r in new]
 need(len(set.union(*old_orbits))==4*n and len(set.union(*new_orbits))==4*n,'no reused closed resource at either endpoint')
 need(len(set.union(*(old_orbits+new_orbits)))==8*n,'distinct old/new black orbits')
 need(Counter(r['target'] for r in old)==Counter(r['target'] for r in new),'actual target margins')
 need(all((X[s]-X[a['target']])%8==(X[s]-X[b['target']])%8 for s,a,b in zip(ss,old,new)),'original frozen per-slot phase')
 accepted=[];mod8=0;balanced=0;separators=Counter()
 count=[0]*(1<<n);mass=[0]*(1<<n);parity=[0]*(1<<n)
 for bits in range(1,1<<n):
  low=bits&-bits;j=low.bit_length()-1;prev=bits^low;count[bits]=count[prev]+sign[j];mass[bits]=mass[prev]+sign[j]*X[ss[j]];parity[bits]=parity[prev]^sum((Y[t][ss[j]]%2)<<t for t in range(4))
 for bits in range(1<<n):
  if count[bits]:continue
  balanced+=1
  if parity[bits]==0 and mass[bits]%4==0:mod8+=1
  if parity[bits]==0 and mass[bits]%128==0:accepted.append(bits)
  else:
   witness=(parity[bits]&-parity[bits]).bit_length()-1 if parity[bits] else 4
   z=sum(sign[j]*chars[witness][ss[j]] for j in range(n) if bits>>j&1)
   need(2*z%256!=0,'actual nonzero target-owned current excludes proper choice')
   separators[witness]+=1
 need(accepted==[0,(1<<n)-1],'exact two-state high-dyadic guarded menu; primitive atomic support')
 cost=sum(len(certs[s]['changed_clocks']) for s in ss)
 need(cost==c['owner_clock_cost'],'actual endpoint owner-clock distance')
 need(all(len(certs[s]['complete_pinned_paths'])==2 for s in ss),'every selected source is raw-pinned binary without target filtering')
 reports.append({'slots':n,'plus':plus,'minus':minus,'cost':cost,'joint_choices':1<<n,'balanced_choices':balanced,'full_H8_accepted_choices':mod8,'full_H256_accepted_choices':len(accepted),'proper_choice_separators':dict(separators),'frozen_X_debt_at_U':2*sum(sg*X[s] for s,sg in zip(ss,sign)),'closed_ports_per_endpoint':4*n,'closed_ports_union':8*n})
print(json.dumps({'status':'PASS','checks':checks,'literal_source_options':list(certs.values()),'circuits':reports,'scope':'Complete fixed-common-clock binary native menu for selected18; full28 inventory and all57 pinned paths authenticated; full nodewise dyadic bank and chosen per-slot X mod8 preserved, not every ordinary per-slot character or a complete-atlas rewiring.'},indent=2,sort_keys=True))
