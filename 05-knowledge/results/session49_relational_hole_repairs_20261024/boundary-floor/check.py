#!/usr/bin/env python3
"""Exact560-subset exclusion and an authenticated native attaining boundary."""
from pathlib import Path
from itertools import combinations
from hashlib import sha256
import argparse,csv,json
P=Path(__file__).resolve().parent
N=0
def need(ok,why):
 global N
 N+=1
 if not ok:raise RuntimeError(why)
def plus(a,b):return 16*((a//16+b//16)%128)+((a^b)&15)
def charge(ss):
 q=0
 for s in ss:q=plus(q,labels[s])
 return q
pairs=list(combinations(range(8),2));chords=[(i,j) for i,j in pairs if j>i+1];full=(1<<21)-1
rev=[chords.index((7-j,7-i)) for i,j in chords]
def matrix(code):
 a=[[0]*8 for _ in range(8)]
 for i,(s,t) in enumerate(pairs):a[s][t]=(code>>(27-i))&1;a[t][s]=1-a[s][t]
 return a
def paths(a):
 def extend(pp,left):
  if not left:yield pp;return
  for v in sorted(left):
   if not pp or a[pp[-1]][v]:yield from extend(pp+(v,),left-{v})
 yield from extend((),set(range(8)))
def word(a,pp):return sum(a[pp[i]][pp[j]]<<k for k,(i,j) in enumerate(chords))
def reflection(w):return sum(((w>>i)&1)<<j for i,j in enumerate(rev))
def orbit(w):return {w,reflection(w),full^w,full^reflection(w)}
parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
for item in json.loads((P/'SOURCE_PINS.json').read_text()):
 raw=(P/item['path']).read_bytes();need(len(raw)==item['bytes'] and sha256(raw).hexdigest()==item['sha256'],'original source bytes')
e=json.loads((P/'inputs/electrical.json').read_text())['order8'];X=e['higher_layers'][-1]['lifts'][0];Y=e['radical_lifts_mod4'];F=[X]+Y;mods=[256]+[4]*5
old=json.loads((P/'inputs/domains.json').read_text())['sources'];records={z['source']:sorted(z['exact_all_shared_clocks_target_domain'],key=lambda r:r['target']) for z in old};pins={z['source']:z['fixed_clocks'] for z in old}
cost={s:sum(a!=b for a,b in zip(*(r['path'] for r in rr))) for s,rr in records.items()}
labels={s:16*(X[s]%128)+sum((Y[j][s]%2)<<(j-1) for j in range(1,5)) for s in records}
A=sorted(s for s in records if cost[s]==2 and s not in {176,1425,2071,1457});triple=[1457,1841,3015]
need(len(A)==16 and all(cost[s]==2 for s in A),'original16cheap base recovered')
need(charge(triple)==1263 and sum(cost[s] for s in triple)==8,'actual triple1263and8clocks')
triples=[{'sources':list(ss),'charge':charge(ss)} for ss in combinations(A,3)]
need(len(triples)==560 and all(r['charge']!=1263 for r in triples),'all560three-base-subsets exclude1263')
# Formal lower-bound arithmetic for all possible left/right intersection sizes.
for common in range(17):
 for left_only in range(17-common):
  for right_only in range(17-common-left_only):
   if left_only-right_only!=3:continue
   h=left_only+right_only
   need(h>=3 and h%2==1,'all existing grade-difference3 disagreements are odd and at least3')
   need(h!=3 or (left_only==3 and right_only==0),'only possible14clock form is a three-subset')
   if h!=3:need(8+2*h>=18 and 3+h>=8,'all other existing crossings pay18clocks8slots')
leftmask,rightmask=11079,2826
left=[s for j,s in enumerate(A) if leftmask>>j&1];right_base=[s for j,s in enumerate(A) if rightmask>>j&1];right=sorted(right_base+triple)
need(len(left)==len(right)==8 and len(right_base)==5,'attaining witness exact V grades')
need(charge(left)==charge(right)==1,'attaining witness exact residual charge1')
changed=set(left)^set(right)
need(len(changed)==8 and sum(cost[s] for s in changed)==18,'attaining actual8slot18clock boundary')
# Independently authenticate every receiver record and guarded menu from canonical
# tournaments; global82-slot capacity and all-dyadic saturation remain inherited.
nodes=list(csv.DictReader((P/'inputs/nodes.tsv').open(),delimiter='\t'));U,V=1965,2234
targets={t:set() for t in (U,V)}
for t in targets:
 for code in {int(nodes[t][k]) for k in ('canon_a','canon_b')}:
  a=matrix(code);targets[t].update(word(a,pp) for pp in paths(a))
need(not targets[U]&targets[V],'actual target classes distinct')
R=A+triple;reserved=set();native={}
for s in R:
 a=matrix(int(nodes[s]['canon_a']));menu=set();kernels=[[],[]]
 for pp in paths(a):
  w=word(a,pp)
  for b,r in enumerate(records[s]):
   if w in orbit(r['mask']):kernels[b].append(list(pp))
  if all(pp[t]==v for t,v in pins[s]):menu.update((pp,w,t) for t in (U,V) if full^w in targets[t])
 need(menu=={(tuple(r['path']),r['mask'],r['target']) for r in records[s]},'exhaustive actual guarded source menu')
 for b,r in enumerate(records[s]):
  pp=r['path'];w=r['mask'];orb=orbit(w);pos={v:i for i,v in enumerate(pp)}
  partner=[[0 if i==j else a[i][j] if abs(pos[i]-pos[j])==1 else a[j][i] for j in range(8)] for i in range(8)]
  need(r['target']==(U,V)[b] and word(a,pp)==w and full^w in targets[r['target']],'actual selected source/target record')
  need(word(partner,pp)==full^w and word(list(zip(*a)),pp[::-1])==reflection(w),'actual inverse and converse paths')
  need(pp in kernels[b] and len(orb)==4 and not reserved&orb,'selected kernel and distinct reserved resource')
  reserved|=orb
 native[str(s)]={'records':records[s],'full_kernels':kernels,'cost':cost[s]}
need(len(reserved)==152,'receiver38reserved inverse-converse units')
# Difference on the receiver leaves all exterior rows unchanged.
def physical(selected):
 d={s:[0]*6 for s in R+[U,V]}
 for s in R:
  t=records[s][int(s in selected)]['target']
  for j,f in enumerate(F):z=2*(f[s]-f[t]);d[s][j]+=z;d[t][j]-=z
 return {s:[z%m for z,m in zip(v,mods)] for s,v in d.items()}
need(physical(left)==physical(right),'all actual original source/target generator rows preserved')
full_witness=[]
for chosen in (set(left),set(right)):
 state={str(s):records[s][int(s in chosen)] for s in sorted(R)}
 need(len({state[str(s)]['target'] for s in triple})==1,'actual triple target equality at both endpoints')
 full_witness.append({'V_sources':sorted(chosen),'records':state})
result={'status':'PASS','checks':N,'base_sources':A,'triple_sources':triple,'triple_charge':1263,'triple_switch_clocks':8,'all_three_base_subsets':triples,'claimed_floor':{'clocks':18,'slots':8,'quantifier':'Every existing opposite-triple-branch receiver crossing at the same grade and charge, with exterior records fixed. No existence claim at other grades.'},'attaining_grade':8,'attaining_charge':1,'left_base_mask':leftmask,'right_base_mask':rightmask,'attaining_full_receiver_records':full_witness,'native_receiver_authentication':native,'inherited_global_authentication':'First native audit manifest c77ee9067796f528f9bb45d3c619d9cdad21f19bcbd70c79a055296ad284dca6; no new full82-slot resource or dyadic saturation theorem.'}
encoded=json.dumps(result,indent=2,sort_keys=True)+'\n'
if args.write:(P/'receipt.json').write_text(encoded)
else:need(json.loads((P/'receipt.json').read_text())==json.loads(encoded),'frozen floor receipt exact')
print(json.dumps({'status':'PASS','checks':N,'three_subsets':560,'floor_clocks':18,'floor_slots':8,'attains_at_grade8_charge1':True,'receipt_sha256':sha256(encoded.encode()).hexdigest()},sort_keys=True))
