#!/usr/bin/env python3
"""Independent Session47 finite helper resilience audit; no assertion statements."""
import argparse,csv,gzip,hashlib,itertools,json
from pathlib import Path
from collections import Counter
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
parser=argparse.ArgumentParser();parser.add_argument('--repo-root',type=Path,default=ROOT/'session35-publication/clean-repo');parser.add_argument('--write',action='store_true');args=parser.parse_args()
RESULTS=args.repo_root/'05-knowledge/results'
SUPPLY=RESULTS/'compatible_charged_helpers_20261023/supply/selected_supply.json'
P53=HERE/'quad_packings_initial9.json'
checks=0
def require(x,msg):
 global checks
 checks+=1
 if not x:raise RuntimeError(msg)
def dump(p,x):p.write_text(json.dumps(x,indent=2,sort_keys=True)+'\n')
def sha(p):return hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
def add(a,b):return (((a//16+b//16)%128)*16)|((a^b)&15)
def qsum(vals):
 ans=0
 for val in vals:ans=add(ans,val)
 return ans
C=json.loads(SUPPLY.read_text());helpers={h['source']:h for h in C['helpers']};ids=sorted(helpers)
E=json.loads((RESULTS/'metagraph_order8_electrical_20261023/receipt.json').read_text())['order8']
X=E['higher_layers'][-1]['lifts'][0];Y=E['radical_lifts_mod4']
labels={s:16*(X[s]%128)+sum((Y[j][s]%2)<<(j-1) for j in range(1,5)) for s in ids}
require(len(ids)==53,'53 distinct actual owners')
for s,h in helpers.items():
 require(h['label']==[X[s]%128]+[Y[j][s]%2 for j in range(1,5)],'original source potentials differ')
 require(h['selected']['cost']==sum(a!=b for a,b in zip(h['selected']['before']['path'],h['selected']['after']['path'])),'actual owner-clock Hamming distance')
 require(h['selected']['cost']==8-len(h['selected']['common_clocks']),'all equal clocks pinned')
# Four normalized endpoints separately counted as words, one orbit as one resource.
FULL=(1<<21)-1;chords=[(i,j) for i in range(8) for j in range(i+2,8)]
rev=[chords.index((7-j,7-i)) for i,j in chords]
def reflect(m):return sum(((m>>k)&1)<<rev[k] for k in range(21))
def orbit(m):return {m,reflect(m),FULL^m,FULL^reflect(m)}
allwords=set();res=set();costs=[]
for h in helpers.values():
 costs.append(h['selected']['cost'])
 for name,target in [('before',1965),('after',2234)]:
  r=h['selected'][name];o=orbit(r['mask'])
  require(r['target']==target and len(o)==4 and sorted(o)==r['orbit'] and min(o)==r['resource'],'literal normalized closure')
  require(not (allwords&o),'normalized endpoint resource collision');allwords|=o;res.add(min(o))
require(len(res)==106 and len(allwords)==424,'helper reservation count')
old=json.loads((RESULTS/C['old_domain']).read_text())['sources'];oldwords=set();oldres=set();oldcosts=[]
for h in old:
 oldcosts.append(sum(a!=b for a,b in zip(h['old_before'],h['old_after'])))
 for r in h['exact_all_shared_clocks_target_domain']:
  o=orbit(r['mask']);require(not(o&allwords),'helper vs old normalized collision');oldwords|=o;oldres.add(min(o))
require(len(old)==28 and len(oldres)==56 and len(oldwords)==224,'old native resource counts')
# Exhaust every distinct triple from all53 actual slots, using ORIGINAL labels.
triples=[[] for _ in range(2048)]
for cc in itertools.combinations(range(53),3):
 q=qsum(labels[ids[i]] for i in cc);triples[q].append(sum(1<<i for i in cc))
require(sum(map(len,triples))==23426 and all(triples),'full exact triple census')
fullmask=(1<<53)-1
single_failures=[]
for q,fs in enumerate(triples):
 intersection=fullmask
 for m in fs:intersection &= m
 if intersection:single_failures.append((q,intersection))
require(not single_failures,'not all one-U restrictions supported')
expected=[(776,1182,2638),(1391,1948,3139),(1400,3111,3139),(2439,2728,3139)]
actual=[tuple(ids[i] for i in range(53) if m>>i&1) for m in triples[305]]
require(actual==expected,'charge305 exact fiber')
bad=(1<<ids.index(776))|(1<<ids.index(3139))
require(all(m&bad for m in triples[305]),'two-U hostile does not hit all305 triples')
# Main40: all single restrictions really fail and miss counts agree with inherited claims.
mainids=C['banks'][C['main_bank_index']]['source_ids'];main_del=[]
for s in mainids:
 image={qsum(labels[v] for v in cc) for cc in itertools.combinations([t for t in mainids if t!=s],3)}
 missing=sorted(set(range(2048))-image);require(bool(missing),'40 main bank is not inclusion-minimal')
 main_del.append({'source':s,'missing_count':len(missing),'first_missing':missing[0]})
# Independent 46-bank verification, no inherited optimizer Boolean accepted.
B46=json.loads((HERE/'robust_bank46.json').read_text());ids46=B46['source_ids']
require(B46['count']==len(ids46)==46 and len(set(ids46))==46 and set(ids46)<=set(ids),'46 distinct authenticated slots')
def fiber_for(bank):
 fs=[[] for _ in range(2048)]
 for cc in itertools.combinations(bank,3):fs[qsum(labels[s] for s in cc)].append(frozenset(cc))
 return fs
F46=fiber_for(ids46)
require(all(F46),'46-bank coverage without restriction')
for q,fs in enumerate(F46):require(not set.intersection(*(set(f) for f in fs)),'46-bank one-U robustness')
controls46=[]
for removed in ids46:
 rem=[s for s in ids46 if s!=removed]
 fs=fiber_for(rem);found=None
 for q,f in enumerate(fs):
  if not f:found={'first_removed':removed,'charge':q,'already_missing':True};break
  inter=set.intersection(*(set(m) for m in f))
  if inter:found={'first_removed':removed,'charge':q,'second_forced_U':min(inter)};break
 require(found is not None,'46 bank not inclusion-minimal for one-U robustness')
 controls46.append(found)
require(all({776,3139}&f for f in F46[305]),'46-bank exact two-U failure')
# Check every supplied four-subset packing against original owner labels, not producer arithmetic.
P=json.loads(P53.read_text());require(P['source_ids']==ids,'packing source IDs mismatch')
require(len(P['packings'])==2048,'must certify every original charge')
packing_counts=[]
for q,ms in enumerate(P['packings']):
 used=0
 for m in ms:
  require(type(m) is int and 0<=m<1<<53 and m.bit_count()==4,'quad mask type/domain/cardinality')
  require(not(m&used),'quad packing not pairwise disjoint');used|=m
  require(qsum(labels[ids[i]] for i in range(53) if m>>i&1)==q,'quad actual-original charge mismatch')
 packing_counts.append(len(ms))
require(min(packing_counts)>=9,'advertised nine-disjoint universal certificate fails')
# Explicit mixed restriction: force776 to V, force3139 to U; others remain eligible.
fv=776;fu=3139;target=add(305,labels[fv]);remaining=[s for s in ids if s not in (fv,fu)]
solutions=[cc for cc in itertools.combinations(remaining,3) if qsum(labels[s] for s in cc)==305]
require(not solutions,'mixed1V1U hostile unexpectedly feasible')
# Independently census quadruple multiplicities; repeated labels never replace actual distinct slots.
quads=[0]*2048
for cc in itertools.combinations(ids,4):quads[qsum(labels[s] for s in cc)]+=1
require(sum(quads)==292825 and min(quads)>0,'quad full census fails')
out={'status':'PASS','checks':checks,'input_sha256':sha(SUPPLY),'packing_sha256':sha(P53),
 'original_labels_sha256':sha(RESULTS/'metagraph_order8_electrical_20261023/receipt.json'),
 'triple_total':sum(map(len,triples)),'triple_multiplicity_range':[min(map(len,triples)),max(map(len,triples))],
 'every_one_U_supported':True,'two_U_failure':{'charge':305,'force_U':[776,3139],'complete_original_triples':actual},
 'main40_single_deletion_controls':main_del,
 'bank46':{'source_sha256':sha(HERE/'robust_bank46.json'),'sources':ids46,'triple_count':sum(map(len,F46)),'triple_multiplicity_range':[min(map(len,F46)),max(map(len,F46))],'every_one_U_supported':True,'inclusion_minimal_for_one_U_robustness':True,'inclusion_minimal_controls':controls46,'two_U_failure_charge':305,'two_U_failure':[776,3139],'core_plus_helper_arrival_slots':74,'reserved_resources':148,'reserved_words':592,'arrival_words':296},
 'quad_total':sum(quads),'quad_multiplicity_range':[min(quads),max(quads)],
 'verified_disjoint_quad_packing_histogram':dict(sorted(Counter(packing_counts).items())),
 'mixed_failure':{'force_V':[fv],'force_U':[fu],'target_charge':target,'residual_charge':305},
 'resources':{'helper_slots':53,'helper_normalized_reservation':106,'helper_closed_words':424,'core_slots':28,'total_slots_per_arrival':81,'total_resources_reserved':162,'total_words_reserved':648,'total_words_per_arrival':324},
 'costs':{'helper_distribution':dict(sorted(Counter(costs).items())),'helper_largest_six_sum':sum(sorted(costs,reverse=True)[:6]),'helper_largest_eight_sum':sum(sorted(costs,reverse=True)[:8]),'core_largest_two_sum':sum(sorted(oldcosts,reverse=True)[:2]),'three_occupancy_support_bound':8,'three_occupancy_clock_bound':sum(sorted(costs,reverse=True)[:6])+sum(sorted(oldcosts,reverse=True)[:2]),'four_occupancy_support_bound':10,'four_occupancy_clock_bound':sum(sorted(costs,reverse=True)[:8])+sum(sorted(oldcosts,reverse=True)[:2])}}
out=json.loads(json.dumps(out))
if args.write:dump(HERE/'receipt.json',out)
else:require(json.loads((HERE/'receipt.json').read_text())==out,'frozen independent audit receipt mismatch')
print(json.dumps(out,indent=2,sort_keys=True))
