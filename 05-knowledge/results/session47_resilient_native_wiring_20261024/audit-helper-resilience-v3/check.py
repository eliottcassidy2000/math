#!/usr/bin/env python3
"""Independent Session47 finite helper resilience audit; no assertion statements."""
import argparse,csv,gzip,hashlib,itertools,json
from pathlib import Path
from collections import Counter
from fractions import Fraction
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
parser=argparse.ArgumentParser();parser.add_argument('--repo-root',type=Path,default=ROOT/'session35-publication/clean-repo');parser.add_argument('--write',action='store_true');args=parser.parse_args()
RESULTS=args.repo_root/'05-knowledge/results'
SUPPLY=RESULTS/'compatible_charged_helpers_20261023/supply/selected_supply.json'
P53=HERE/'quad_packings12.json'
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
require(min(packing_counts)>=12,'advertised twelve-disjoint universal certificate fails')
# Explicit mixed restriction: force776 to V, force3139 to U; others remain eligible.
fv=776;fu=3139;target=add(305,labels[fv]);remaining=[s for s in ids if s not in (fv,fu)]
solutions=[cc for cc in itertools.combinations(remaining,3) if qsum(labels[s] for s in cc)==305]
require(not solutions,'mixed1V1U hostile unexpectedly feasible')
# Independently census quadruple multiplicities; repeated labels never replace actual distinct slots.
quads=[0]*2048;quadfibers=[[] for _ in range(2048)]
for cc in itertools.combinations(ids,4):
 q=qsum(labels[s] for s in cc);quads[q]+=1;quadfibers[q].append(frozenset(cc))
require(sum(quads)==292825 and min(quads)>0,'quad full census fails')
# All mixed restrictions, independently reduced to complete residual triple fibers.
mixed_failures=[];mixed_pair_counts=Counter()
for vi,v in enumerate(ids):
 for q,fs in enumerate(triples):
  inter=(1<<53)-1;alive=False
  for m in fs:
   if not(m>>vi&1):inter &= m;alive=True
  require(alive,'one-V occupancy4 positive padding claim')
  for ui,u in enumerate(ids):
   if inter>>ui&1:
    require(u!=v,'mixed restriction pairs disjoint')
    mixed_pair_counts[v,u]+=1
    mixed_failures.append({'force_V':v,'force_U':u,'residual_charge':q,'total_charge':add(q,labels[v])})
require(len(mixed_pair_counts)==32 and len(mixed_failures)==32,'mixed pair census32 failures')
require(set(mixed_pair_counts.values())=={1},'each hostile directed pair misses exactlyone charge')
# A 13-packing uses52 owners, so13q=T-label(omitted). Inverse13 modulo128 is69.
T=qsum(labels.values())
require(len(set(labels.values()))==53,'53 distinct owner labels for obstruction count')
eligible={((((T//16-l//16)*69)%128)*16)|((T^l)&15) for l in labels.values()}
require(len(eligible)==53 and len(set(range(2048))-eligible)==1995,'uniform13 packing obstruction')
for q in eligible:
 thirteen_q=16*((13*(q//16))%128)|(q&15)
 require(any(add(thirteen_q,l)==T for l in labels.values()),'explicit13-packing necessary equation')
require(min(packing_counts)==12,'sharp uniform packing size12')
# Independently validate48 bank too.
B48=json.loads((HERE/'robust_bank48.json').read_text());ids48=B48['source_ids'];F48=fiber_for(ids48)
require(len(ids48)==48 and len(set(ids48))==48 and set(ids48)<=set(ids),'48 native owner set')
require(all(F48),'48-bank all charges')
for f in F48:require(not set.intersection(*(set(t) for t in f)),'48-bank one-U robustness')
# Direct sparse signed-Laplacian replay, not producer's accumulated-row routine.
core=sorted(h['source'] for h in old);cost_by={s:h['selected']['cost'] for s,h in helpers.items()}
for h in old:cost_by[h['source']]=sum(a!=b for a,b in zip(h['old_before'],h['old_after']))
all_sources=set(core)|set(ids)
require(not(set(core)&set(ids)),'native core and helper indices disjoint')
core_labels={s:16*(X[s]%128)+sum((Y[j][s]%2)<<(j-1) for j in range(1,5)) for s in core}
all_labels={**labels,**core_labels}
def charge_native(ss):return qsum(all_labels[s] for s in ss)
def subtract(a,b):return (16*((a//16-b//16)%128))|((a^b)&15)
pots=[X]+Y[:5];moduli=[256]+[4]*5
start=set(core[:14]);finish=set(core[14:]);moves=list(zip(sorted(start),sorted(finish)));moves += [(b,a) for a,b in reversed(moves)]
route_events=0;evacuations=0;peak_support=0;peak_clocks=0
for total in range(2048):
 initial_helpers=F46[subtract(total,charge_native(start))][0];curcore=set(start);curhelpers=initial_helpers
 for step,(a,b) in enumerate(moves):
  nxtcore=(curcore-{a})|{b};debt=subtract(total,charge_native(nxtcore))
  if step==27:nexthelpers=initial_helpers;forbidden=next(s for s in ids46 if s not in initial_helpers)
  else:
   forbidden=min(curhelpers);nexthelpers=next(t for t in F46[debt] if forbidden not in t);evacuations+=1
  before=curcore|set(curhelpers);after=nxtcore|set(nexthelpers)
  require(forbidden not in nexthelpers and len(nexthelpers)==3 and len(nxtcore)==14,'causal destination restriction and margins')
  require(charge_native(after)==total,'total exact marked charge in independently selected arrival')
  signed_edges=[]
  for s in sorted(before^after):
   delta=1 if s in after else -1
   signed_edges.extend([(s,1965,-2*delta),(s,2234,2*delta)])
  rows={}
  for s,t,w in signed_edges:
   for j,phi in enumerate(pots):
    d=w*(phi[s]-phi[t]);rows[s,j]=rows.get((s,j),0)+d;rows[t,j]=rows.get((t,j),0)-d
  require(all(val%moduli[j]==0 for (s,j),val in rows.items()),'direct signed edge Laplacian all physical changing rows')
  support=len(before^after);clocks=sum(cost_by[s] for s in before^after)
  require(support<=8 and clocks<=50,'independently replayed primary slot/clock budgets')
  peak_support=max(peak_support,support);peak_clocks=max(peak_clocks,clocks);route_events+=1;curcore=nxtcore;curhelpers=nexthelpers
 require(curcore==start and curhelpers==initial_helpers,'full literal endpoint return at fixed records')
require(route_events==57344 and evacuations==55296,'all2048 totals and28 exchanges')
# Exact native fractional-matching certificates: all arithmetic here is integer/rational.
W=json.loads((HERE/'fractional_weights.json').read_text())
require(W['source_ids']==ids and len(W['charges'])==2048,'weighted certificate full original owner universe')
weight_count=0;weighted_stats=[];ratios=[]
for q,certificate in enumerate(W['charges']):
 d=certificate['denominator'];require(type(d) is int and d>0,'exact positive integer vertex capacity')
 loads={s:0 for s in ids};total=0;seen=set()
 for mask,weight in certificate['edges']:
  require(type(mask) is int and 0<=mask<1<<53 and mask.bit_count()==4,'weighted edge four distinct native owner indices')
  require(mask not in seen,'no duplicate weighted edge');seen.add(mask)
  require(type(weight) is int and weight>0,'strictly positive integer edge weight')
  edge=[ids[i] for i in range(53) if mask>>i&1]
  require(qsum(labels[s] for s in edge)==q,'each weighted quadruple original physical owner label')
  for vertex in edge:loads[vertex]+=weight
  total+=weight;weight_count+=1
 require(all(load<=d for load in loads.values()),'integer weighted load bound for each actual owner')
 require(total>13*d,'weighted matching strictly exceeds thirteen')
 require(sum(loads.values())==4*total,'exact four-uniform incidence sum')
 require(4*total<=53*d,'fractional method ceiling53/4')
 ratios.append(Fraction(total,d))
 weighted_stats.append({'charge':q,'edge_count':len(seen),'total_integer_weight':total,'denominator':d,'maximum_integer_vertex_load':max(loads.values()),'minimum_integer_vertex_load':min(loads.values()),'strict_integer_margin_over13':total-13*d})
minq=min(range(2048),key=lambda q:ratios[q]);minratio=ratios[minq]
require(weight_count==108309,'complete108309 weighted native quadruples')
require(minq==540 and minratio==Fraction(13145902,1000001),'minimum exact weighted ratio')
upper=json.loads((HERE/'upper_transversal.json').read_text());bad=upper['forced_U'];badset=set(bad)
require(len(bad)==len(badset)==upper['size']==15 and badset<=set(ids),'15 distinct actual forced-U owners')
require(upper['charge']==353 and len(quadfibers[353])==upper['full_fibre_count']==139,'independent complete139quad fiber353')
require(all(edge&badset for edge in quadfibers[353]),'15-set hits every actual353 quadruple')
# Independently replay both directions of all2048 occupancy4 two-event returns.

quadsets=[[frozenset(ids[i] for i in range(53) if m>>i&1) for m in ms] for ms in P['packings']]
quad_events=0;quad_peak_support=0;quad_peak_clocks=0
nextcore=(start-{core[0]})|{core[14]}
for total in range(2048):
 h0=quadsets[subtract(total,charge_native(start))][0]
 D0=set(h0)
 for owner in ids:
  if len(D0)==13:break
  D0.add(owner)
 h1=next(edge for edge in quadfibers[subtract(total,charge_native(nextcore))] if not(edge&D0))
 D1=set(sorted(set(ids)-set(h0))[:13])
 require(len(D0)==len(D1)==13 and set(h0)<=D0,'thirteen forced-U owners include all four current V owners')
 before=start|set(h0)
 for corestate,helperstate,D in [(nextcore,h1,D0),(start,h0,D1)]:
  require(not(set(helperstate)&D) and len(corestate)==14 and len(helperstate)==4,'occupancy4 admissible arrival and margins')
  after=corestate|set(helperstate);require(charge_native(after)==total,'occupancy4 exact marked total')
  signed_edges=[]
  for s in sorted(before^after):
   delta=1 if s in after else -1;signed_edges.extend([(s,1965,-2*delta),(s,2234,2*delta)])
  rows={}
  for a,b,w in signed_edges:
   for j,phi in enumerate(pots):
    d=w*(phi[a]-phi[b]);rows[a,j]=rows.get((a,j),0)+d;rows[b,j]=rows.get((b,j),0)-d
  require(all(val%moduli[j]==0 for (s,j),val in rows.items()),'occupancy4 both direction signed Laplacian current test')
  support=len(before^after);clocks=sum(cost_by[s] for s in before^after)
  require(support<=10 and clocks<=64,'occupancy4 both direction support/work bounds')
  quad_events+=1;quad_peak_support=max(quad_peak_support,support);quad_peak_clocks=max(quad_peak_clocks,clocks);before=after
 require(before==start|set(h0),'occupancy4 exact literal full record return')
require(quad_events==4096,'all2048 occupancy4 routes both directions')
out={'status':'PASS','checks':checks,'input_sha256':sha(SUPPLY),'packing_sha256':sha(P53),
 'original_labels_sha256':sha(RESULTS/'metagraph_order8_electrical_20261023/receipt.json'),
 'fractional_certificate':{'weighted_edges':weight_count,'minimum_ratio_charge':minq,'minimum_ratio_numerator':minratio.numerator,'minimum_ratio_denominator':minratio.denominator,'guaranteed_minimum_transversal':14,'certificate_sha256':sha(HERE/'fractional_weights.json'),'per_charge_exact_statistics':weighted_stats},
 'upper_transversal':{'certificate':upper,'certificate_sha256':sha(HERE/'upper_transversal.json'),'complete_quadruple_fiber':[sorted(t) for t in quadfibers[353]]},
 'occupancy4_deletion_radius_interval':[13,14],
 'sharp_uniform_quad_packing_size':12,'guaranteed_occupancy4_U_restrictions':13,
 'thirteen_packing_necessary_equation':{'all53_charge':T,'eligible_charges':sorted(eligible),'ineligible_charge_count':1995},
 'mixed_census':{'directed_pairs':53*52,'complete_pairs':53*52-len(mixed_pair_counts),'hostile_pairs':len(mixed_pair_counts),'failures':mixed_failures},
 'bank48':{'triple_count':sum(map(len,F48)),'multiplicity_range':[min(map(len,F48)),max(map(len,F48))],'every_one_U_supported':True,'sha256':sha(HERE/'robust_bank48.json')},
 'independent_occupancy4_route_replay':{'total_charges':2048,'events':quad_events,'forced_U_each_destination':13,'peak_support':quad_peak_support,'peak_primary_clocks':quad_peak_clocks,'literal_endpoint_return':True},
 'independent_causal_route_replay':{'total_charges':2048,'events':route_events,'actual_forced_evacuations':evacuations,'peak_support':peak_support,'peak_primary_clocks':peak_clocks,'direct_sparse_signed_Laplacian_test':True,'literal_endpoint_return':True},
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
