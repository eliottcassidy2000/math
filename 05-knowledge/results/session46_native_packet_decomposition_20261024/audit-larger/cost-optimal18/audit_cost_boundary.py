#!/usr/bin/env python3
"""Independent exhaustive cost36/37 support audit, no producer-code import.
Supports are enumerated by the omitted cheap-source set and (at cost37) the
unique cost-three source. Partitions anchor the expensive source on the plus
side, unlike the producer's smallest-source anchoring. All proper balanced
subrelations are tested through side-size four, using original graph labels.
"""
import argparse,csv,hashlib,itertools,json,math
from pathlib import Path
from collections import Counter
p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,required=True);p.add_argument('--packet',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();checks=0

def need(ok,msg):
 global checks
 checks+=1
 if not ok:raise RuntimeError(msg)
D=a.repo/'05-knowledge/results';electrical=D/'metagraph_order8_electrical_20261023/receipt.json';E=json.loads(electrical.read_text())['order8'];Y=E['radical_lifts_mod4'];X=E['higher_layers'][-1]['lifts'][0]
records=json.loads((a.packet/'pool.json').read_text())['options'];cost={r['source']:sum(x!=y for x,y in zip(r['before']['path'],r['after']['path']))for r in records}
need(all(cost[r['source']]==r['cost']for r in records),'cost independently derived from literal owner-clock changes')
label={s:(sum((Y[h][s]%2)<<h for h in range(4)),X[s]%128)for s in cost}
need(len(set(label.values()))==len(label)==28,'distinct original source charge labels eliminate side-size1 relations')
cheap=sorted(s for s,c in cost.items()if c==2);extra=sorted(s for s,c in cost.items()if c==3)
need(len(cheap)==20 and len(extra)==5 and Counter(cost.values())=={2:20,3:5,4:3},'full cost strata')

def code(S):
 par=0;mass=0
 for s in S:par^=label[s][0];mass+=label[s][1]
 return (par,mass%128)
# Cache small subset charges keyed by exact owner set. This cache is generated
# independently of the producer's per-partition stamped-array routine.
small={}
for j in range(1,5):
 for comb in itertools.combinations(cheap+extra,j):small[frozenset(comb)]=code(comb)

def obstruction(A,B):
 for j in range(1,5):
  positive={small[frozenset(t)]:t for t in itertools.combinations(A,j)}
  for t in itertools.combinations(B,j):
   q=small[frozenset(t)]
   if q in positive:return j,positive[q],t
 return None
results=[];witness_hash=hashlib.sha256();owner_index={s:j for j,s in enumerate(sorted(cost))}
for target in (36,37):
 supports=[]
 if target==36:
  for omitted in itertools.combinations(cheap,2):supports.append(tuple(s for s in cheap if s not in omitted))
 else:
  for e in extra:
   for omitted in itertools.combinations(cheap,3):supports.append((e,)+tuple(s for s in cheap if s not in omitted))
 relations=0;primitive=[];admissible=0;partitions=0;minimum=Counter();checked_supports=0;canonical_relations=[]
 for support in supports:
  need(len(support)==len(set(support))==18 and sum(cost[s]for s in support)==target,'distinct literal support and exact cost')
  checked_supports+=1
  par,mass=code(support)
  if par or mass%2:continue
  admissible+=1;anchor=support[0];remaining=support[1:];support_set=set(support)
  for c in itertools.combinations(remaining,8):
   A=(anchor,)+c;partitions+=1
   # Full group equality is directly tested, rather than assuming parity after a prune.
   B=tuple(sorted(support_set-set(A)))
   if code(A)!=code(B):continue
   relations+=1;canonical_relations.append(sorted([sum(1<<owner_index[s]for s in A),sum(1<<owner_index[s]for s in B)]));bad=obstruction(A,B)
   if bad is None:
    primitive.append(sorted([sorted(A),sorted(B)]))
   else:
    j,pa,pb=bad;need(len(pa)==len(pb)==j and code(pa)==code(pb)and set(pa)<=set(A)and set(pb)<=set(B),'exact proper balanced witness')
    minimum[j]+=1;witness_hash.update(json.dumps([sorted(A),sorted(B),j,sorted(pa),sorted(pb)],separators=(',',':')).encode())
 need(checked_supports==(math.comb(20,18)if target==36 else 5*math.comb(20,17)),'complete cheap-plus-one-expensive support count')
 need(partitions==admissible*math.comb(17,8),'complete unordered sign partition count')
 result={'cost':target,'support_sets':checked_supports,'parity_admissible_support_sets':admissible,'balanced_partitions_examined':partitions,'equal_charge_relations':relations,'primitive_relations':len(primitive),'canonical_primitive_partitions':primitive,'minimal_proper_subrelation_side_size_counts':dict(minimum),'complete_relation_set_sha256':hashlib.sha256(json.dumps(sorted(canonical_relations),separators=(',',':')).encode()).hexdigest()}
 results.append(result)
need(results[0]['primitive_relations']==0 and results[1]['primitive_relations']==1,'sharp attainable electrical support18 cost boundary')
claimed=json.loads((a.packet/'cost37.json').read_text());actual=results[1]
for k in ['cost','support_sets','parity_admissible_support_sets','balanced_partitions_examined','equal_charge_relations','primitive_relations']:need(actual[k]==claimed[k],'producer exact-cost count '+k)
need(actual['canonical_primitive_partitions']==[sorted([sorted(claimed['first_plus']),sorted(claimed['first_minus'])])],'producer unique primitive support and sign partition')
report={'status':'PASS','checks':checks,'input_sha256':{'pool.json':hashlib.sha256((a.packet/'pool.json').read_bytes()).hexdigest(),'cost37.json':hashlib.sha256((a.packet/'cost37.json').read_bytes()).hexdigest(),'original_electrical':hashlib.sha256(electrical.read_bytes()).hexdigest()},'results':results,'proper_subrelation_witness_stream_sha256':witness_hash.hexdigest(),'scope':'Exact minimum ordinary owner-clock cost37 for an18-source primitive equal-charge binary switch in this fixed28 one-slot-per-source supply, once its independent native and full-image endpoint certificate is supplied. Cost<36 impossible since each nonidentity option changes at least2 clocks. No other menus, conventions, supports or all-order claim.'}
a.out.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');print(json.dumps({'status':'PASS','checks':checks,'results':results},sort_keys=True))
