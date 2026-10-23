#!/usr/bin/env python3
"""Independent fixed-pool upper bound: enumerate all 9-subsets, not omitted pairs."""
import argparse,json,hashlib,itertools
from pathlib import Path
from collections import defaultdict,Counter
P=argparse.ArgumentParser();P.add_argument('--repo',type=Path,required=True);P.add_argument('--packet',type=Path,required=True);P.add_argument('--out',type=Path,required=True);a=P.parse_args()
checks=0
def require(v,why):
 global checks
 checks+=1
 if not v:raise RuntimeError(why)
Epath=a.repo/'05-knowledge/results/metagraph_order8_electrical_20261023/receipt.json';Ppath=a.packet/'pool.json';E=json.loads(Epath.read_text())['order8'];data=json.loads(Ppath.read_text());X=E['higher_layers'][-1]['lifts'][0];Y=E['radical_lifts_mod4']
ids=sorted(row['source']for row in data['options']);require(len(ids)==len(set(ids))==20,'distinct fixed menu source count')
labels={s:(sum((Y[i][s]%2)*(1<<i)for i in range(4)),X[s]%128)for s in ids}
for row in data['options']:
 s=row['source'];require(labels[s]==tuple(data['labels'][str(s)]),'pool compressed charge labels match literal original electrical data');require(row['source_label']==[Y[i][s]%2 for i in range(4)]+[X[s]%128],'full source labels match literal originals')
def add(subset):
 parity=0;mass=0
 for j in subset:parity^=labels[ids[j]][0];mass+=labels[ids[j]][1]
 return parity,mass%128
fivebit=[labels[s][0]|((labels[s][1]%2)<<4)for s in ids];totalparity=0
for v in fivebit:totalparity^=v
require(totalparity==6,'twenty-source total five-bit parity forbids every 10+10 partition')
groups=defaultdict(list);subsetcount=0
for subset in itertools.combinations(range(20),9):
 mask=sum(1<<i for i in subset);groups[add(subset)].append(mask);subsetcount+=1
require(subsetcount==167960,'all nine-subsets independently enumerated')
relations=[];pairtests=0;omitted=Counter();proper_sizes=Counter();full=(1<<20)-1
for charge,group in groups.items():
 for left,right in itertools.combinations(group,2):
  pairtests+=1
  if left&right:continue
  A=[i for i in range(20)if left>>i&1];B=[i for i in range(20)if right>>i&1]
  require(len(A)==len(B)==9 and add(A)==add(B),'actual disjoint nine-nine charge relation')
  # Search every size1..4 proper balanced subrelation independently. Since
  # the complete relation is balanced, complements cover possible sizes5..8.
  witness=None
  for k in range(1,5):
   index=defaultdict(list)
   for C in itertools.combinations(B,k):index[add(C)].append(C)
   for C in itertools.combinations(A,k):
    if add(C)in index:witness=(C,index[add(C)][0]);break
   if witness:break
  require(witness is not None,'every 18-support relation has a proper balanced subrelation')
  C,D=witness;require(0<len(C)==len(D)<9 and add(C)==add(D),'explicit proper balanced witness rechecked')
  om=tuple(ids[i]for i in range(20)if not(left|right)>>i&1);omitted[om]+=1;proper_sizes[len(C)]+=1
  relations.append({'plus':[ids[i]for i in A],'minus':[ids[i]for i in B],'proper_plus':[ids[i]for i in C],'proper_minus':[ids[i]for i in D]})
require(len(relations)==1129 and len(omitted)==3,'independent exhaustive 18-support relation count')
require(dict(proper_sizes)=={2:288,3:832,4:9},'minimum-side proper-subrelation histogram')
# Authenticate the exhibited lower bound directly from its labels, without
# depending on another search or receipt.
cases=json.loads((a.packet/'candidates.json').read_text())['found'];c=cases['8'];A=[ids.index(s)for s in c['plus']];B=[ids.index(s)for s in c['minus']]
require(len(A)==len(B)==8 and not set(A)&set(B)and add(A)==add(B),'exhibited sixteen-support relation')
for k in range(1,5):
 values={add(C)for C in itertools.combinations(A,k)}
 require(all(add(C)not in values for C in itertools.combinations(B,k)),'exhibited sixteen-support relation has no proper balanced subrelation')
report={'status':'PASS','checks':checks,'scope':'Exact maximum support 16 among balanced primitive one-option-per-source circuits in this fixed 20-source menu; not a global order-eight or all-order bound. Native authentication is the separate main audit.','algorithm':'All 167,960 nine-subsets grouped by actual original charge, then all unordered same-charge pairs tested for disjointness, without omitted-pair or five-bit pruning.','nine_subsets':subsetcount,'occupied_charge_groups':len(groups),'same_charge_pairs_tested':pairtests,'eighteen_relations':len(relations),'omitted_pair_counts':[{'omitted':list(k),'relations':v}for k,v in sorted(omitted.items())],'proper_subrelation_side_size_counts':dict(proper_sizes),'twenty_source_fivebit_parity':totalparity,'sharp_maximum_support':16,'witness_digest_sha256':hashlib.sha256(json.dumps(sorted(relations,key=lambda r:r['plus']),sort_keys=True,separators=(',',':')).encode()).hexdigest(),'source_hashes':{'electrical':hashlib.sha256(Epath.read_bytes()).hexdigest(),'pool':hashlib.sha256(Ppath.read_bytes()).hexdigest()}}
a.out.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');print(json.dumps(report,sort_keys=True))
