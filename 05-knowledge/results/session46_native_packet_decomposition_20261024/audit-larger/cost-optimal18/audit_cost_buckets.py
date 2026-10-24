#!/usr/bin/env python3
"""Third-strategy exact cost-boundary relation enumeration, no support parity prune.
Group every cheap nine-subset by its complete original electrical charge.
For cost36 test every unordered same-charge pair for disjointness.
For cost37 test every expensive-plus-eight-cheap subset against its entire
matching cheap-nine charge bucket. This never enumerates or prunes supports.
"""
import argparse,hashlib,itertools,json
from collections import defaultdict
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,required=True);p.add_argument('--packet',type=Path,required=True);p.add_argument('--census',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();checks=0

def need(ok,msg):
 global checks
 checks+=1
 if not ok:raise RuntimeError(msg)
electrical=a.repo/'05-knowledge/results/metagraph_order8_electrical_20261023/receipt.json';E=json.loads(electrical.read_text())['order8'];Y=E['radical_lifts_mod4'];X=E['higher_layers'][-1]['lifts'][0]
rows=json.loads((a.packet/'pool.json').read_text())['options'];cost={r['source']:sum(x!=y for x,y in zip(r['before']['path'],r['after']['path']))for r in rows};idx={s:i for i,s in enumerate(sorted(cost))};cheap=sorted(s for s,c in cost.items()if c==2);extra=sorted(s for s,c in cost.items()if c==3)
label={s:(sum((Y[h][s]%2)<<h for h in range(4)),X[s]%128)for s in cost}
need(len(cheap)==20 and len(extra)==5,'full cheap and cost-three strata')
def code(S):
 p=x=0
 for s in S:p^=label[s][0];x+=label[s][1]
 return (p<<7)|(x%128)
buckets=defaultdict(list)
for C in itertools.combinations(cheap,9):buckets[code(C)].append(sum(1<<idx[s]for s in C))
need(sum(map(len,buckets.values()))==167960,'all cheap nine-subsets enumerated')
results=[];previous={r['cost']:r for r in json.loads(a.census.read_text())['results']}
for target in (36,37):
 relations=[];examined=0;positive_subsets=0
 if target==36:
  for group in buckets.values():
   for A,B in itertools.combinations(group,2):
    examined+=1
    if not A&B:relations.append(sorted([A,B]))
 else:
  for C in itertools.combinations(cheap,8):
   cp,cx=code(C)>>7,code(C)&127;cm=sum(1<<idx[s]for s in C)
   for e in extra:
    positive_subsets+=1;ep,ex=label[e];A=cm|(1<<idx[e]);q=((cp^ep)<<7)|((cx+ex)%128)
    for B in buckets[q]:
     examined+=1
     if not A&B:relations.append(sorted([A,B]))
 need(len({tuple(r)for r in relations})==len(relations),'each unordered signed relation occurs exactly once')
 sha=hashlib.sha256(json.dumps(sorted(relations),separators=(',',':')).encode()).hexdigest()
 need(len(relations)==previous[target]['equal_charge_relations']and sha==previous[target]['complete_relation_set_sha256'],'exact relation-set equality with independent support/partition census')
 results.append({'cost':target,'matching_charge_pairs_examined':examined,'expensive_plus_eight_subsets':positive_subsets,'disjoint_relations':len(relations),'complete_relation_set_sha256':sha})
report={'status':'PASS','checks':checks,'source_sha256':{'original_electrical':hashlib.sha256(electrical.read_bytes()).hexdigest(),'pool.json':hashlib.sha256((a.packet/'pool.json').read_bytes()).hexdigest(),'census':hashlib.sha256(a.census.read_bytes()).hexdigest()},'cheap_nine_subsets':167960,'occupied_charge_buckets':len(buckets),'results':results,'scope':'Independent exhaustive relation-set enumeration with no support-parity pruning. Uniqueness and primitiveness use the separately verified per-relation proper-subset census, whose complete relation set matches exactly.'}
a.out.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');print(json.dumps(report,sort_keys=True))
