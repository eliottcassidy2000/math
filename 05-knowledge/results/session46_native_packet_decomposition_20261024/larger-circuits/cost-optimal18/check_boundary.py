#!/usr/bin/env python3
"""Exhaustive clock-cost36/37 boundary, independent of the C++ discovery loop."""
from pathlib import Path
from itertools import combinations
from collections import Counter
import json
H=Path(__file__).resolve().parent
P=json.loads((H/'pool.json').read_text());L={int(s):tuple(q) for s,q in P['labels'].items()};C={r['source']:r['cost'] for r in P['options']}
checks=0
def need(ok,why):
 global checks
 checks+=1
 if not ok:raise RuntimeError(why)
need(len(C)==28 and len(set(L.values()))==28,'exact distinct-label full28 supplied options')
need(Counter(C.values())=={2:20,3:5,4:3},'exact actual clock-cost inventory')
cheap=sorted(s for s,c in C.items() if c==2);three=sorted(s for s,c in C.items() if c==3)
def charge(seq):
 p=x=0
 for s in seq:p^=L[s][0];x+=L[s][1]
 return p,x%128
def first_proper(A,B):
 for j in range(2,5):
  bank={charge(z) for z in combinations(A,j)}
  for z in combinations(B,j):
   if charge(z) in bank:return j
 return None
reports=[]
for cost in [36,37]:
 supports=(tuple(z) for z in combinations(cheap,18)) if cost==36 else (tuple(sorted((*z,s))) for s in three for z in combinations(cheap,17))
 counts=Counter();primitives=[];witness_sizes=Counter()
 for S in supports:
  counts['support_sets']+=1;p,total=charge(S)
  if p or total%2:continue
  counts['parity_admissible_support_sets']+=1
  for tail in combinations(S[1:],8):
   counts['balanced_partitions_examined']+=1
   A=(S[0],*tail)
   if (2*sum(L[s][1] for s in A)-total)%128:continue
   counts['equal_charge_relations']+=1;aset=set(A);B=tuple(s for s in S if s not in aset)
   need(charge(A)==charge(B),'candidate balanced charges agree exactly')
   bad=first_proper(A,B)
   if bad is None:primitives.append({'plus':list(A),'minus':list(B)})
   else:witness_sizes[bad]+=1
 saved=json.loads((H/f'cost{cost}.json').read_text())
 for name,value in counts.items():need(saved[name]==value,'independent boundary count '+str(cost)+' '+name)
 need(saved['primitive_relations']==len(primitives),'independent primitive count')
 if primitives:need(primitives[0]=={'plus':saved['first_plus'],'minus':saved['first_minus']},'independent unique primitive support')
 reports.append({'cost':cost,**counts,'primitive_relations':len(primitives),'proper_relation_half_sizes':dict(witness_sizes),'primitives':primitives})
need(reports[0]['primitive_relations']==0 and reports[1]['primitive_relations']==1,'sharp minimum37 and uniqueness in supplied28options')
print(json.dumps({'status':'PASS','checks':checks,'minimum_18_slot_clock_cost_in_supplied_28_options':37,'reports':reports,'scope':'One of each supplied fixed native binary option per distinct source; not all order8 pinned menus.'},indent=2,sort_keys=True))
