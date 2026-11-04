"""Exact rooted path comparison for the pinned incoming four-arc route."""
import itertools,json
from pathlib import Path
HERE=Path(__file__).resolve().parent
# Publication uses spectrum/spectrum; workspace uses session190-spectrum.
producer=HERE/'spectrum'
if not producer.exists(): producer=HERE.parent/'session190-spectrum'
ours=json.loads((producer/'receipt.json').read_text())['patch']
incoming=[8,73,3,4,15,31,61]
steps=[incoming[:]]
for u,v in [(0,1),(1,2),(0,3),(3,4)]:
 incoming=incoming[:];incoming[u]^=1<<v;incoming[v]^=1<<u;steps.append(incoming)
permutations=[]
for p in itertools.permutations(range(7)):
 if p[2]!=0: continue
 if all(all(((r['rows'][i]>>j)&1)==((step[p[i]]>>p[j])&1) for i in range(7) for j in range(7)) for r,step in zip(ours,steps)):
  permutations.append(p)
if not permutations: raise ValueError('No single root-preserving relabeling matches all five path vertices')
result={'incoming_commit':'9036e4dcb9e1016ec42c11802ca193e6a09e6f34','incoming_path':'04-computation/experiments/plus_minus_two_port_routes_20261104.out','incoming_root':0,'session190_root':2,'incoming_arc_flips':[[0,1],[1,2],[0,3],[3,4]],'permutation_convention':'p[i] is the incoming name corresponding to Session190 name i','permutations':permutations,'all_five_tournaments_match':True}
print(json.dumps(result,indent=2)+'\n',end='')
