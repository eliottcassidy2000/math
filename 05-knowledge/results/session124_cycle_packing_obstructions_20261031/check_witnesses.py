"""Independently count literal vertex permutations for the nine relative holes."""
from itertools import permutations,combinations
from pathlib import Path
import json
from check import hp,oddcycles
P=Path(__file__).resolve().parent
checks=0
def need(x,msg):
 global checks
 checks+=1
 if not x:raise RuntimeError(msg)
rows=[]
for key,out in json.load(open(P/'global_witnesses.json')).items():
 h=int(key);n=len(out)
 for i in range(n):
  need(not(out[i]>>i&1),'no loops')
  for j in range(i):need((out[i]>>j&1)+(out[j]>>i&1)==1,'tournament')
 literal=sum(all(out[p[i]]>>p[i+1]&1 for i in range(n-1)) for p in permutations(range(n)))
 need(literal==h==hp(out),'literal permutation count')
 cycles,supports=oddcycles(out);pairs=sum(not(a&b) for a,b in combinations(supports,2))
 need(pairs>0,'requires positive pair coefficient')
 need(1+2*sum(cycles.values())+4*pairs==h,'OCF at order8')
 rows.append(dict(H=h,n=n,literal_paths=literal,cycle_lengths=dict(cycles),disjoint_cycle_pairs=pairs))
print(json.dumps(dict(status='PASS',checks=checks,witnesses=rows),indent=2,sort_keys=True))
