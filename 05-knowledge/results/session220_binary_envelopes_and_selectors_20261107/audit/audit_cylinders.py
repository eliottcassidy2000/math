import importlib.util,itertools,json
from pathlib import Path
p=Path(__file__).resolve().parent.parent/'extrema/verify.py'
# Source-workspace fallback; the published packet uses the portable path above.
if not p.is_file():
 p=Path(__file__).resolve().parent.parent/'session220-three-core/verify.py'
s=importlib.util.spec_from_file_location('v',p);v=importlib.util.module_from_spec(s);s.loader.exec_module(v)
checks=0
for c in v.states[::7]:
 for a,b in ((0,4),(1,3),(2,2),(3,1),(4,0)):
  for fixed in itertools.product((None,0,1),repeat=4):
   if a and fixed[a-1]!=1:continue
   pb,rb=fixed[:a],fixed[a:]
   poss=lambda bs:[n for n in range(1<<len(bs)) if all(z is None or ((n>>i)&1)==z for i,z in enumerate(bs))]
   expected=min(v.value(c,p,r,b) for p in poss(pb) for r in poss(rb))
   assert expected==v.mincut(c,pb,rb)
   checks+=1
large_checks=0
for c in v.states[::7]:
 scaled=tuple(z*10**150 for z in c)
 pb,rb=[None,1],[None,None]
 expected=min(v.value(scaled,p,r,2) for p in range(2,4) for r in range(4))
 assert expected==v.mincut(scaled,pb,rb)
 large_checks+=1
print(json.dumps({'status':'PASS','independent_fixed_bit_domains':checks,'large_capacity_checks':large_checks,'capacity_scale':'10^150','tested_core_indices':list(range(0,len(v.states),7))},indent=2))
