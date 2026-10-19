from probe_norm import orbits,norm_monomial,rref_mod,basis_mod,poly
import numpy as np,json
from pathlib import Path
out=[]
for p in (23,29,31):
 for n in range(3,2*p+13):
  mons=[(a,b,n-a-b) for a in range(n+1) for b in range(n-a+1)];idx={t:i for i,t in enumerate(mons)};o=orbits(n,True)
  a=np.zeros((len(mons),len(o)),dtype=np.int64)
  for j,orb in enumerate(o):
   for t in orb:
    for u,c in norm_monomial(t).items():a[idx[u],j]+=c%p
  rank=len(rref_mod(a,p)[1]);null=len(o)-rank
  print(p,n,null,flush=True)
  if null:
   row={'p':p,'degree':n,'nullity':null,'witness':poly(basis_mod(a,p)[0],o,p)};out.append(row);break
(Path(__file__).resolve().parent/'next_thresholds.json').write_text(json.dumps(out,indent=2)+'\n')
