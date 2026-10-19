"""Independent sparse-column verification of the mod-p^2 lift obstruction."""
from verify_periods_independent import independent_norm,sparse_rank,require
from math import comb
from pathlib import Path
import json
out=[]
for p in [5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61]:
 cols=[]
 for a in range(p+1):
  for b in range(p-a+1):
   c=p-a-b;t=(a,b,c);orb={t,(b,c,a),(c,a,b)}
   if t!=min(orb):continue
   col={}
   for u in orb:
    for v,d in independent_norm(u,p).items():col[v]=(col.get(v,0)+d)%p
   cols.append({v:d for v,d in col.items() if d})
 target={}
 for k in range(1,p):
  require(comb(p,k)%p==0,'non-integral divided norm')
  target[(k,p-k,0)]=(comb(p,k)//p*(-1)**(p-k))%p
 _,r=sparse_rank((v.copy() for v in cols),p)
 _,rr=sparse_rank((v.copy() for v in cols+[target]),p)
 require(rr==r+1,f'lift obstruction failed at p={p}')
 row={'p':p,'columns':len(cols),'rank':r,'augmented_rank':rr,'lifts_mod_p2':False};out.append(row);print(row,flush=True)
(Path(__file__).resolve().parent/'lift_independent.json').write_text(json.dumps(out,indent=2)+'\n')
print('All 16 lift obstructions independently verified.',flush=True)
