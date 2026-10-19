from probe_norm import matrix,rref_mod
from math import comb
from pathlib import Path
import numpy as np,json
out=[]
for p in [5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61]:
 a,o=matrix(p);mons=[(i,j,p-i-j) for i in range(p+1) for j in range(p-i+1)];idx={t:i for i,t in enumerate(mons)}
 t=np.zeros(len(mons),dtype=np.int64)
 for k in range(1,p):t[idx[(k,p-k,0)]]=comb(p,k)//p*(-1)**(p-k)
 ra=len(rref_mod(a,p)[1]);rb=len(rref_mod(np.column_stack((a,t)),p)[1]);row={'p':p,'columns':len(o),'rank':ra,'augmented_rank':rb,'lifts_mod_p2':ra==rb};print(row,flush=True);out.append(row)
(Path(__file__).resolve().parent/'lift_census.json').write_text(json.dumps(out,indent=2)+'\n')
