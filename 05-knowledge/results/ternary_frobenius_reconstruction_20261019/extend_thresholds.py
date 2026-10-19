from probe_norm import matrix,rref_mod,basis_mod,poly
from pathlib import Path
import json
first={}; rows=[]
for n in range(3,61):
 a,o=matrix(n,True); row={'degree':n}
 for p in (13,17,19):
  if p in first: continue
  r=len(rref_mod(a,p)[1]);dim=len(o)-r; row[p]=dim
  if dim:
   v=basis_mod(a,p)[0];first[p]={'degree':n,'kernel_dim':dim,'witness':poly(v,o,p)}
 print(row,flush=True);rows.append(row)
 if len(first)==3:break
(Path(__file__).resolve().parent / 'extended_thresholds.json').write_text(json.dumps({'rows':rows,'first_interior':first},indent=2)+'\n')
