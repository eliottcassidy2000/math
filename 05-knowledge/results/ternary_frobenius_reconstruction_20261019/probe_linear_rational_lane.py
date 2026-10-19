from rational_kernel import getmatrix
from probe_norm import rref_mod
from pathlib import Path
import json
out=[]
for d in range(31):
 for eps in (1,-1):
  a,o=getmatrix(d,1,eps);null=len(o)-len(rref_mod(a,100003)[1]);row={'original_degree':d,'auxiliary_degree':1,'symmetry':eps,'rows':len(a),'columns':len(o),'nullity_mod_100003':null};print(row,flush=True);out.append(row)
(Path(__file__).resolve().parent/'rational_linear_lane.json').write_text(json.dumps(out,indent=2)+'\n')
