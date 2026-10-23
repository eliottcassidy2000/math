from math import isqrt
from pathlib import Path
import json
D=lambda n:isqrt(2*n-1)-isqrt(n)
rows=[]
for j in range(101):
 # floor((1+sqrt2)*(j+1)) = (j+1)+floor(sqrt(2*(j+1)^2)).
 rmax=(j+1)+isqrt(2*(j+1)**2)
 ranges=[]; values=set()
 for r in range(1,rmax+1):
  L=max(r*r,((r+j)**2+2)//2);U=min((r+1)**2-1,(r+j+1)**2//2)
  if L<=U:ranges.append([L,U]);values.update(range(L,U+1))
 direct={n for n in range(1,(rmax+1)**2) if D(n)==j}
 if values!=direct:raise RuntimeError(j)
 rows.append({'degree':j,'count':len(values),'intervals':ranges})
if [r['count'] for r in rows[:10]]!=[3,12,24,34,47,59,68,84,92,105]:raise RuntimeError('prefix')
out={'status':'PASS exact integer controls','j_max':100,'proof_scope':'all integer j>=0 by floor-interval equivalence','rows':rows}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print('PASS degrees0..100, complete finite windows')
