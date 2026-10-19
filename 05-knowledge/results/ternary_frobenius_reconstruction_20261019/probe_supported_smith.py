from probe_norm import matrix
from sympy import Matrix,ZZ
from sympy.matrices.normalforms import smith_normal_form
from pathlib import Path
import json
out=[]
for n in range(1,32,2):
 M,o=matrix(n);js=[j for j,orb in enumerate(o) if min(orb[0])==0]
 D=smith_normal_form(Matrix(M[:,js].tolist()),domain=ZZ);diag=[abs(int(D[i,i])) for i in range(len(js))]
 if any(x==0 for x in diag):raise RuntimeError('unexpected rational supported odd kernel')
 for p in [5,7,11,13,17,19,23,29,31]:
  count=0
  for x in diag:
   if x%p==0:count+=1
   if x%(p*p)==0:raise RuntimeError(('higher p torsion',n,p,x))
  k=0;pw=p
  while pw<=n:k+=1;pw*=p
  if count!=k:raise RuntimeError(('wrong p torsion rank',n,p,count,k))
 row={'degree':n,'smith_nonzero_diagonal':diag};print(row,flush=True);out.append(row)
(Path(__file__).resolve().parent/'supported_smith.json').write_text(json.dumps(out,indent=2)+'\n')
