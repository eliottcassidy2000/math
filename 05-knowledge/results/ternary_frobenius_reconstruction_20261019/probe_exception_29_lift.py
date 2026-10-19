"""Finite interior lift obstruction for the explicit p=29, n=42 certificate."""
from pathlib import Path
import json
import numpy as np
from probe_norm import orbits,norm_monomial,rref_mod
from verify_exception_29 import polynomial
p,n=29,42
mons=[(a,b,n-a-b) for a in range(n+1) for b in range(n-a+1)]
idx={t:i for i,t in enumerate(mons)};o=orbits(n,True)
A=np.zeros((len(mons),len(o)),dtype=np.int64)
for j,orb in enumerate(o):
 for t in orb:
  for u,c in norm_monomial(t).items():A[idx[u],j]+=c%p
NE={}
for t,c in polynomial().items():
 for u,v in norm_monomial(t).items():NE[u]=NE.get(u,0)+c*v
if any(v%p for v in NE.values()):raise RuntimeError('integral norm not divisible by p')
b=np.array([-(NE.get(t,0)//p)%p for t in mons],dtype=np.int64)
r=len(rref_mod(A,p)[1]);aug=len(rref_mod(np.column_stack([A,b]),p)[1])
receipt={'prime':p,'degree':n,'input':'cyclic all-interior homogeneous polynomials',
'columns':len(o),'rows':len(mons),'rank':r,'augmented_rank':aug,
'kernel_dimension':len(o)-r,'lift_mod_prime_squared_exists':r==aug}
(Path(__file__).resolve().parent/'exception_29_lift.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
