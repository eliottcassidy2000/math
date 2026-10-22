import itertools,json
import numpy as np
from decorated import compiled
out=[]
# Five basis: R+,R-,Q++,Q+-,Q--. IdentityQ--=-R+-R--Q++.
change=np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1],[-1,-1,-1,0]],dtype=np.int64)
for w in (1,2):
 for d,c in compiled(w):
  rounded=np.rint(c.Wb).astype(np.int64)
  if not np.array_equal(c.Wb,rounded):raise RuntimeError('noninteger symbolic weight')
  W=rounded
  for _ in d['edges']:W=np.tensordot(W,change,axes=([0],[0]))
  terms=[{'indices':list(idx),'numerator':int(W[idx])} for idx in np.ndindex(W.shape) if W[idx]]
  out.append(dict(d,basis=['R+','R-','Q++','Q+-'],denominator=c.G,terms=terms))
print(json.dumps(out,indent=2))
