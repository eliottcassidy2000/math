from fractions import Fraction as F
from collections import defaultdict
import json
# Independent truth-table basis: affine polynomials on three bits.
X=list(range(16)); ev=lambda x,y: ((x&1)+sum(((x>>(i+1))&1)*((y>>i)&1) for i in range(3)))%2
fs=defaultdict(list)
for x in X:fs[ev(x,0),ev(x,1)].append(x)
# Vector statistics explicitly mix visible, invisible, and identity directions.
f=lambda x:(ev(x,0),ev(x,2),int(x==0))
mu=[sum(F(f(x)[i],16) for x in X) for i in range(3)]
full=[[sum((F(f(x)[i])-mu[i])*(F(f(x)[j])-mu[j])/16 for x in X) for j in range(3)] for i in range(3)]
obs=[[F(0) for j in range(3)] for i in range(3)]
loss=[[F(0) for j in range(3)] for i in range(3)]
for xs in fs.values():
 m=[sum(F(f(x)[i],len(xs)) for x in xs) for i in range(3)]
 for i in range(3):
  for j in range(3):
   obs[i][j]+=F(len(xs),16)*(m[i]-mu[i])*(m[j]-mu[j])
   loss[i][j]+=sum((F(f(x)[i])-m[i])*(F(f(x)[j])-m[j])/16 for x in xs)
for i in range(3):
 for j in range(3):
  if full[i][j]-obs[i][j]!=loss[i][j]:raise RuntimeError('matrix identity')
# All 16 Fourier coefficients, independently built from dot-product characters.
for i in range(3):
 for j in range(3):
  z=F(0)
  for ell in range(16):
   chi=lambda x:(-1)**((ell&x).bit_count()%2)
   visible=all(len({chi(x) for x in xs})==1 for xs in fs.values())
   if not visible:z+=sum(F(f(x)[i]*chi(x),16) for x in X)*sum(F(f(x)[j]*chi(x),16) for x in X)
  if z!=loss[i][j]:raise RuntimeError('matrix Fourier')
print(json.dumps({'matrix_identity':'PASS','full':[[str(x) for x in r] for r in full],'observed':[[str(x) for x in r] for r in obs],'lost':[[str(x) for x in r] for r in loss]},indent=2))
