"""Independent full enumeration of two finite binary-bank kernels, not E(Q)/2E(Q)."""
import collections,copy,json
from pathlib import Path
D=Path(__file__).resolve().parent
checks=0
def need(b):
 global checks
 checks+=1
 if not b:raise ValueError('kernel rejected')
def rank(vectors):
 a=list(vectors);r=0
 for j in range(30):
  pivot=next((k for k in range(r,len(a))if a[k]>>j&1),None)
  if pivot is None:continue
  a[r],a[pivot]=a[pivot],a[r]
  for k in range(len(a)):
   if k!=r and a[k]>>j&1:a[k]^=a[r]
  r+=1
 return r
def verify(c,k):
 rows=[r['bits']for r in c['character_rows']];basis=k['kernel_basis_masks'];need(k['curve_id']==c['curve']['id'])
 need(rank(rows)==k['bank_rank']==17);need(rank(basis)==k['kernel_dimension']==13)
 need(all(0<v<2**30 for v in basis))
 for v in basis:
  for r in rows:need((v&r).bit_count()%2==0)
 # Independent mask-to-sum enumeration rather than iterative list doubling.
 vectors=[]
 for mask in range(2**13):
  v=0
  for j in range(13):
   if mask>>j&1:v^=basis[j]
  vectors.append(v)
 need(len(set(vectors))==8192)
 counts=collections.Counter(v.bit_count()for v in vectors)
 need(dict(counts)=={int(w):n for w,n in k['weight_enumerator'].items()})
 minimum=min(w for w in counts if w);need(minimum==k['minimum_nonzero_weight'])
 actual={tuple(i for i in range(30)if v>>i&1)for v in vectors if v.bit_count()==minimum}
 need(actual==set(map(tuple,k['minimum_supports_zero_based'])))
 return {'curve':k['curve_id'],'dimension':13,'minimum_weight':minimum,'minimum_support_count':len(actual)}
cs=json.loads((D/'halving_certificate.json').read_text())['curves'];ks=json.loads((D/'binary_kernel_certificate.json').read_text())
results=[verify(c,k)for c,k in zip(cs,ks)]
rejected=0
for what in ('basis','weight','support'):
 k=copy.deepcopy(ks[0])
 if what=='basis':k['kernel_basis_masks'][0]^=1
 if what=='weight':k['minimum_nonzero_weight']-=1
 if what=='support':k['minimum_supports_zero_based'][0][0]=0
 try:verify(cs[0],k)
 except ValueError:rejected+=1
 else:raise ValueError('hostile passed')
need(rejected==3)
print(json.dumps({'results':results,'checks':checks,'hostiles_rejected':rejected,'warning':'Finite-bank kernels are NOT certified rational halving kernels'},sort_keys=True))
