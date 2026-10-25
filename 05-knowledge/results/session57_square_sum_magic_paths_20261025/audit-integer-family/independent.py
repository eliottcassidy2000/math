import sys, math, json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent/'integer'))
from primary_construction import family

def require(x,m):
 if not x: raise ValueError(m)
# Exhaustive affine-vector check, independent of the supplied projective table.
coeffs=[[16608,102220,205548,135135],[5760,50208,168500,252948,141345],[19680,113996,216300,135135]]
checks=0
for u in range(47):
 for v in range(47):
  if u==v==0:continue
  for C in coeffs:
   n=len(C)-1
   require(sum(c*pow(u,n-i,47)*pow(v,i,47) for i,c in enumerate(C))%47!=0,'finite field zero')
   checks+=1
best=None; count=valid=0
for u in range(-200,201):
 for v in range(-200,201):
  if math.gcd(u,v)!=1:continue
  count+=1; A,B=family(u,v)
  if len(set(A))<4 or len(set(B))<4:continue
  sums=[a+b for a in A for b in B]; roots=[math.isqrt(s) for s in sums]
  require(all(r*r==s for r,s in zip(roots,sums)),'sum not square')
  h=math.gcd(*roots); b0=min(B)
  require(all((a+b0)%(h*h)==0 for a in A) and all((b-b0)%(h*h)==0 for b in B),'nonintegral normalization')
  A=sorted((a+b0)//(h*h) for a in A);B=sorted((b-b0)//(h*h) for b in B)
  lo,hi=1,min(A)-1
  if hi<lo:continue
  forbidden={(a-b)//2 for a in A for b in B if (a-b)%2==0}
  k=(max(A)-max(B))//2
  # Convex objective: first permitted point to each side of real minimizer.
  left=min(hi,k);right=max(lo,k+1)
  while left>=lo and left in forbidden:left-=1
  while right<=hi and right in forbidden:right+=1
  candidates=[t for t in (left,right) if lo<=t<=hi]
  if not candidates:continue
  valid+=1
  N,t=min((max(max(A)-t,max(B)+t),t) for t in candidates)
  item=(N,u,v,t)
  if best is None or item<best:best=item
require((count,valid,best)==(97856,97780,(1112833,-11,5,1)),'box optimum mismatch')
print(json.dumps({'status':'PASS','nonzero_F47_vectors_times_forms':checks,'primitive_pairs':count,'valid_pairs':valid,'optimum_N_u_v_t':best},indent=2))
