from probe import make
from itertools import product
import json,time
init,step,val=make(2,[[0,1],[0,0]])
def H(w):
 v=init
 for c in w:v=step(v,c)
 return val(v)
t=time.time();nontrivial=0;pure=[];small=[];best=None
for l in range(0,7):
 for seed in product(range(4),repeat=l):
  for cut in range(l+1):
   for c in range(4):
    f=[H(seed[:cut]+(c,)*m+seed[cut:]) for m in (1,2,3)]
    A=(f[2]-3*f[1]+2*f[0])//6;B=(4*f[1]-f[2]-3*f[0])//2;C=(f[2]-5*f[1]+6*f[0])//2
    assert all(A*3**m+B*2**m+C==f[m-1] for m in (1,2,3))
    if A:
     nontrivial+=1
     if B==0:pure.append((seed,cut,c,A,B,C))
     if best is None or abs(B)<best[0]:best=(abs(B),seed,cut,c,A,B,C)
print(json.dumps({'context_max_length':6,'nonzero_3_mode':nontrivial,'pure_3_mode':pure,'pure_3_count':len(pure),'smallest_abs_2_mode':best,'seconds':time.time()-t},indent=2))
