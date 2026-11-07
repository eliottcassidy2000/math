from probe import make
from itertools import product,combinations
from functools import lru_cache
import json,time
init,step,val=make(2,[[0,1],[0,0]])
@lru_cache(maxsize=200000)
def H(w):
 v=init
 for m in w:v=step(v,m)
 return val(v)
counts=[0]*4;survivors=[];t=time.time()
# Full four-letter central language, all256 global letter maps, insertion≤2.
wrap=[(w[:cut],w[cut:]) for l in range(1,3) for w in product(range(4),repeat=l) for cut in range(l+1)]
for l in range(4,7):
 for seed in product(range(4),repeat=l):
  h=H(seed)
  if h<23:continue
  for cut in range(l+1):
   P,S=seed[:cut],seed[cut:];singles=[H(P+(c,)+S) for c in range(4)]
   sols={-1:[],1:[]}
   for u,v in wrap:
    hnew=H(P+u+v+S);sgn=hnew-2*h
    if sgn not in (-1,1):continue
    counts[0]+=1;out=[H(P+u+(d,)+v+S) for d in range(4)]
    imgs=[[d for d in range(4) if out[d]==2*singles[c]+sgn] for c in range(4)]
    if not all(imgs):continue
    counts[1]+=1
    for phi in product(*imgs):
     ok=True
     for d in range(2,4):
      for w in product(range(4),repeat=d):
       if H(P+u+tuple(phi[c] for c in w)+v+S)!=2*H(P+w+S)+sgn:ok=False;break
      if not ok:break
     if ok:sols[sgn].append((u,v,phi));counts[2]+=1
   if all(sols.values()):survivors.append((seed,cut,sols));counts[3]+=1
print(json.dumps({'seed_max_length':6,'alphabet_size':4,'insertion_max_length':2,'test_depth':3,'counts':counts,'survivors':survivors,'seconds':time.time()-t},indent=2))
