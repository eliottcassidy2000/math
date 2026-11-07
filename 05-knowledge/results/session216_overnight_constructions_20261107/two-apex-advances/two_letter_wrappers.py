from probe import make
from itertools import product,combinations
from functools import lru_cache
import json,time
init,step,val=make(2,[[0,1],[0,0]])
@lru_cache(maxsize=300000)
def H(w):
 v=init
 for m in w:v=step(v,m)
 return val(v)
counts=[0]*4;survivors=[];t=time.time()
for alpha in combinations(range(4),2):
 wrap=[(w[:cut],w[cut:]) for l in range(1,4) for w in product(alpha,repeat=l) for cut in range(l+1)]
 for l in range(4,7):
  for seed in product(range(4),repeat=l):
   h=H(seed)
   if h<23:continue
   for cut in range(l+1):
    P,S=seed[:cut],seed[cut:];singles=[H(P+(c,)+S) for c in alpha];sols={-1:[],1:[]}
    for u,v in wrap:
     sgn=H(P+u+v+S)-2*h
     if sgn not in (-1,1):continue
     counts[0]+=1;out={d:H(P+u+(d,)+v+S) for d in alpha}
     imgs=[[d for d in alpha if out[d]==2*singles[c]+sgn] for c in range(2)]
     if not all(imgs):continue
     counts[1]+=1
     for image in product(*imgs):
      phi=dict(zip(alpha,image));ok=True
      for d in range(2,5):
       for w in product(alpha,repeat=d):
        if H(P+u+tuple(phi[c] for c in w)+v+S)!=2*H(P+w+S)+sgn:ok=False;break
       if not ok:break
      if ok:sols[sgn].append((u,v,image));counts[2]+=1
    if all(sols.values()):survivors.append((seed,cut,alpha,sols));counts[3]+=1
print(json.dumps({'seed_max_length':6,'alphabet_size':2,'insertion_max_length':3,'test_depth':4,'counts':counts,'survivors':survivors,'seconds':time.time()-t},indent=2))
