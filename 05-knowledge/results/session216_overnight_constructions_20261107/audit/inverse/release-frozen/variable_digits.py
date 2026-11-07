from digits_boundary import H,dual
from itertools import product
from functools import lru_cache
import json,time
macros=[w for l in range(1,4) for w in product(range(4),repeat=l)]
@lru_cache(maxsize=50000)
def viable(w,depth):
 if depth==0:return True
 h=H(w)
 children={-1:set(),1:set()}
 for u in macros:
  for q in (w+u,u+w):
   hq=H(q)
   for s in [-1,1]:
    if hq==2*h+s:children[s].add(min(q,dual(q)))
 return all(any(viable(q,depth-1) for q in children[s]) for s in [-1,1])
counts=[0]*5;examples={};t=time.time()
for l in range(4,7):
 for w in product(range(4),repeat=l):
  if H(w)<23:continue
  for d in range(1,5):
   if viable(w,d):counts[d]+=1;examples.setdefault(d,w)
   else:break
print(json.dumps({'macro_max_length':3,'seed_max_length':6,'counts':counts,'examples':examples,'seconds':time.time()-t},indent=2))
