from probe import make
from itertools import product
from functools import lru_cache
import json,time
init,step,val=make(2,[[0,1],[0,0]])
@lru_cache(maxsize=100000)
def H(w):
 v=init
 for m in w:v=step(v,m)
 return val(v)
def dual(w):return tuple((3,1,2,0)[m] for m in w[::-1])
actions=list(product(range(2),range(2),range(4),range(4)))
def act(w,a):
 side,end,flip,m=a;q=list(w);q[-1 if end else 0]^=flip;q=tuple(q)
 return (m,)+q if side else q+(m,)
@lru_cache(maxsize=100000)
def children(w,s):
 target=2*H(w)+s
 return tuple(set(min(q,dual(q)) for a in actions if H(q:=act(w,a))==target))
@lru_cache(maxsize=100000)
def viable(w,d):
 if d==0:return True
 return all(any(viable(q,d-1) for q in children(w,s)) for s in [-1,1])
t=time.time();counts=[0]*6;examples={}
for l in range(4,7):
 for w in product(range(4),repeat=l):
  if H(w)<23:continue
  for d in range(1,6):
   if viable(w,d):counts[d]+=1;examples.setdefault(d,w)
   else:break
print(json.dumps({'counts':counts,'examples':examples,'seconds':time.time()-t},indent=2))
