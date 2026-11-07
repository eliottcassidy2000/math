from probe import make
from itertools import product
from functools import lru_cache
import time,json
init,step,val=make(2,[[0,1],[0,0]])
@lru_cache(maxsize=100000)
def H(w):
 v=init
 for m in w:v=step(v,m)
 return val(v)
def dual(w):return tuple((3,1,2,0)[m] for m in w[::-1])
macros=[w for l in range(1,3) for w in product(range(4),repeat=l)]
actions=[(u,side,d) for u in macros for side in range(2) for d in range(2)]
def act(w,a):
 u,side,d=a;r=u+w if side else w+u
 return dual(r) if d else r
counts=[0]*5;survivors=[];t=time.time()
for l in range(4,6):
 for w in product(range(4),repeat=l):
  h=H(w)
  if h<23:continue
  a=[x for x in actions if H(act(w,x))==2*h-1]
  b=[x for x in actions if H(act(w,x))==2*h+1]
  for x in a:
   for y in b:
    counts[0]+=1;front=[w];ok=True
    for depth in range(4):
     nxt=[]
     for q in front:
      for action,s in [(x,-1),(y,1)]:
       r=act(q,action)
       if H(r)!=2*H(q)+s:ok=False;break
       nxt.append(r)
      if not ok:break
     if not ok:break
     counts[depth+1]+=1;front=nxt
    if ok:survivors.append((w,x,y))
print(json.dumps({'seed_max_length':5,'macro_max_length':2,'counts':counts,'survivors':survivors,'seconds':time.time()-t},indent=2))
