from probe import make
from itertools import product
import json,time
init,step,val=make(2,[[0,1],[0,0]])
macros=[w for l in range(1,4) for w in product(range(4),repeat=l)]
def apply(v,w):
 for m in w:v=step(v,m)
 return v
seeds=[w for l in range(4,7) for w in product(range(4),repeat=l)]
counts=[0]*5;survivors=[];t=time.time()
for w in seeds:
 v=apply(init,w);h=val(v)
 if h<23:continue
 a=[u for u in macros if val(apply(v,u))==2*h-1]
 b=[u for u in macros if val(apply(v,u))==2*h+1]
 for u in a:
  for z in b:
   counts[0]+=1;front=[v];ok=True
   for d in range(4):
    nxt=[]
    for x in front:
     for macro,sgn in [(u,-1),(z,1)]:
      y=apply(x,macro)
      if val(y)!=2*val(x)+sgn:ok=False;break
      nxt.append(y)
     if not ok:break
    if not ok:break
    counts[d+1]+=1;front=nxt
   if ok:survivors.append((w,u,z))
print(json.dumps({'seed_max_length':6,'macro_max_length':3,'counts':counts,'survivors':survivors,'seconds':time.time()-t},indent=2))
