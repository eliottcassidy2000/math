from functools import lru_cache
from time import monotonic

def paley(p):
 sq={i*i%p for i in range(1,p)}
 return [[0 if i==j else 1 if (j-i)%p in sq else -1 for j in range(p)]for i in range(p)]
def solve(K):
 n=len(K);N=1<<n;sizes=[i.bit_count() for i in range(N)]
 rows=[]
 for row in K:
  r=[0]*N
  for s in range(1,N):
   bit=s&-s;r[s]=r[s^bit]+row[bit.bit_length()-1]
  rows.append(r)
 @lru_cache(None)
 def cross(a,b):
  total=0
  while a:
   bit=a&-a;a^=bit;total+=rows[bit.bit_length()-1][b]
  return total
 f=[0]*N;tree=[None]*N;parts=[0,0]
 for s in range(1,N):
  if sizes[s]<2:continue
  anchor=s&-s;rest=s^anchor;a0=rest
  best=sizes[s]*(sizes[s]-1)//2
  while True:
   a=a0|anchor;bc=s^a
   if bc:
    z=cross(a,bc);val=f[a]+f[bc]+(sizes[a]*sizes[bc]-abs(z))//2;parts[0]+=1
    if val<best:best=val;tree[s]=('join',a,bc,1 if z>=0 else -1)
    anchor2=bc&-bc;rest2=bc^anchor2;b0=rest2
    while True:
     b=b0|anchor2;c=bc^b
     if c:
      z=cross(a,b)+cross(b,c)+cross(c,a);val=f[a]+f[b]+f[c]+(sizes[a]*sizes[b]+sizes[b]*sizes[c]+sizes[c]*sizes[a]-abs(z))//2;parts[1]+=1
      if val<best:best=val;tree[s]=('cycle',a,b,c,1 if z>=0 else -1)
     if b0==0:break
     b0=(b0-1)&rest2
   if a0==0:break
   a0=(a0-1)&rest
  f[s]=best
 return f,tree,parts
def expand(s, tree):
 if s.bit_count()==1:return {"vertex":s.bit_length()-1}
 item=tree[s]
 return {"kind":item[0],"orientation":item[-1],"children":[expand(t,tree) for t in item[1:-1]]}
if __name__=="__main__":
 import json,pathlib
 result={"normalization":"K_uv=+1 iff u->v; individual unordered arc flips","grammar":"arbitrary-depth disjoint ordered-join/C3 substitution; named singleton leaves","cases":[]}
 for n in (3,7,11):
  f,tr,parts=solve(paley(n))
  result["cases"].append({"n":n,"minimum":f[-1],"subset_lower_bounds":f,"root_partition_counts":parts,"witness":expand((1<<n)-1,tr)})
 print(json.dumps(result,sort_keys=True,indent=2))
