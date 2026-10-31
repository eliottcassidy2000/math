from itertools import combinations
from functools import lru_cache
import json
@lru_cache(None)
def h(out):
 n=len(out);d=[[0]*n for _ in range(1<<n)]
 for i in range(n):d[1<<i][i]=1
 for m in range(1,1<<n):
  for i in range(n):
   if d[m][i]:
    todo=out[i]&~m
    while todo:
     b=todo&-todo;j=b.bit_length()-1;d[m|b][j]+=d[m][i];todo-=b
 return sum(d[-1])
def extend(out,cut):
 n=len(out)
 return tuple(out[i]|(0 if cut>>i&1 else 1<<n) for i in range(n))+(cut,)
def children(out,target):
 for cut in range(1<<len(out)):
  q=extend(out,cut)
  if h(q)==target:yield cut,q

def viable(out,d):
 if not d:return []
 H=h(out);proof=[]
 for delta in [-1,1]:
  for a,q in children(out,2*H+delta):
   tail=viable(q,d-1)
   if tail is not None:proof.append([a,tail]);break
  else:return None
 return proof
if __name__=='__main__':
 for n in [4,5]:
  edges=list(combinations(range(n),2)); tested=0
  for mask in range(1<<len(edges)):
   out=[0]*n
   for j,(a,b) in enumerate(edges):
    if mask>>j&1:out[a]|=1<<b
    else:out[b]|=1<<a
   out=tuple(out);H=h(out)
   if H<9:continue
   tested+=1
   ans=viable(out,2)
   if ans is not None:
    print(json.dumps({'n':n,'out':out,'H':H,'depth2':ans}),flush=True)
    raise SystemExit
  print('n',n,'tested',tested,'cache',h.cache_info(),flush=True)
