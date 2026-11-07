import functools,itertools,json
from pathlib import Path
# Load only independently implemented compiler definitions, not their test driver.
p=Path(__file__).parent
ns={'__file__':str(p/'independent.py')};exec((p/'independent.py').read_text().split('rng=random.Random(216);controls=0')[0],ns)
init,step,val,_=ns['compiler'](2,[[0,1],[0,0]])
@functools.lru_cache(maxsize=150000)
def H(w):
 v=init
 for m in w:v=step(v,m)
 return val(v)
macros=[w for n in range(1,4) for w in itertools.product(range(4),repeat=n)]
@functools.lru_cache(maxsize=100000)
def kids(w,kind):
 out=[set(),set()];h=H(w)
 if kind=='boundary':
  candidates=(q for u in macros for q in (u+w,w+u))
 else:
  candidates=[]
  for end in (0,len(w)-1):
   for letter in range(4):
    q=list(w);q[end]=letter;q=tuple(q)
    for m in range(4):candidates.extend(((m,)+q,q+(m,)))
 for q in candidates:
  count=H(q)
  if count==2*h-1:out[0].add(q)
  if count==2*h+1:out[1].add(q)
 return tuple(tuple(sorted(s)) for s in out)
@functools.lru_cache(maxsize=100000)
def viable(w,d,kind):
 return d==0 or all(any(viable(q,d-1,kind) for q in branch) for branch in kids(w,kind))
result={}
for kind,maxdepth in [('boundary',2),('recolor',3)]:
 counts=[0]*(maxdepth+1);depth2=[]
 for n in range(4,7):
  for w in itertools.product(range(4),repeat=n):
   if H(w)<23:continue
   for d in range(1,maxdepth+1):
    if not viable(w,d,kind):break
    counts[d]+=1
    if d==2:depth2.append(w)
 result[kind]={'counts':counts,'depth2_seeds':depth2}
if result['boundary']['counts']!=[0,8,0] or result['recolor']['counts']!=[0,267,2,0]:raise RuntimeError(result)
print(json.dumps({'method':'Independent compressed transfer; raw words, no dual quotient','results':result},indent=2))
# Literal subset-DP verification of one depth-two positive control tree.
def tree(w,d):
 if ns['mod'].hp(w)!=H(w):raise RuntimeError('literal tree count')
 out={'word':w,'H':H(w)}
 if d:out['children']=[tree(next(q for q in branch if viable(q,d-1,'recolor')),d-1) for branch in kids(w,'recolor')]
 return out
(p/'ADAPTIVE_TREE.json').write_text(json.dumps(tree((1,1,2,2),2),indent=2)+'\n')
