from common import *
from kernel import *
from numba import njit
import numpy as np,gzip,time
@njit
def taxmat(rows):
 n=len(rows);size=1<<n;full=size-1
 a=np.zeros((size,n),dtype=np.int64);b=np.zeros((size,n),dtype=np.int64)
 for s in range(1,size):
  for v in range(n):
   if s>>v&1:
    t=s^(1<<v)
    if not t:a[s,v]=b[s,v]=1
    else:
     for u in range(n):
      if rows[u]>>v&1:a[s,v]+=a[t,u]
      if rows[v]>>u&1:b[s,v]+=b[t,u]
 q=np.zeros((n,n),dtype=np.int64)
 for s in range(1,full):
  for u in range(n):
   for v in range(n):q[u,v]+=a[s,u]*b[full^s,v]
 return np.sum(a[full]),q+q.T

def delta(tab,O):
 h,w=tab;return int(h)-sum(int(w[u,v]) for u,v in itertools.combinations(O,2))
if __name__=='__main__':
 bank=load_bank()
 ans=[];total=0;t=time.time()
 for r in bank:
  if r['n']>8 or not r['strong'] or not r['seeds']:continue
  rows=r['rows'];n=len(rows);sels=sorted({r['attachment_orbit'][q[0]] for q in r['seeds']})
  tabs={}
  for p in itertools.combinations(range(n),2):
   masks=[m for m in sels if all(not m>>v&1 for v in p)]
   if not masks:continue
   child,_=expand(rows,[],p,[2,2]);tb=taxmat(np.array(child,dtype=np.int64))
   for m in masks:
    O=[v for v in range(n) if m>>v&1];_,Oc=expand(rows,O,p,[2,2]);total+=1
    if delta(tb,Oc)!=1:continue
    k=exp(skeleton(rows,O,p));diag=defaultdict(int)
    for (a,b),c in k.items():diag[a*b]+=c
    diag={a:c for a,c in diag.items() if c}
    # true compensated diagonal, not jointly constant
    if diag=={1:1} and k!={(1,1):1}:
     q={'rows':rows,'O':O,'ports':p,'kernel':[[a,b,c] for (a,b),c in sorted(k.items())]};ans.append(q);print('FOUND',q,flush=True)
  if r['rep']%200==0:print('progress',n,r['rep'],total,time.time()-t,flush=True)
 (HERE/'screen.json').write_text(json.dumps({'tested':total,'compensated_diagonals':ans},indent=2));print('DONE',total,len(ans))
