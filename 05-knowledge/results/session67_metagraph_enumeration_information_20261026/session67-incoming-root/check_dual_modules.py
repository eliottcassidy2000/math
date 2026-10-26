from itertools import combinations
from collections import Counter
import json

def a(dual=False):
 rows=[14,4,8,2]
 return [((15^(1<<i))^r) if dual else r for i,r in enumerate(rows)]
def make(duals):
 blocks=[a(d) for d in duals]; rows=[0]*12
 for i in range(12):
  b,u=divmod(i,4);rows[i]=(blocks[b][u]<<(4*b)) | (15<<(4*((b+1)%3)))
 return rows

def aut(rows):
 n=len(rows);degs=[r.bit_count() for r in rows];p=[-1]*n;cnt=0
 # Complete constraint search: only exact necessary degree/arc constraints prune.
 def go(used):
  nonlocal cnt
  if used.bit_count()==n:cnt+=1;return
  best=None
  for u in range(n):
   if p[u]>=0:continue
   cand=[v for v in range(n) if not(used>>v&1) and degs[u]==degs[v] and all(p[w]<0 or ((rows[u]>>w&1)==(rows[v]>>p[w]&1)) for w in range(n))]
   if best is None or len(cand)<len(best[1]):best=(u,cand)
  u,cand=best
  for v in cand:p[u]=v;go(used|(1<<v));p[u]=-1
 go(0);return cnt

def cyc(rows):
 return [t for t in combinations(range(len(rows)),3) if all(sum(rows[u]>>v&1 for v in t if v!=u)==1 for u in t)]
T,U=make([0,0,0]),make([1,0,0])
if cyc(T)!=cyc(U):raise RuntimeError('cyclic field changed')
autT,autU=aut(T),aut(U)
if (autT,autU)!=(81,27):raise RuntimeError('aut mismatch')
print(json.dumps({'status':'PASS','T':T,'U':U,'cyclic_triples':len(cyc(T)),'identical_labelled_C3':True,'automorphism_orders':[autT,autU],'degree_histograms':[dict(sorted(Counter(r.bit_count() for r in X).items())) for X in [T,U]]},sort_keys=True))
