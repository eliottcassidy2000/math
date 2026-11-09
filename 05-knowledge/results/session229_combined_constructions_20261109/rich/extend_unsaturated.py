# Explicit sibling loading supports isolated Python (-I) without cwd imports.
import importlib.util
from pathlib import Path
_spec = importlib.util.spec_from_file_location("rich_module_subset_dp", Path(__file__).resolve().with_name("search.py"))
_search = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_search)
data = _search.data
import json, math
from collections import defaultdict
cores=[[64,121,35,37,13,81,28],[20,85,72,195,204,159,33,71]]
fib=defaultdict(dict);count=0
for ci,base in enumerate(cores):
 for orient in range(1<<len(base)):
  a=[w|((1<<len(base)) if not orient>>i&1 else 0) for i,w in enumerate(base)]+[orient]
  H,Q,aw,en,st,h=data(a)
  if math.gcd(H,Q)!=1:continue
  count+=1;n=len(a);N=1<<n
  gw=[sum(en[s][v]*h[(N-1)^s] for s in range(1,N)) for v in range(n)]
  K=[[sum(en[s][i]*st[(N-1)^s][j] for s in range(1,N-1)) for j in range(n)] for i in range(n)]
  for m in range(N):
   vs=[i for i in range(n) if m>>i&1];sig=(H,Q)+tuple(sum(w[i] for i in vs) for w in [aw,gw,st[-1],en[-1]])
   T=sum(K[i][j] for i in vs for j in vs);fib[sig].setdefault(T,dict(core=ci,orientation=orient,adj=a,mask=m,T=T))
best=0;wins=[]
for sig,ts in fib.items():
 for t in ts:
  if t-2 in ts:continue
  run=[];u=t
  while u in ts:run.append(ts[u]);u+=2
  if len(run)>best:best=len(run);wins=[]
  if len(run)==best:wins.append(dict(signature=sig,run=run))
print(json.dumps(dict(status='FINITE-EXACT',universe='all 384 one-vertex extensions of the specified B95 (128) and B389 (256) cores, all cut masks, gcd(H,Q)=1 filter',eligible_extensions=count,interfaces=len(fib),max_run=best,maximal_run_count=len(wins),max_fiber_size=max(map(len,fib.values())),maximal_run_seed_orientations=sorted({(w["core"],w["orientation"]) for z in wins for w in z["run"]}),witnesses=wins),indent=2))
