"""Independent graph-search replay over inherited positive transfer closure.
Unlike producer's topological bitset DP, uses explicit graph reachability per state.
Transfer construction itself is inherited and has separate literal graph controls.
"""
from pathlib import Path
import sys,contextlib,io,json
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'wildcard'))
from capped import build,minimize
with contextlib.redirect_stdout(io.StringIO()):
 vs,ts,ws,hs=build(127);tr,out,part=minimize(ts,hs)
R=[]
for s in range(len(tr)):
 seen={s};q=[s]
 while q:
  u=q.pop()
  for v in tr[u]:
   if v not in seen:seen.add(v);q.append(v)
 R.append(seen)
classes={h:{i for i,x in enumerate(out) if x==h} for h in set(out)}
def ready(s,W):
 h=out[s]
 return h>=3 and 2*h+1<=127 and bool(R[s]&classes.get(2*h-1,set())&W) and bool(R[s]&classes.get(2*h+1,set())&W)
W0={i for i,h in enumerate(out) if h>0};W1={i for i in W0 if ready(i,W0)};W2={i for i in W0 if ready(i,W1)}
expected=json.loads((Path(__file__).resolve().parents[1]/'wildcard/READINESS.json').read_text())
for row in expected['rows']:
 h=row['h'];ss=classes.get(h,set());masks=set()
 for s in ss:masks.add(sum(1<<b for b,t in enumerate((2*h-1,2*h+1)) if R[s]&classes.get(t,set())))
 assert sorted(masks)==row['available_branch_masks']
 assert bool(ss&W1)==row['common_parent_exists']
 if h<=31:assert bool(ss&W2)==row['two_level_common_parent_exists']
print(json.dumps({'status':'PASS','raw_live_states':len(vs),'count_Moore_states':len(tr),'source_h_range':'odd3..63','common_parent_h':sorted({out[i] for i in W1 if out[i]<=63}),'two_level_h_through31':sorted({out[i] for i in W2 if out[i]<=31}),'method':'explicit DFS reachability independently replays producer topological-bitset readiness on inherited transfer closure'},indent=2))
