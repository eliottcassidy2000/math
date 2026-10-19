"""Audit the exact two-phase committed prefix, without claiming full coverage."""
from pathlib import Path
from collections import Counter
import json,gzip,hashlib
ROOT=Path(__file__).resolve().parent;OLD=ROOT.parent/'canonical331-stage'
def ck(x,m):
 if not x:raise RuntimeError(m)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
state=json.loads((ROOT/'progress.json').read_text());old=json.loads((OLD/'progress.json').read_text());cover=json.loads(gzip.decompress((OLD/'analytic_cover.json.gz').read_bytes()))
for name,h in state['dependencies'].items():ck(sha(ROOT/name)==h,'continuation source')
for name,h in state['phase_one_sources'].items():ck(sha(OLD/name)==h,'old source')
ck(state['status']==old['status']=='EXACT_CONE_OBSTACLE','stopped state')
seen=[];total=Counter();phases=[]
for folder,expected in [(OLD,old),(ROOT,state)]:
 count=0;h=Counter();ck((folder/'certificates.jsonl.gz').stat().st_size==expected['committed_bytes'],'exact committed bytes')
 with gzip.open(folder/'certificates.jsonl.gz','rt') as f:
  for line in f:
   r=json.loads(line);i=r['core_index'];ck(not seen or i>seen[-1],'disjoint increasing ids');ck(cover['cores'][i]=={'pred':r['pred']},'core identity');seen.append(i);count+=1
   for k in ['residual_orders','analytic_directed_orders','short_span_orders','canonical_excluded_orders','canonical_fully_removed']:h[k]+=r.get(k,0)
   h.update(r.get('hist',{}))
 ck(count==expected['completed'] and dict(h)==expected['hist'],'phase accounting');total.update(h);phases.append({'count':count,'record_sha256':sha(folder/'certificates.jsonl.gz')})
expected=[i for i,item in enumerate(cover['cores']) if len(item)==1];ck(seen==expected[:len(seen)],'complete prefix');ck(expected[len(seen)]==state['obstacle_core'],'first unprocessed obstacle')
out={'status':'EXACT_TWO_PHASE_PREFIX_COVERAGE_AUDIT_PASSED','scope':'Partial scan only; no full morphology theorem','completed':len(seen),'residual_target':len(expected),'next_obstacle_core':state['obstacle_core'],'phases':phases,'hist':dict(total)}
(ROOT/'stopped_prefix_audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
