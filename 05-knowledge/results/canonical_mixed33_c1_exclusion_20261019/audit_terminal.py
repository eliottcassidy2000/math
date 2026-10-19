#!/usr/bin/env python3
"""Fast complete three-phase coverage and immutable-source audit."""
from pathlib import Path
from collections import Counter
import json,gzip,hashlib
ROOT=Path(__file__).resolve().parent
PHASES=[ROOT/'canonical331-stage',ROOT/'canonical331-continuation-stage',ROOT/'canonical331-cap-continuation-stage']
def ck(x,m):
 if not x:raise RuntimeError(m)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def leaves(c):
 if c.get('leaf') or 'tree' not in c and 'split_vertex' not in c:return [c]
 t=c.get('tree',c)
 if t.get('leaf'):return[t]
 return leaves(t['A_below'])+leaves(t['B_above'])
cover=json.loads(gzip.decompress((PHASES[0]/'analytic_cover.json.gz').read_bytes()));seen=[];h=Counter();phases=[];certs=0;maxcomparison=maxcap=0
for j,folder in enumerate(PHASES):
 state=json.loads((folder/'progress.json').read_text());ck(state['status']==('ALL_CANONICAL_CORES_EXACT_CERTIFIED_AND_REPLAYED' if j==2 else 'EXACT_CONE_OBSTACLE'),'phase terminal')
 for name,digest in state['dependencies'].items():ck(sha(folder/name)==digest,'immutable phase source '+name)
 if (folder/'source_lock.json').exists():
  for name,digest in json.loads((folder/'source_lock.json').read_text()).items():ck(sha(ROOT/name)==digest,'prior source lock')
 ck((folder/'certificates.jsonl.gz').stat().st_size==state['committed_bytes'],'committed stream size');local=Counter();count=0
 with gzip.open(folder/'certificates.jsonl.gz','rt') as f:
  for line in f:
   r=json.loads(line);i=r['core_index'];ck(not seen or i>seen[-1],'disjoint order');ck(cover['cores'][i]=={'pred':r['pred']},'core identity');seen.append(i);count+=1
   for key in ['residual_orders','analytic_directed_orders','short_span_orders','canonical_excluded_orders','canonical_fully_removed']:local[key]+=r.get(key,0)
   local.update(r.get('hist',{}));certs+=len(r.get('certificates',[]))
   for cert in r.get('certificates',[]):
    for leaf in leaves(cert):
     for a,b,w in leaf['terms']:ck(type(w) is int and w>0,'comparison integer');maxcomparison=max(maxcomparison,w)
     for cap in leaf.get('caps',[]):ck(type(cap[-1]) is int and cap[-1]>0,'cap integer');maxcap=max(maxcap,cap[-1])
 ck(count==state['completed'] and dict(local)==state['hist'],'phase accounting');h.update(local);phases.append({'phase':j+1,'cores':count,'hist':dict(local),'record_sha256':sha(folder/'certificates.jsonl.gz')})
ck(seen==[i for i,item in enumerate(cover['cores']) if len(item)==1] and len(seen)==80914,'complete coverage')
out={'status':'COMPLETE_THREE_PHASE_COVERAGE_AUDIT_PASSED','scope':'canonical(maximum span,maximum common count), mixed(3,3,1)','natural_cores':96428,'degree_valid':91604,'analytic_balanced':10690,'residual_cores':len(seen),'hist':dict(h),'phases':phases,'distinct_certificates':certs,'maximum_comparison_weight':maxcomparison,'maximum_cap_multiplier':maxcap};(ROOT/'audit_receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
