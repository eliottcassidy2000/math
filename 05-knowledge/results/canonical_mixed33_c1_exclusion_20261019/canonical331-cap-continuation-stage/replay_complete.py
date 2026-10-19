#!/usr/bin/env python3
"""Complete independent three-phase replay with immutable source locks."""
from pathlib import Path
from collections import Counter
import gzip,json,hashlib,time,importlib.util
from catalogue import catalogue
import exact_replay_engine as R
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent
PHASES=[BASE/'canonical331-stage',BASE/'canonical331-continuation-stage',ROOT]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 start=time.time();state=json.loads((ROOT/'progress.json').read_text());R.check(state['status']=='ALL_CANONICAL_CORES_EXACT_CERTIFIED_AND_REPLAYED','not complete')
 for name,h in state['dependencies'].items():R.check(sha(ROOT/name)==h,'phase-three sources')
 for name,h in json.loads((ROOT/'source_lock.json').read_text()).items():R.check(sha(BASE/name)==h,'prior source lock')
 cover=json.loads(gzip.decompress((PHASES[0]/'analytic_cover.json.gz').read_bytes()));cat=catalogue();valid={p for p in cat if all(8-p[v].bit_count()-sum(bool(q>>v&1) for q in p)<=7 for v in range(9))}
 R.check(len(cat)==96428 and len(valid)==91604 and valid=={tuple(x['pred']) for x in cover['cores']},'independent catalogue')
 for item in cover['cores']:
  if 'balanced_pair' in item:R.replay(tuple(item['pred']),{'analytic_balanced':item['balanced_pair']})
  elif 'cycle_arcs' in item:R.replay(tuple(item['pred']),{'analytic_cycle':item['cycle_arcs']})
 seen=[];total=Counter();phases=[]
 for k,folder in enumerate(PHASES):
  spec=importlib.util.spec_from_file_location('replay_phase_'+str(k),folder/'exact_replay_engine.py');verifier=importlib.util.module_from_spec(spec);spec.loader.exec_module(verifier)
  expected=json.loads((folder/'progress.json').read_text());h=Counter();count=0
  R.check((folder/'certificates.jsonl.gz').stat().st_size==expected['committed_bytes'],'committed bytes')
  with gzip.open(folder/'certificates.jsonl.gz','rt') as f:
   for line in f:
    rec=json.loads(line);i=rec['core_index'];R.check(not seen or i>seen[-1],'ordered disjoint phases');R.check(cover['cores'][i]=={'pred':rec['pred']},'core identity');verifier.replay(tuple(rec['pred']),rec);seen.append(i);count+=1
    for key in ['residual_orders','analytic_directed_orders','short_span_orders','canonical_excluded_orders','canonical_fully_removed']:h[key]+=rec.get(key,0)
    h.update(rec.get('hist',{}))
    if len(seen)%5000==0:print(json.dumps({'replayed':len(seen),'seconds':time.time()-start}),flush=True)
  R.check(count==expected['completed'] and dict(h)==expected['hist'],'phase totals');total.update(h);phases.append({'records':count,'sha256':sha(folder/'certificates.jsonl.gz')})
 R.check(seen==[i for i,x in enumerate(cover['cores']) if len(x)==1] and len(seen)==80914,'complete union coverage');R.check(phases[-1]['sha256']==state['record_sha256'],'terminal hash')
 out={'status':'FULL_THREE_PHASE_INDEPENDENT_EXACT_REPLAY_PASSED','natural':len(cat),'valid':len(valid),'residual':len(seen),'phases':phases,'hist':dict(total),'seconds':time.time()-start};(ROOT/'complete_replay.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)
if __name__=='__main__':main()
