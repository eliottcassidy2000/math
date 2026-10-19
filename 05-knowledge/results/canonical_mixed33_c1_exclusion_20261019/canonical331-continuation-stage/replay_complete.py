#!/usr/bin/env python3
"""Independent two-phase exact union replay; no numerical dependencies."""
from pathlib import Path
from collections import Counter
import gzip,json,hashlib,time,importlib.util
from catalogue import catalogue
import exact_replay_engine as V2
ROOT=Path(__file__).resolve().parent;OLD=ROOT.parent/'canonical331-stage'
spec=importlib.util.spec_from_file_location('phase_one_replay',OLD/'exact_replay_engine.py');V1=importlib.util.module_from_spec(spec);spec.loader.exec_module(V1)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 t=time.time();state=json.loads((ROOT/'progress.json').read_text());old=json.loads((OLD/'progress.json').read_text())
 V2.check(state['status']=='ALL_CANONICAL_CORES_EXACT_CERTIFIED_AND_REPLAYED','continuation incomplete')
 for name,h in state['dependencies'].items():V2.check(sha(ROOT/name)==h,'new source changed')
 for name,h in state['phase_one_sources'].items():V2.check(sha(OLD/name)==h,'old source changed')
 V2.check(sha(ROOT/'certificates.jsonl.gz')==state['record_sha256'],'new record hash')
 cover=json.loads(gzip.decompress((OLD/'analytic_cover.json.gz').read_bytes()));cat=catalogue();valid={p for p in cat if all(8-p[v].bit_count()-sum(bool(q>>v&1) for q in p)<=7 for v in range(9))}
 V2.check(len(cat)==96428 and len(valid)==91604 and valid=={tuple(x['pred']) for x in cover['cores']},'independent catalogue')
 for item in cover['cores']:
  if 'balanced_pair' in item:V1.replay(tuple(item['pred']),{'analytic_balanced':item['balanced_pair']})
  elif 'cycle_arcs' in item:V1.replay(tuple(item['pred']),{'analytic_cycle':item['cycle_arcs']})
 seen=[];total=Counter();phases=[]
 for folder,replayer,expected in [(OLD,V1,old),(ROOT,V2,state)]:
  h=Counter();count=0
  with gzip.open(folder/'certificates.jsonl.gz','rt') as f:
   for line in f:
    r=json.loads(line);i=r['core_index'];V2.check(not seen or i>seen[-1],'phase order/disjointness');V2.check(cover['cores'][i]=={'pred':r['pred']},'core identity');replayer.replay(tuple(r['pred']),r);seen.append(i);count+=1
    for k in ['residual_orders','analytic_directed_orders','short_span_orders','canonical_excluded_orders','canonical_fully_removed']:h[k]+=r.get(k,0)
    h.update(r.get('hist',{}))
    if len(seen)%5000==0:print(json.dumps({'replayed':len(seen),'seconds':time.time()-t}),flush=True)
  V2.check(count==expected['completed'] and dict(h)==expected['hist'],'phase accounting');total.update(h);phases.append({'records':count,'sha256':sha(folder/'certificates.jsonl.gz')})
 V2.check(seen==[i for i,x in enumerate(cover['cores']) if len(x)==1] and len(seen)==80914,'complete exact union')
 out={'status':'FULL_TWO_PHASE_INDEPENDENT_EXACT_REPLAY_PASSED','natural':len(cat),'valid':len(valid),'residual':len(seen),'phases':phases,'hist':dict(total),'seconds':time.time()-t};(ROOT/'complete_replay.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)
if __name__=='__main__':main()
