#!/usr/bin/env python3
"""Fresh complete stdlib replay, independent of the numerical search."""
from pathlib import Path
from collections import Counter
import json,gzip,hashlib,time
from catalogue import catalogue
import exact_replay_engine as R
ROOT=Path(__file__).resolve().parent

def main():
 start=time.time();cover=json.loads(gzip.decompress((ROOT/'analytic_cover.json.gz').read_bytes()));cat=catalogue();valid={p for p in cat if all(8-p[v].bit_count()-sum(bool(q>>v&1) for q in p)<=7 for v in range(9))};R.check(len(cat)==96428 and len(valid)==91604,'catalogue counts');R.check(valid=={tuple(x['pred']) for x in cover['cores']},'catalogue equality')
 for item in cover['cores']:
  if 'balanced_pair' in item:R.replay(tuple(item['pred']),{'analytic_balanced':item['balanced_pair']})
  elif 'cycle_arcs' in item:R.replay(tuple(item['pred']),{'analytic_cycle':item['cycle_arcs']})
 seen=[];hist=Counter();last=-1
 with gzip.open(ROOT/'certificates.jsonl.gz','rt') as f:
  for line in f:
   rec=json.loads(line);i=rec['core_index'];R.check(i>last,'record order');last=i;R.check(cover['cores'][i]=={'pred':rec['pred']},'record identity');R.replay(tuple(rec['pred']),rec);seen.append(i)
   for k in ['residual_orders','analytic_directed_orders','short_span_orders','canonical_excluded_orders','canonical_fully_removed']:hist[k]+=rec.get(k,0)
   hist.update(rec.get('hist',{}))
   if len(seen)%5000==0:print(json.dumps({'replayed':len(seen),'hist':dict(hist),'seconds':time.time()-start}),flush=True)
 R.check(seen==[i for i,x in enumerate(cover['cores']) if len(x)==1],'full coverage');scan=json.loads((ROOT/'progress.json').read_text());R.check(scan['status']=='ALL_CANONICAL_CORES_EXACT_CERTIFIED_AND_REPLAYED' and scan['completed']==len(seen),'scan terminal');R.check(scan['hist']==dict(hist),'histogram');digest=hashlib.sha256((ROOT/'certificates.jsonl.gz').read_bytes()).hexdigest();R.check(digest==scan['record_sha256'],'record hash')
 out={'status':'FULL_INDEPENDENT_CANONICAL_REPLAY_PASSED','natural_cores':len(cat),'degree_valid_cores':len(valid),'residual_cores':len(seen),'hist':dict(hist),'record_sha256':digest,'seconds':time.time()-start};(ROOT/'complete_replay.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)
if __name__=='__main__':main()
