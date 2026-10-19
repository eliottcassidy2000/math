#!/usr/bin/env python3
"""Check full deterministic coverage and authenticated per-record replay receipts.
This is a union/provenance audit; run replay_full.py to repeat mathematical replay.
"""
from pathlib import Path
from collections import Counter
import json,gzip,hashlib
ROOT=Path(__file__).resolve().parent

def check(x,msg):
 if not x:raise RuntimeError(msg)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def lines(p):return [json.loads(t) for t in gzip.decompress(p.read_bytes()).splitlines()]
def main():
 cover=json.loads(gzip.decompress((ROOT/'analytic_cover.json.gz').read_bytes()));plan=json.loads((ROOT/'shard_plan.json').read_text());prefix=json.loads((ROOT/'progress.json').read_text());replay=json.loads((ROOT/'independent_replay.json').read_text())
 check(sha(ROOT/'analytic_cover.json.gz')==plan['cover_sha256']==replay['cover_sha256'],'cover identity')
 check(replay['status']=='SCAN_STOPPED_PARTIAL_REPLAY_PASSED' and replay['verified_residual_cores']==plan['prefix_completed'],'prefix replay incomplete')
 check(sha(ROOT/'residual_certificates.jsonl.gz')==plan['prefix_record_sha256'],'prefix identity');records=lines(ROOT/'residual_certificates.jsonl.gz');check(len(records)==plan['prefix_completed'],'prefix count')
 residual=[i for i,item in enumerate(cover['cores']) if len(item)==1];check([r['core_index'] for r in records]==residual[:len(records)],'prefix IDs')
 check(plan['shards']==[[i for pos,i in enumerate(residual[len(records):]) if pos%4==j] for j in range(4)],'shard partition')
 for j,ids in enumerate(plan['shards']):
  folder=ROOT/f'shard_{j}';state=json.loads((folder/'progress.json').read_text());check(state['status']=='ALL_SHARD_CERTIFICATES_REPLAYED','shard incomplete');check(state['completed']==len(ids),'shard count')
  for name,h in state['dependencies'].items():check(sha(ROOT/name)==h,'changed source '+name)
  check(sha(folder/'certificates.jsonl.gz')==state['record_sha256'],'record hash');rr=lines(folder/'certificates.jsonl.gz');check([r['core_index'] for r in rr]==ids,'shard IDs');records.extend(rr)
 records.sort(key=lambda r:r['core_index']);check([r['core_index'] for r in records]==residual,'union coverage')
 for rec in records:check(cover['cores'][rec['core_index']]=={'pred':rec['pred']},'core record')
 hist=Counter();orders=directed=short=0;weights=[];distinct=0
 for rec in records:
  short+=rec.get('short_span_orders',0);hist.update(rec.get('hist',{}));orders+=rec.get('residual_orders',0);directed+=rec.get('analytic_directed_orders',0);distinct+=len(rec.get('certificates',[]));weights.extend(w for cert in rec.get('certificates',[]) for a,b,w in cert['terms'])
 check(all(type(w) is int and w>0 for w in weights),'certificate weight is not a positive integer')
 out={'status':'COMPLETE_EXACT_CERTIFIED_AND_INDEPENDENTLY_REPLAYED_UNION','core_total':cover['natural_nonroot_count'],'degree_valid':cover['degree_valid'],'analytic_hist':cover['hist'],'residual_cores':len(records),'residual_orders':orders,'directed_screen_orders':directed,'known_six_thin_short_span_orders':short,'certificate_order_kinds':dict(hist),'distinct_certificates':distinct,'maximum_integer_weight':max(weights,default=0),'cover_sha256':sha(ROOT/'analytic_cover.json.gz'),'prefix_sha256':sha(ROOT/'residual_certificates.jsonl.gz'),'shard_record_sha256':[sha(ROOT/f'shard_{j}'/'certificates.jsonl.gz') for j in range(4)]}
 (ROOT/'union_receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)
if __name__=='__main__':main()
