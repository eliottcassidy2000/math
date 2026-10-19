#!/usr/bin/env python3
"""Fast union/provenance audit. Mathematical replay lives in replay_complete.py."""
from pathlib import Path
from collections import Counter
import json,gzip,hashlib
ROOT=Path(__file__).resolve().parent

def check(x,msg):
 if not x:raise RuntimeError(msg)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 state=json.loads((ROOT/'progress.json').read_text());cover=json.loads(gzip.decompress((ROOT/'analytic_cover.json.gz').read_bytes()));proof=json.loads((ROOT/'independent_cover.json').read_text());check(state['status']=='ALL_CANONICAL_CORES_EXACT_CERTIFIED_AND_REPLAYED','not complete')
 for name,digest in state['dependencies'].items():check(sha(ROOT/name)==digest,'source changed '+name)
 check(sha(ROOT/'analytic_cover.json.gz')==proof['cover_sha256'],'cover identity');check(proof['status']=='COMPLETE_INDEPENDENT_CATALOGUE_AND_ANALYTIC_REPLAY_PASSED','cover replay');check(sha(ROOT/'certificates.jsonl.gz')==state['record_sha256'],'record hash')
 ids=[];h=Counter();max_weight=certs=0
 with gzip.open(ROOT/'certificates.jsonl.gz','rt') as f:
  for line in f:
   rec=json.loads(line);i=rec['core_index'];ids.append(i);check(cover['cores'][i]=={'pred':rec['pred']},'core identity')
   for k in ['residual_orders','analytic_directed_orders','short_span_orders','canonical_excluded_orders','canonical_fully_removed']:h[k]+=rec.get(k,0)
   h.update(rec.get('hist',{}));certs+=len(rec.get('certificates',[]))
   for cert in rec.get('certificates',[]):
    for a,b,w in cert['terms']:check(type(w) is int and w>0,'integer weight');max_weight=max(max_weight,w)
 check(ids==[i for i,item in enumerate(cover['cores']) if len(item)==1],'coverage');check(len(ids)==state['completed']==70103,'count');check(dict(h)==state['hist'],'hist')
 out={'status':'COMPLETE_CANONICAL_COVERAGE_AND_REPLAY_PROVENANCE_PASSED','degree_valid':91604,'analytic_balanced':21501,'residual_cores':len(ids),'hist':dict(h),'maximum_weight':max_weight,'distinct_certificates':certs,'record_sha256':state['record_sha256']};(ROOT/'audit_receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)
if __name__=='__main__':main()
