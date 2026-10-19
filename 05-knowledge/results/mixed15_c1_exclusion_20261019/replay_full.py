#!/usr/bin/env python3
"""Full independent stdlib replay; takes longer than provenance-only union audit."""
from pathlib import Path
import json,gzip,importlib.util,time
ROOT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('V',ROOT/'replay_c1.py');V=importlib.util.module_from_spec(spec);spec.loader.exec_module(V)

def main():
 start=time.time();data=json.loads(gzip.decompress((ROOT/'analytic_cover.json.gz').read_bytes()));cat=V.catalogue();V.check(len(cat)==96428,'catalogue')
 valid={p for p in cat if all(8-p[v].bit_count()-sum(bool(q>>v&1) for q in p)<=7 for v in range(9))};V.check(len(valid)==91604 and valid=={tuple(x['pred']) for x in data['cores']},'valid cores')
 for item in data['cores']:
  if 'balanced_pair' in item:V.R.replay(tuple(item['pred']),{'analytic_balanced':item['balanced_pair']})
  elif 'cycle_arcs' in item:V.R.replay(tuple(item['pred']),{'analytic_cycle':item['cycle_arcs']})
 seen=set();orders=0
 for path in [ROOT/'residual_certificates.jsonl.gz']+[ROOT/f'shard_{j}'/'certificates.jsonl.gz' for j in range(4)]:
  for line in gzip.decompress(path.read_bytes()).splitlines():
   rec=json.loads(line);i=rec['core_index'];V.check(i not in seen,'duplicate');seen.add(i);V.check(data['cores'][i]=={'pred':rec['pred']},'core mismatch');V.R.replay(tuple(rec['pred']),rec);orders+=rec.get('residual_orders',0)
  print(json.dumps({'source':str(path.relative_to(ROOT)),'verified_cores':len(seen),'verified_orders':orders}),flush=True)
 V.check(seen=={i for i,x in enumerate(data['cores']) if len(x)==1},'missing core')
 print(json.dumps({'status':'FULL_INDEPENDENT_REPLAY_PASSED','cores':len(seen),'orders':orders,'seconds':time.time()-start}),flush=True)
if __name__=='__main__':main()
