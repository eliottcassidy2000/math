from pathlib import Path
from collections import Counter
import time,json,gzip,hashlib
import analytic_census as A
import exact_replay_engine as R
ROOT=Path(__file__).resolve().parent
for scope in [(2,4,0),(3,3,0)]:
 label='%s%s%s'%scope;data=A.prepare(scope);start=time.time();hist=Counter();records=[];status='COMPLETE_EXACT_REPLAY_PASSED'
 for i,item in enumerate(data['cores']):
  p=tuple(item['pred'])
  if 'balanced_pair' in item:r={'analytic_balanced':item['balanced_pair']}
  elif 'cycle_arcs' in item:r={'analytic_cycle':item['cycle_arcs']}
  else:r=A.F.scan(p)
  r.update(pred=p,core_index=i)
  if 'obstruction' in r:
   (ROOT/f'full_obstacle_{label}.json').write_text(json.dumps(r,indent=2));status='EXACT_CONE_OBSTACLE';break
  R.replay(p,r);records.append(r)
  hist['cores']+=1;hist['balanced']+=int('analytic_balanced' in r);hist['cycle']+=int('analytic_cycle' in r);hist['residual_orders']+=r.get('residual_orders',0);hist['directed_orders']+=r.get('analytic_directed_orders',0);hist['short_span_orders']+=r.get('short_span_orders',0)
  for k,v in r.get('hist',{}).items():hist[k]+=v
  if i%100==0:
   (ROOT/f'full_progress_{label}.json').write_text(json.dumps({'scope':scope,'status':'RUNNING','hist':dict(hist),'seconds':time.time()-start},indent=2))
   print(json.dumps({'scope':scope,'hist':dict(hist),'seconds':time.time()-start}),flush=True)
 file=ROOT/f'full_certificates_{label}.json.gz';file.write_bytes(gzip.compress(json.dumps(records,separators=(',',':')).encode(),mtime=0));out={'scope':scope,'status':status,'hist':dict(hist),'seconds':time.time()-start,'record_sha256':hashlib.sha256(file.read_bytes()).hexdigest()};(ROOT/f'full_progress_{label}.json').write_text(json.dumps(out,indent=2));print(json.dumps(out),flush=True)
