#!/usr/bin/env python3
"""Independent standard-library cover and streaming certificate replay."""
from pathlib import Path
from itertools import product,combinations
from collections import Counter
import json,gzip,hashlib,importlib.util,time
ROOT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('R',ROOT/'exact_replay_engine.py');R=importlib.util.module_from_spec(spec);spec.loader.exec_module(R)
def check(ok,msg):
 if not ok:raise RuntimeError(msg)
def write(path,val):
 t=path.with_suffix(path.suffix+'.tmp');t.write_text(json.dumps(val,indent=2)+'\n');t.replace(path)
def small_posets():
 out=set();edges=[(0,1),(0,2),(1,2)]
 for bits in range(8):
  rel={e for i,e in enumerate(edges) if bits>>i&1}
  if (0,1) in rel and (1,2) in rel:rel.add((0,2))
  out.add(tuple(sum(1<<a for a,b in rel if b==j) for j in range(3)))
 return sorted(out)
def catalogue():
 out=set();small=small_posets()
 for zp,ap in product(small,repeat=2):
  ideals=[m for m in range(8) if all(not(m>>i&1) or not(zp[i]&~m) for i in range(3))]
  upsets=[7^m for m in ideals]
  aupsets=[m for m in range(8) if all(not(m>>i&1) or m>>j&1 for j in range(3) for i in range(3) if ap[j]>>i&1)]
  for lower,d0,d1,d2,ba in product(upsets,ideals,ideals,ideals,aupsets):
   ds=(d0,d1,d2)
   if any(ds[i]&~ds[j] for j in range(3) for i in range(3) if ap[j]>>i&1):continue
   if any(lower&ds[i] and not(ba>>i&1) for i in range(3)):continue
   p=(0,4,0)+tuple((zp[i]<<3)|(4 if lower>>i&1 else 0) for i in range(3))+tuple(1|(ds[i]<<3)|(ap[i]<<6)|(4 if ba>>i&1 else 0) for i in range(3))
   out.add(p)
 return out

def main():
 start=time.time();raw=(ROOT/'analytic_cover.json.gz').read_bytes();data=json.loads(gzip.decompress(raw));cover_sha=hashlib.sha256(raw).hexdigest();target=ROOT/'independent_replay.json'
 if target.exists():
  state=json.loads(target.read_text());check(state['cover_sha256']==cover_sha,'cover changed')
 else:
  cat=catalogue();check(len(cat)==96428,'independent total');valid={p for p in cat if all(8-p[v].bit_count()-sum(bool(q>>v&1) for q in p)<=7 for v in range(9))};check(len(valid)==82987,'independent valid');check(valid=={tuple(r['pred']) for r in data['cores']},'catalogue mismatch')
  hist=Counter()
  for item in data['cores']:
   p=tuple(item['pred'])
   if 'balanced_pair' in item:R.replay(p,{'analytic_balanced':item['balanced_pair']});hist['balanced_pair']+=1
   elif 'cycle_arcs' in item:R.replay(p,{'analytic_cycle':item['cycle_arcs']});hist['cycle']+=1
   else:hist['residual']+=1
  check(dict(hist)==data['hist'],'cover histogram')
  state={'status':'COVER_VERIFIED_REPLAY_RUNNING','cover_sha256':cover_sha,'cover_hist':dict(hist),'independent_catalogue_total':len(cat),'independent_degree_valid':len(valid),'verified_residual_cores':0,'verified_residual_orders':0,'verified_bytes':0,'last_core_index':-1,'start_unix':time.time()};write(target,state);print(json.dumps(state),flush=True)
 while True:
  scan=json.loads((ROOT/'progress.json').read_text());end=scan['committed_bytes']
  if end>state['verified_bytes']:
   with (ROOT/'residual_certificates.jsonl.gz').open('rb') as f:f.seek(state['verified_bytes']);payload=f.read(end-state['verified_bytes'])
   records=[json.loads(line) for line in gzip.decompress(payload).splitlines()]
   for rec in records:
    i=rec['core_index'];check(i>state['last_core_index'],'record order');check(data['cores'][i]=={'pred':rec['pred']},'record catalogue mismatch');R.replay(tuple(rec['pred']),rec);state['verified_residual_cores']+=1;state['verified_residual_orders']+=rec.get('residual_orders',0);state['last_core_index']=i
   state['verified_bytes']=end;state['elapsed']=time.time()-state['start_unix'];write(target,state);print(json.dumps(state),flush=True)
  elif scan['status']!='RUNNING':
   if scan['status']=='ALL_CORES_CERTIFIED':
    check(state['verified_residual_cores']==data['hist']['residual'],'coverage missing');state['status']='ALL_EXACT_REPLAY_PASSED';state['certificate_sha256']=hashlib.sha256((ROOT/'residual_certificates.jsonl.gz').read_bytes()).hexdigest()
   else:state['status']='SCAN_STOPPED_PARTIAL_REPLAY_PASSED'
   write(target,state);print(json.dumps(state),flush=True);return
  else:time.sleep(5)
if __name__=='__main__':main()
