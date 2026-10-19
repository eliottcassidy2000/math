#!/usr/bin/env python3
"""Deterministic residual shard; each record independently replayed before commit."""
from pathlib import Path
import json,gzip,hashlib,importlib.util,time,sys
from collections import Counter
ROOT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('S',ROOT/'scan_c2.py');S=importlib.util.module_from_spec(spec);spec.loader.exec_module(S)
spec=importlib.util.spec_from_file_location('R',ROOT/'exact_replay_engine.py');R=importlib.util.module_from_spec(spec);spec.loader.exec_module(R)
F=S.F

def main():
 shard=int(sys.argv[1]);plan=json.loads((ROOT/'shard_plan.json').read_text());ids=plan['shards'][shard];data=json.loads(gzip.decompress((ROOT/'analytic_cover.json.gz').read_bytes()));out=ROOT/f'shard_{shard}';out.mkdir(exist_ok=True);progress=out/'progress.json';records=out/'certificates.jsonl.gz'
 deps={p.name:S.sha(p) for p in (Path(__file__),ROOT/'scan_c2.py',ROOT/'cone_engine.py',ROOT/'exact_replay_engine.py',ROOT/'shard_plan.json',ROOT/'analytic_cover.json.gz')}
 if progress.exists():
  state=json.loads(progress.read_text());F.require(state['dependencies']==deps,'dependencies');F.require(state['status']=='RUNNING','terminal')
  with records.open('r+b') as f:f.truncate(state['committed_bytes'])
 else:
  state={'status':'RUNNING','shard':shard,'target':len(ids),'completed':0,'last_position':-1,'last_core_index':None,'orders':0,'screen_orders':0,'hist':{},'committed_bytes':0,'dependencies':deps,'start_unix':time.time()};records.write_bytes(b'');S.write(progress,state)
 for pos in range(state['last_position']+1,len(ids)):
  i=ids[pos];p=tuple(data['cores'][i]['pred']);F.require(len(data['cores'][i])==1,'not residual');rec=F.scan(p);rec['core_index']=i;rec['pred']=p
  if 'obstruction' in rec:
   S.write(out/'first_obstacle.json',rec);state['status']='EXACT_CONE_OBSTACLE';state['obstacle_core']=i;S.write(progress,state);print(json.dumps(state),flush=True);return
  R.replay(p,rec)
  payload=gzip.compress((json.dumps(rec,separators=(',',':'))+'\n').encode(),mtime=0)
  with records.open('ab') as f:f.write(payload);f.flush()
  state['committed_bytes']+=len(payload);state['completed']+=1;state['last_position']=pos;state['last_core_index']=i;state['orders']+=rec.get('residual_orders',0);state['screen_orders']+=rec.get('analytic_directed_orders',0);h=Counter(state['hist']);h.update(rec.get('hist',{}));state['hist']=dict(h);state['elapsed']=time.time()-state['start_unix'];S.write(progress,state)
  if state['completed']%100==0:print(json.dumps(state),flush=True)
 state['status']='ALL_SHARD_CERTIFICATES_REPLAYED';state['record_sha256']=S.sha(records);S.write(progress,state);print(json.dumps(state),flush=True)
if __name__=='__main__':main()
