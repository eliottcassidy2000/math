from pathlib import Path
from functools import lru_cache
from collections import Counter
import json,gzip,hashlib,time
import cone_engine as F
import exact_replay_engine as R
from catalogue import catalogue
ROOT=Path(__file__).resolve().parent
OLD=ROOT.parent/'canonical331-stage'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,d):
 t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(d,indent=2)+'\n');t.replace(p)
@lru_cache(maxsize=1)
def bounds(p):
 n=len(p);up=tuple(sum(1<<j for j in range(n) if p[j]>>i&1) for i in range(n));deg=F.degrees(p);seed=3|sum(1<<i for i,d in enumerate(deg) if d==7)
 alpha=[0 if ((1<<v)|p[v])&seed else 7-deg[v] for v in range(n)];beta=[0 if ((1<<v)|up[v])&seed else 7-deg[v] for v in range(n)]
 return {(u,v):(p[v]&~p[u]).bit_count()+(up[u]&~up[v]).bit_count()+alpha[u]+beta[v] for u in range(n) for v in range(n) if u!=v and not(p[v]>>u&1 or p[u]>>v&1)}
F.defect_bound=lambda p,u,v:bounds(p)[u,v]
def main():
 data=json.loads(gzip.decompress((OLD/'analytic_cover.json.gz').read_bytes()))
 prior=ROOT.parent/'canonical331-continuation-stage';old=json.loads((prior/'progress.json').read_text())
 F.require(old['status']=='EXACT_CONE_OBSTACLE' and old['completed']==2876,'phase-two stop')
 dep={name:sha(ROOT/name) for name in ['run.py','cone_engine.py','exact_replay_engine.py','catalogue.py','source_lock.json','replay_complete.py']}
 locked=json.loads((ROOT/'source_lock.json').read_text())
 for name,digest in locked.items():F.require(sha(ROOT.parent/name)==digest,'immutable prior source changed')
 progress=ROOT/'progress.json';records=ROOT/'certificates.jsonl.gz'
 if progress.exists():
  state=json.loads(progress.read_text());F.require(state['dependencies']==dep,'source changed');F.require(state['status']=='RUNNING','terminal')
  with records.open('r+b') as f:f.truncate(state['committed_bytes'])
 else:
  state={'schema':'139-ancestor-caps-v3','status':'RUNNING','scope':'canonical(maximum span,maximum common count), mixed(3,3,1)','dependencies':dep,'prior_completed':57231,'target':23683,'completed':0,'last_core_index':old['last_core_index'],'committed_bytes':0,'hist':{},'start_unix':time.time()};records.write_bytes(b'');write(progress,state)
 for i,item in enumerate(data['cores']):
  if i<=state['last_core_index'] or len(item)>1:continue
  p=tuple(item['pred']);r=F.scan(p);r.update(pred=p,core_index=i,schema='139-ancestor-caps-v3')
  if 'obstruction' in r:
   write(ROOT/'first_obstacle.json',r);state['status']='EXACT_CONE_OBSTACLE';state['obstacle_core']=i;write(progress,state);print(json.dumps(state),flush=True);return
  R.replay(p,r);payload=gzip.compress((json.dumps(r,separators=(',',':'))+'\n').encode(),mtime=0)
  with records.open('ab') as f:f.write(payload);f.flush()
  state['committed_bytes']+=len(payload);state['completed']+=1;state['last_core_index']=i;h=Counter(state['hist'])
  for k in ['residual_orders','analytic_directed_orders','short_span_orders','canonical_excluded_orders','canonical_fully_removed']:h[k]+=r.get(k,0)
  h.update(r.get('hist',{}));state['hist']=dict(h);state['seconds']=time.time()-state['start_unix'];write(progress,state)
  if state['completed']%100==0:print(json.dumps(state),flush=True)
 state['status']='ALL_CANONICAL_CORES_EXACT_CERTIFIED_AND_REPLAYED';state['record_sha256']=sha(records);write(progress,state);print(json.dumps(state),flush=True)
if __name__=='__main__':main()
