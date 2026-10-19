from pathlib import Path
from functools import lru_cache
from collections import Counter
import json,gzip,hashlib,time
import cone_engine as F
import exact_replay_engine as R
from catalogue import catalogue
ROOT=Path(__file__).resolve().parent

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
 data=json.loads(gzip.decompress((ROOT/'analytic_cover.json.gz').read_bytes()));dep={name:sha(ROOT/name) for name in ['run.py','cone_engine.py','exact_replay_engine.py','catalogue.py','analytic_cover.json.gz']};progress=ROOT/'progress.json';records=ROOT/'certificates.jsonl.gz'
 if progress.exists():
  state=json.loads(progress.read_text());F.require(state['dependencies']==dep,'sources changed');F.require(state['status']=='RUNNING','terminal');
  with records.open('r+b') as f:f.truncate(state['committed_bytes'])
 else:
  cat=catalogue();valid={p for p in cat if max(F.degrees(p))<=7};F.require(len(cat)==96428 and len(valid)==91604,'catalogue counts');F.require(valid=={tuple(item['pred']) for item in data['cores']},'independent catalogue match')
  h=Counter()
  for item in data['cores']:
   if 'balanced_pair' in item:R.replay(tuple(item['pred']),{'analytic_balanced':item['balanced_pair']});h['balanced']+=1
   elif 'cycle_arcs' in item:R.replay(tuple(item['pred']),{'analytic_cycle':item['cycle_arcs']});h['cycle']+=1
   else:h['residual']+=1
  F.require(h['balanced']==21501 and h['residual']==70103 and h['cycle']==0,'analytic counts')
  write(ROOT/'independent_cover.json',{'status':'COMPLETE_INDEPENDENT_CATALOGUE_AND_ANALYTIC_REPLAY_PASSED','natural':len(cat),'valid':len(valid),'hist':dict(h),'cover_sha256':dep['analytic_cover.json.gz']})
  state={'status':'RUNNING','scope':'canonical(maximum span,maximum common count), mixed(2,4,1)','dependencies':dep,'target':70103,'completed':0,'last_core_index':-1,'committed_bytes':0,'hist':{},'start_unix':time.time()};records.write_bytes(b'');write(progress,state)
 for i,item in enumerate(data['cores']):
  if i<=state['last_core_index'] or len(item)>1:continue
  p=tuple(item['pred']);r=F.scan(p);r.update(pred=p,core_index=i)
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
