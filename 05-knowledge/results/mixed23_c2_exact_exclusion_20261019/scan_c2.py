#!/usr/bin/env python3
"""Checkpointed mixed(2,3)c2 morphology scan; stop at first exact cone obstacle."""
from pathlib import Path
from itertools import combinations
from functools import lru_cache
from collections import Counter
import json,gzip,hashlib,time,importlib.util
ROOT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('engine',ROOT/'cone_engine.py');F=importlib.util.module_from_spec(spec);spec.loader.exec_module(F)

def write(path,value):
 t=path.with_suffix(path.suffix+'.tmp');t.write_text(json.dumps(value,indent=2)+'\n');t.replace(path)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def natural_ideals(n):
 level=[()]
 for size in range(n):
  nxt=[]
  for p in level:
   for mask in range(1<<size):
    if all(not(mask>>i&1) or not(p[i]&~mask) for i in range(size)):nxt.append(p+(mask,))
  level=nxt
 return level

@lru_cache(maxsize=1)
def all_bounds(p):
 n=len(p);up=tuple(sum(1<<j for j in range(n) if p[j]>>i&1) for i in range(n));deg=F.degrees(p);seed=3|sum(1<<i for i,d in enumerate(deg) if d==7)
 alpha=[0 if ((1<<v)|p[v])&seed else 7-deg[v] for v in range(n)]
 beta=[0 if ((1<<v)|up[v])&seed else 7-deg[v] for v in range(n)]
 return {(u,v):(p[v]&~p[u]).bit_count()+(up[u]&~up[v]).bit_count()+alpha[u]+beta[v] for u in range(n) for v in range(n) if u!=v and not(p[v]>>u&1 or p[u]>>v&1)}
F.defect_bound=lambda p,u,v:all_bounds(p)[u,v]

def prepare():
 target=ROOT/'analytic_cover.json.gz'
 if target.exists():return json.loads(gzip.decompress(target.read_bytes()))
 start=time.time();nonroot=natural_ideals(7);F.require(len(nonroot)==96428,'natural count');kept=[];discarded=0
 for q in nonroot:
  p=[0,12]+[v<<2 for v in q]
  for a in (6,7,8):p[a]|=1
  p=tuple(p);F.require(F.closure(p)==p,'closure')
  if max(F.degrees(p))>7:discarded+=1;continue
  kept.append(p)
 F.require(len(kept)==87137 and discarded==9291,'degree-valid count');kept.sort();F.require(len(set(kept))==len(kept),'duplicates')
 data={'status':'ANALYTIC_COVER_COMPLETE','natural_nonroot_count':len(nonroot),'degree_invalid':discarded,'degree_valid':len(kept),'cores':[],'hist':Counter()}
 for i,p in enumerate(kept):
  bounds=all_bounds(p);wit=next(((a,b) for (a,b),h in bounds.items() if a<b and h<=1 and bounds[b,a]<=1),None)
  if wit:record={'pred':p,'balanced_pair':wit};data['hist']['balanced_pair']+=1
  else:
   cyc=F.forced_cycle(p,bounds)
   if cyc:record={'pred':p,'cycle_arcs':cyc};data['hist']['cycle']+=1
   else:record={'pred':p};data['hist']['residual']+=1
  data['cores'].append(record)
  if i%10000==0:print(json.dumps({'phase':'analytic','processed':i+1,'hist':data['hist'],'seconds':time.time()-start}),flush=True)
 data['seconds']=time.time()-start;target.write_bytes(gzip.compress(json.dumps(data,separators=(',',':')).encode(),mtime=0));print(json.dumps({k:v for k,v in data.items() if k!='cores'}),flush=True);return data

def main():
 data=prepare();deps={p.name:sha(p) for p in (Path(__file__),ROOT/'cone_engine.py',ROOT/'analytic_cover.json.gz')};progress=ROOT/'progress.json';records=ROOT/'residual_certificates.jsonl.gz'
 if progress.exists():
  state=json.loads(progress.read_text());F.require(state['dependencies']==deps,'changed dependency');F.require(state['status']=='RUNNING','terminal state')
  with records.open('r+b') as h:h.truncate(state['committed_bytes'])
 else:
  state={'status':'RUNNING','scope':'Partial coverage only until terminal success; exact cone obstacles are relaxation limitations, not actual posets.','dependencies':deps,'analytic_covered':data['hist'].get('balanced_pair',0)+data['hist'].get('cycle',0),'residual_target':data['hist']['residual'],'completed_residual_cores':0,'last_core_index':-1,'committed_bytes':0,'residual_orders':0,'directed_screen_orders':0,'hist':{},'start_unix':time.time()};records.write_bytes(b'');write(progress,state)
 for i,item in enumerate(data['cores']):
  if i<=state['last_core_index'] or len(item)>1:continue
  p=tuple(item['pred']);rec=F.scan(p);rec['core_index']=i;rec['pred']=p
  if 'obstruction' in rec:
   write(ROOT/'first_obstacle.json',rec);state['status']='EXACT_CONE_OBSTACLE';state['first_obstacle_core']=i;write(progress,state);print(json.dumps({'status':state['status'],'core':i}),flush=True);return
  payload=gzip.compress((json.dumps(rec,separators=(',',':'))+'\n').encode(),mtime=0)
  with records.open('ab') as h:h.write(payload);h.flush()
  state['committed_bytes']+=len(payload);state['completed_residual_cores']+=1;state['last_core_index']=i;state['residual_orders']+=rec.get('residual_orders',0);state['directed_screen_orders']+=rec.get('analytic_directed_orders',0);hist=Counter(state['hist']);hist.update(rec.get('hist',{}));state['hist']=dict(hist);state['elapsed']=time.time()-state['start_unix'];write(progress,state)
  if state['completed_residual_cores']%100==0:print(json.dumps(state),flush=True)
 state['status']='ALL_CORES_CERTIFIED';state['record_sha256']=sha(records);write(progress,state);print(json.dumps(state),flush=True)
if __name__=='__main__':main()
