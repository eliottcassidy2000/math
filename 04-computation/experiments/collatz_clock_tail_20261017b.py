#!/usr/bin/env python3
"""Full finite-graph clock-potential audit and physical cone-tail receipts.
No ROOT oracle. Every observed transition is the literal shortcut T.
"""
import argparse,hashlib,json
from collections import Counter
from pathlib import Path
import numpy as np
import collatz_interval_atlas_20261016e as A
from collatz_compact_resets_20261015f import Graph,T
checks=0
def need(ok,label):
 global checks
 checks+=1
 if not ok:raise RuntimeError(label)

def audit(g):
 need(not g.seen[1] and int(g.endpoint[1])==1 and int(g.clock[1])==0,'stopped root')
 memo={};outside_steps=0
 def value(x):
  nonlocal outside_steps
  path=[];where=set()
  while x>g.M and g.has(x) and x not in memo:
   need(x not in where,'outside graph has no cycle');where.add(x);path.append(x);x=T(x)
  if x<=g.M:e,d=int(g.endpoint[x]),int(g.clock[x])
  elif x in memo:e,d=memo[x]
  else:e,d=x,0
  for x in reversed(path):d+=1;memo[x]=(e,d);outside_steps+=1
  return e,d
 # These are the whole observed outside edges, not selected random paths.
 for x in g.outside:value(x)
 edges=0
 for lo in range(1,g.M+1,1<<18):
  hi=min(g.M+1,lo+(1<<18));n=np.arange(lo,hi,dtype=np.uint64);seen=g.seen[lo:hi]
  need(bool(np.all(g.endpoint[lo:hi][~seen]==n[~seen])) and bool(np.all(g.clock[lo:hi][~seen]==0)), 'every unobserved interior terminal has zero clock and its own label')
  src=n[seen];y=np.where(src&1,(3*src+1)//2,src//2);inside=y<=g.M
  need(bool(np.all(g.endpoint[src[inside]]==g.endpoint[y[inside]])) and
       bool(np.all(g.clock[src[inside]].astype(np.uint64)==g.clock[y[inside]].astype(np.uint64)+1)),
       'all observed interior-to-interior edges have exact endpoint and clock increment')
  for s,t in zip(src[~inside],y[~inside]):
   e,d=value(int(t));need(int(g.endpoint[s])==e and int(g.clock[s])==d+1,'all boundary-crossing edges have exact endpoint and clock increment')
  edges+=len(src)
 need(outside_steps==len(g.outside),'every outside observed edge audited once')
 return dict(interior_edges=edges,outside_edges=outside_steps,original_sources=g.M,
             endpoint_sha256=hashlib.sha256(g.endpoint.tobytes()).hexdigest(),clock_sha256=hashlib.sha256(g.clock.tobytes()).hexdigest())

def tail(values):
 histogram=Counter(values);count=len(values);energy=sum(v*v*m for v,m in histogram.items());maximum=max(histogram)
 remaining=count;weak=(0,0,0);layer=0;previous=0
 for v,m in sorted(histogram.items()):
  # Tail at v-1 includes exactly the remaining values >=v; v^2 is the
  # matching discrete weak-L2 weight. This convention also works at v=0.
  score=v*v*remaining
  if score>weak[0]:weak=(score,v,remaining)
  layer+=(v*v-previous*previous)*remaining
  remaining-=m;previous=v
 need(layer==energy and remaining==0,'exact layer-cake identity')
 return dict(energy=energy,maximum=maximum,weak_numerator=weak[0],critical_height=weak[1],critical_rows=weak[2],histogram=sorted(histogram.items()))

def selection(X,values,b,stage):
 anchors=[];cost=0
 for lo in range(0,X,b):
  j=min(range(lo,lo+b),key=values.__getitem__);B=values[j];anchors.append([X+j,B,stage])
  cost+=sum((B+4*abs(i-j))**2 for i in range(lo,lo+b))
 cone=A.cone_envelope(X,anchors);need(all(v<=c for v,c in zip(values,cone)),'source-preserving anchor cones dominate supplied fine residuals')
 need(sum(c*c for c in cone)<=cost,'all-anchor cone improves assigned-bin invoice')
 return dict(bin_width=b,anchor_count=len(anchors),anchors=anchors,assigned_energy=cost,cone_tail=tail(cone))

def controls():
 g=Graph(512);audit(g)
 for _ in range(4):
  cols,inv=np.unique(g.endpoint[1:],return_inverse=True);g.extend(cols);g.close(cols,inv);audit(g);g.replay(range(1,513))
 old=int(g.clock[27]);g.clock[27]+=2
 try:audit(g)
 except RuntimeError:pass
 else:raise RuntimeError('corrupted full clock escaped audit')
 g.clock[27]=old
 old=int(g.endpoint[27]);g.endpoint[27]=old+2
 try:audit(g)
 except RuntimeError:pass
 else:raise RuntimeError('corrupted endpoint escaped audit')
 g.endpoint[27]=old;audit(g)
 for v in ([0,0,0],[0,1,2,3],[7]*12,[0,8,0,8,32]):
  t=tail(v);need(t['energy']==sum(x*x for x in v),'direct tail energy positive/hostile control')
 return dict(full_original_replay_sources=512,stages=4,corruptions_rejected=['clock','endpoint'])

def run(k):
 X=1<<k;base=k//2+2;targets={base,base+2,base+6};g=Graph(8*X);initial=audit(g);stages=[]
 for stage in range(1,max(targets)+1):
  cols,inv=np.unique(g.endpoint[1:],return_inverse=True);old=g.clock.copy();g.extend(cols);g.close(cols,inv)
  need(bool(np.all(g.clock>=old)),'every original clock preserved through extension');del old
  receipt=audit(g)
  if stage in targets:
   values=A.norms(X,g);b=1<<(k//2)
   records=[selection(X,values,max(1,b//q),stage) for q in (1,4)]
   stages.append(dict(stage=stage,audit=receipt,fine_tail=tail(values),selections=records,max_clock=int(g.clock.max())))
 return dict(k=k,X=X,initial=initial,stages=stages)

def inherited(paths):
 out=[]
 for path in paths:
  raw=Path(path).read_bytes();obj=json.loads(raw)
  for r in obj['runs']:
   values=A.cone_envelope(r['X'],r['anchors']);q=tail(values)
   need(q['energy']==r.get('native_cone_energy',r.get('cone_energy')),'native anchor projection agrees with saved cone energy')
   out.append(dict(k=r['k'],X=r['X'],tail=q,source_sha256=hashlib.sha256(raw).hexdigest()))
 return out

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--k',type=int,nargs='*',default=[]);p.add_argument('--receipt',nargs='*',default=[]);args=p.parse_args()
 control=controls();runs=[run(k) for k in args.k];old=inherited(args.receipt)
 print(json.dumps(dict(status='FINITE-EXACT; uniform weak-tail estimate OPEN',controls=control,runs=runs,inherited=old,checks=checks),sort_keys=True,separators=(',',':')))
