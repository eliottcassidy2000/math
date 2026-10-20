"""Independent cyclic orbit routing, explicit invariant and inverse checks."""
import json
from functools import lru_cache

def need(ok,msg):
 if not ok:raise RuntimeError(msg)

reports=[]
for m in range(1,15):
 mask=(1<<m)-1;e=m&-m
 def rot(x):return (x>>1)|((x&1)<<(m-1))
 @lru_cache(None)
 def canon(x):
  best=x;y=rot(x)
  for _ in range(m-1):best=min(best,y);y=rot(y)
  return best
 def D(x):return canon(x^rot(x))
 def stable(x):
  r=0
  while x:r^=x&((1<<e)-1);x>>=e
  return r==0
 def preds(y):
  if y.bit_count()%2:return []
  z=0;bit=0
  for i in range(m-1):
   bit^=(y>>i)&1;z|=bit<<(i+1)
  return sorted({canon(z),canon(z^mask)})
 reps=sorted({canon(x) for x in range(1<<m)})
 actual={x:[] for x in reps}
 for x in reps:actual[D(x)].append(x)
 for y in reps:need(preds(y)==sorted(actual[y]),'preimage solver')
 src=[x for x in reps if len(actual[x])==0];branch=[x for x in reps if len(actual[x])==2]
 need(src==[x for x in reps if x.bit_count()%2],'source parity')
 def period(x):
  y=rot(x);d=1
  while y!=x:y=rot(y);d+=1
  return d
 def evenroot(x):
  d=period(x);return ((x&((1<<d)-1)).bit_count()%2)==0
 need(branch==[x for x in reps if evenroot(x)],'branch primitive parity')
 def cut(y):
  pp=preds(y);need(len(pp)==2,'branch required')
  if stable(y):
   qq=[x for x in pp if not stable(x)];need(len(qq)==1,'one unstable core predecessor');return qq[0]
  return pp[-1]
 def forward(x):
  need(x in src,'source required');steps=0
  while True:
   y=D(x);steps+=1
   if len(preds(y))==2 and cut(y)==x:return y,steps
   x=y;need(steps<=e,'forward termination')
 def inverse(y):
  x=cut(y);steps=1
  while True:
   pp=[z for z in preds(x) if not(len(preds(x))==2 and cut(x)==z)]
   if not pp:return x,steps
   need(len(pp)==1,'unique remaining predecessor');x=pp[0];steps+=1;need(steps<=e,'inverse termination')
 pairs={};maxsteps=0
 for s in src:
  b,steps=forward(s);back,isteps=inverse(b)
  need(back==s and steps==isteps and steps<=e,'inverse or depth')
  ds,db=period(s),period(b)
  need(ds//(ds&-ds)==db//(db&-db),'odd periods')
  pairs[s]=b;maxsteps=max(maxsteps,steps)
 need(set(pairs.values())==set(branch) and len(pairs)==len(branch),'bijection')
 reports.append(dict(m=m,nilpotent_depth=e,necklaces=len(reps),sources=len(src),branches=len(branch),max_steps=maxsteps))
print(json.dumps(dict(status='PASS',reports=reports),indent=2))
