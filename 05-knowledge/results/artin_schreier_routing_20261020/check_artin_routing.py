import json
MODS={1:3,2:7,3:11,4:19,5:37,6:67,7:131,8:285}
checks=0;rows=[]
def need(ok,msg):
 global checks
 checks+=1
 if not ok:raise RuntimeError(msg)
for m,mod in MODS.items():
 q=1<<m;e=m&-m
 def mul(a,b):
  z=0
  while b:
   if b&1:z^=a
   b>>=1;a<<=1
   if a&q:a^=mod
  return z
 def power(a,k):
  z=1
  while k:
   if k&1:z=mul(z,a)
   a=mul(a,a);k>>=1
  return z
 def sq(a):return mul(a,a)
 def trace(a,d,step=1):
  z=0
  for _ in range(d):z^=a;a=power(a,1<<step)
  return z
 def orb(a):
  out=[a];b=sq(a)
  while b!=a:out.append(b);b=sq(b);need(len(out)<=m,'orbitfinite')
  return out
 for a in range(1,q):need(power(a,q-1)==1,'irreducible modulus field test')
 orbits={min(orb(a)):orb(a) for a in range(q)};canon={a:r for r,o in orbits.items() for a in o}
 poly={}
 for r,o in orbits.items():
  co=[1]
  for a in o:
   nxt=[0]*(len(co)+1)
   for i,b in enumerate(co):nxt[i]^=mul(a,b);nxt[i+1]^=b
   co=nxt
  need(all(c in (0,1) for c in co),'minimal polynomial binary')
  poly[r]=sum(c<<i for i,c in enumerate(co))
 def D(r):return canon[sq(r)^r]
 def stable(r):return trace(r,m//e,e)==0
 pred={r:[] for r in orbits}
 for r in orbits:pred[D(r)].append(r)
 src=[r for r in orbits if not pred[r]];br=[r for r in orbits if len(pred[r])==2]
 for r,o in orbits.items():
  d=len(o);t=trace(r,d);ambient=trace(r,m)
  need(t in (0,1) and ambient in (0,1),'binary trace')
  need((r in src)==(ambient==1),'source trace')
  need((r in br)==(t==0),'branch trace')
  need(len(pred[r])<=2,'indegree')
  v=0;period=0
  while True:
   v=sq(v)^r;period+=1
   if v==0:break
   need(period<=2*d,'critical termination')
  need(period==d*(1+t),'critical period formula')
  y=r
  for _ in range(e):y=D(y)
  need(stable(y),'Fitting bound')
 def cut(y):
  need(len(pred[y])==2,'branch')
  if stable(y):
   off=[r for r in pred[y] if not stable(r)];need(len(off)==1,'unique offcore');return off[0]
  return max(pred[y],key=lambda r:poly[r])
 images={};maxsteps=0
 for s in src:
  v=s;steps=0
  while True:
   y=D(v);steps+=1
   if len(pred[y])==2 and cut(y)==v:break
   v=y;need(steps<=e,'forwardbound')
  images[s]=y;maxsteps=max(maxsteps,steps)
  v=cut(y);backsteps=1
  while pred[v]:
   pp=[w for w in pred[v] if not(len(pred[v])==2 and cut(v)==w)]
   need(len(pp)==1,'unique inverse');v=pp[0];backsteps+=1;need(backsteps<=e,'inversebound')
  need(v==s and steps==backsteps,'inverse')
  ds,dt=len(orbits[s]),len(orbits[y]);need(ds//(ds&-ds)==dt//(dt&-dt),'odddegree')
 need(set(images.values())==set(br) and len(images)==len(br),'bijection')
 rows.append(dict(m=m,modulus=mod,orbits=len(orbits),source_count=len(src),branch_count=len(br),depth=e,max_steps=maxsteps,pairs=[dict(source_polynomial=poly[s],target_polynomial=poly[t]) for s,t in sorted(images.items())]))
print(json.dumps(dict(status='PASS',checks=checks,fields=rows),indent=2))
