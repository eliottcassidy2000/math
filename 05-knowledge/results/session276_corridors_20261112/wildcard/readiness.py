"""All-order append readiness, preserving the actual capped continuation state."""
from capped import build,minimize,count
import collections,json
N=127
vs,ts,ws,hs=build(N)
mt,mh,part=minimize(ts,hs)
rep={}
for i,p in enumerate(part[:-1]):rep.setdefault(p,ws[i])
adj=[set(j for j in row if j!=i) for i,row in enumerate(mt)]
indeg=[0]*len(mt)
for r in adj:
 for j in r:indeg[j]+=1
q=collections.deque(i for i,d in enumerate(indeg) if not d);order=[]
while q:
 i=q.popleft();order.append(i)
 for j in adj[i]:
  indeg[j]-=1
  if not indeg[j]:q.append(j)
if len(order)!=len(mt):raise RuntimeError('Moore nontrivial SCC')
reach=[0]*len(mt)
for i in reversed(order):
 reach[i]=1<<i
 for j in adj[i]:reach[i]|=reach[j]

def nodes(mask):
 while mask:
  b=mask&-mask;yield b.bit_length()-1;mask-=b
by=collections.defaultdict(list)
for i,h in enumerate(mh):by[h].append(i)
# winning[d] bit mask: states supporting complete binary +/-1 fan of depth d.
winning=[sum(1<<i for i,h in enumerate(mh) if h>0)]
trees=[]
for d in range(1,3):
 win=0
 for i,h in enumerate(mh):
  if h<3 or 2*h+1>N:continue
  lo=sum(1<<j for j in by[2*h-1]);hi=sum(1<<j for j in by[2*h+1])
  if reach[i]&lo&winning[-1] and reach[i]&hi&winning[-1]:win|=1<<i
 winning.append(win)

def path(i,j):
 q=collections.deque([(i,())]);seen=set()
 while q:
  k,w=q.popleft()
  if k==j:return w
  if k in seen:continue
  seen.add(k)
  for l,t in enumerate(mt[k]):q.append((t,w+(l,)))
 raise RuntimeError(('path',i,j))

def witness(i,d,prefix=None):
 if prefix is None:prefix=rep[i]
 h=mh[i]
 if count(prefix)!=h:raise RuntimeError("physical prefix mismatch")
 res={'H':h,'prefix':''.join(map(str,prefix))}
 if not d:return res
 ch=[]
 for h1 in (2*h-1,2*h+1):
  eligible=reach[i]&winning[d-1]&sum(1<<j for j in by[h1])
  j=next(nodes(eligible));s=path(i,j)
  if count(prefix+s)!=h1:raise RuntimeError('physical child mismatch')
  ch.append({'suffix':''.join(map(str,s)),'H':h1,'state':j,'physical_child_word':''.join(map(str,prefix+s)),'continuation':witness(j,d-1,prefix+s) if d>1 else None})
 res['children']=ch;return res
rows=[]
for h in range(3,64,2):
 states=by[h];masks=set();example=None
 for i in states:
  mask=sum(1<<b for b,targ in enumerate((2*h-1,2*h+1)) if any(reach[i]>>j&1 for j in by[targ]))
  masks.add(mask)
  if mask==3 and example is None:example=witness(i,1)
 rows.append({'h':h,'source_continuation_classes':len(states),'available_branch_masks':sorted(masks),'common_parent_exists':bool(set(states)&set(nodes(winning[1]))),'two_level_common_parent_exists':bool(set(states)&set(nodes(winning[2]))) if h<=31 else None,'example':example})
 print(h,sorted(masks),rows[-1]['two_level_common_parent_exists'],flush=True)
separate=[]
for target in (45,47):
 for i in by[23]:
  js=[j for j in by[target] if reach[i]>>j&1]
  if js:
   j=js[0];suffix=path(i,j);prefix=rep[i]
   if count(prefix)!=23 or count(prefix+suffix)!=target:raise RuntimeError('separate branch')
   separate.append({'source_word':''.join(map(str,prefix)),'source_count':23,'suffix':''.join(map(str,suffix)),'target_word':''.join(map(str,prefix+suffix)),'target_count':target});break
with open('05-knowledge/results/session276_corridors_20261112/wildcard/READINESS.json','w') as f:json.dump({'separate_H23_branches':separate,'cap':N,'scope':'ALL source-word and suffix lengths; two-level only h<=31 fully inside cap','rows':rows,'two_level_examples':[witness(i,2) for i in nodes(winning[2]) if mh[i]>=23][:3]},f,indent=2)
