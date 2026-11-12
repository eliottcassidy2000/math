"""All-length subsequence obstruction and a physical old-arc repair at H13."""
from capped import build,minimize,count
from fibres import graph,ham
from collections import deque
import json
vs,ts,ws,hs=build(27)
tr,out,part=minimize(ts,[int(h==27) for h in hs],0)
cores=['120','212','230','302','310','321'];cert=[]
for text in cores:
 core=tuple(map(int,text));q=deque([(part[0],0)]);seen=set();accepted=[]
 while q:
  s,j=q.popleft()
  if (s,j) in seen:continue
  seen.add((s,j))
  if out[s] and j==len(core):accepted.append((s,j))
  for l,t in enumerate(tr[s]):q.append((t,j+int(j<len(core) and l==core[j])))
 if accepted:raise RuntimeError(('unexpected insertion',text,accepted))
 cert.append({'H13core':text,'reachable_product_states':len(seen),'accepting_product_states':0})
source=graph((1,2,0));lower=graph((1,1,2,0));padded=graph((0,1,2,0));mid=padded.copy();mid[1]^=1<<2;mid[2]^=1<<1;upper=mid.copy();upper[0]^=1<<2;upper[2]^=1<<0
checks=[ham(g) for g in (source,lower,padded,mid,upper)]
if checks!=[13,25,13,15,27]:raise RuntimeError(('repair',checks))
# new upper fixed spine order is old labels(2,0,1,3), same apexes4,5
perm=[2,0,1,3,4,5]
renamed=[sum(1<<j for j,v in enumerate(perm) if upper[u]>>v&1) for u in perm]
if renamed!=graph((2,0,1,0)):raise RuntimeError('reframe')
hostile=graph((0,1,1,2,0))
for i,j in [(2,3),(0,3)]:hostile[i]^=1<<j;hostile[j]^=1<<i
if ham(hostile)!=71 or count((1,1,2,0))!=25:raise RuntimeError('hostile repeat')
out={'hostile_repeat':{'source_word':'1120','source_count':25,'neutral_prefix_word':'01120','flips':[[2,3],[0,3]],'result_rows':hostile,'result_count':71,'desired_upper_binary_count':51},'target':27,'scope':'all source word lengths and all interior insertions retaining the apex pair, old arcs and spine subsequence','core_certificates':cert,'repair':{'source_word':'120','source_rows':source,'lower_word':'1120','lower_rows':lower,'upper_word_after_reframe':'2010','padded_rows':padded,'old_arc_flip':[1,2],'intermediate_rows':mid,'new_arc_flip':[0,2],'upper_rows':upper,'counts':checks,'upper_spine_reframe':[2,0,1,3]}}
with open('05-knowledge/results/session276_corridors_20261112/wildcard/INSERTION.json','w') as f:json.dump(out,f,indent=2)
print(json.dumps(out,indent=2))
