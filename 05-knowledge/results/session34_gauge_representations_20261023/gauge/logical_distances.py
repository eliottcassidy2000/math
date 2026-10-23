from itertools import combinations
from collections import Counter
import json
E=list(combinations(range(7),2));A=[4,1,3,7,6,5,2];L=[{i,(i+1)%7,(i+3)%7} for i in range(7)]
line=[sum(1<<j for j,e in enumerate(E) if set(e)<=l) for l in L]
cuts=[sum((((g>>a)^(g>>b))&1)<<j for j,(a,b) in enumerate(E)) for g in range(64)]
rows=[];kernel=set();counter=Counter();dist_hist=Counter()
for m in range(128):
 s=0;x=0
 for i in range(7):
  if m>>i&1:s^=A[i];x^=line[i]
 if s:continue
 coset={x^b for b in cuts};kernel|=coset;wh=Counter(y.bit_count() for y in coset);counter.update(wh)
 d=min(wh);dist_hist[d]+=1
 rows.append({'line_mask':m,'line_weight':m.bit_count(),'minimal_edge_weight':d,'minimizer_count':wh[d],'edge_weight_enumerator':dict(sorted(wh.items()))})
if len(kernel)!=1024 or sum(counter.values())!=1024 or len(rows)!=16:raise RuntimeError('count failure')
print(json.dumps({'status':'PASS','kernel_weight_enumerator':dict(sorted(counter.items())),'logical_distance_histogram':dict(sorted(dist_hist.items())),'logical_cosets':rows},indent=2))
