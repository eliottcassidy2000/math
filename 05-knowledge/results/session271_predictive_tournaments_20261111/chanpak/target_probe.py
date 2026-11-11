#!/usr/bin/env python3
"""Equal-budget heuristic target search: exact transfer scores, not certificates of absence."""
from predict import spine,HERE
import random,math,json
rng=random.Random(27191)
rows=[]
for d in [12,30,100]:
 N=10**d+1;m=math.ceil(math.log(8*N,3));budget=1200
 for method in ['iid','local_abs','residue_tiebreak']:
  best=None;cur=None;cw=None
  for t in range(budget):
   if method=='iid' or cw is None or t%120==0:
    w=[rng.randrange(4) for _ in range(m)]
   else:
    w=cw.copy()
    for _ in range(1 if t%5 else 3):w[rng.randrange(m)]=rng.randrange(4)
   h=spine.value(w);err=abs(h-N)
   # Integer relative-distance bucket + residue tie-break retains scale.
   if method=='residue_tiebreak':score=(err*1000//N,sum(min((h-N)%q,(N-h)%q)/q for q in [9,25,49]),err)
   else:score=(err,)
   if best is None or err<best[0]:best=(err,h,w)
   if method=='iid' or cur is None or t%120==0 or score<cur:cur=score;cw=w
  rows.append({'digits':d,'target':str(N),'m':m,'queries':budget,'method':method,'best_H':str(best[1]),'relative_error':best[0]/N,'exact_hit':best[0]==0,'word':best[2]})
(HERE/'target_probes.json').write_text(json.dumps(rows,indent=2)+'\n')
print(json.dumps([{k:v for k,v in x.items() if k not in ['word','best_H','target']} for x in rows],indent=2))
