"""Count all disjoint balanced squarefree moves of half-size at most four."""
import argparse,json
from pathlib import Path
from itertools import combinations
from collections import defaultdict,Counter
P=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--repo',type=Path,default=Path('session35-publication/clean-repo'));args=ap.parse_args()
r=sorted(json.loads((args.repo/'05-knowledge/results/native_packet_decomposition_20261023/native/domains.json').read_text())['sources'],key=lambda t:t['source'])
x=[t['code'][0]for t in r];y=[sum(v<<j for j,v in enumerate(t['code'][1:]))for t in r];w=[sum(a!=b for a,b in zip(t['old_before'],t['old_after']))for t in r]
def sub(ids):
 yy=0
 for i in ids:yy^=y[i]
 return sum(x[i]for i in ids)%128,yy
report=[]
for k in (2,3,4):
 g=defaultdict(list)
 for ss in combinations(range(28),k):g[sub(ss)].append(ss)
 hist=Counter();prim=Counter();allowed=Counter();pa=Counter()
 for states in g.values():
  for a,b in combinations(states,2):
   if set(a)&set(b):continue
   cost=sum(w[i]for i in a+b);hist[cost]+=1
   isp=not any(sub(c)==sub(d)for c in combinations(a,2)for d in combinations(b,2))if k==4 else True
   if isp:prim[cost]+=1
   if cost<=22:
    allowed[cost]+=1
    if isp:pa[cost]+=1
 report.append({'changed_slots':2*k,'all_moves':sum(hist.values()),'all_weight_histogram':dict(sorted(hist.items())),'primitive_moves':sum(prim.values()),'primitive_weight_histogram':dict(sorted(prim.items())),'allowed_moves':sum(allowed.values()),'allowed_primitive_moves':sum(pa.values()),'allowed_primitive_weight_histogram':dict(sorted(pa.items()))})
(P/'move_bank.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
