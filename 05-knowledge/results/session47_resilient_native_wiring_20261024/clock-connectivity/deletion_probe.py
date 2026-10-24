#!/usr/bin/env python3
"""All single-source-deletion margin-four fibres; exact cheap sheet gate test."""
import argparse,json
from pathlib import Path
from itertools import combinations
from collections import defaultdict,Counter
HERE=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,default=Path('session35-publication/clean-repo'));a=p.parse_args()
r=sorted(json.loads((a.repo/'05-knowledge/results/native_packet_decomposition_20261023/native/domains.json').read_text())['sources'],key=lambda z:z['source']);x=[z['code'][0]for z in r];y=[sum(b<<j for j,b in enumerate(z['code'][1:]))for z in r];w=[sum(a!=b for a,b in zip(z['old_before'],z['old_after']))for z in r]
three=sum(1<<i for i,v in enumerate(w)if v==3);four=sum(1<<i for i,v in enumerate(w)if v==4)
def work(m):return 2*m.bit_count()+(m&three).bit_count()+2*(m&four).bit_count()
g=defaultdict(list)
for ii in combinations(range(28),4):
 yy=0
 for i in ii:yy^=y[i]
 g[sum(x[i]for i in ii)%128,yy].append(sum(1<<i for i in ii))
reports=[];hostiles=[]
for h in range(28):
 hist=Counter();bad=0;nf=0
 for q,ss in g.items():
  states=[m for m in ss if not m>>h&1]
  if not states:continue
  nf+=1;left=set(states);cc=[]
  while left:
   seed=min(left);left.remove(seed);todo=[seed]
   for m in todo:
    hits=[n for n in left if work(m^n)<=22]
    left.difference_update(hits);todo.extend(hits)
   cc.append(sorted(todo))
  hist[len(cc)]+=1
  if len(cc)>1:
   bad+=1;hostiles.append({'deleted_source':r[h]['source'],'charge':q,'states':states,'components':cc,'least_cross_component_clock':min(work(a^b)for a in cc[0]for b in cc[1])})
 reports.append({'deleted_source':r[h]['source'],'nonempty_fibres':nf,'disconnected_at22':bad,'component_counts':dict(sorted(hist.items()))})
if len(hostiles)!=8 or sum(z['disconnected_at22']>0 for z in reports)!=7:raise RuntimeError('single-deletion counts changed')
if any(len(z['components'])!=2 or z['least_cross_component_clock']!=23 for z in hostiles):raise RuntimeError('sharp single-deletion23 crossing failed')
res={'status':'FINITE-EXACT single-deletion sheet gate','margin':4,'states_per_deleted_bank':17550,'reports':reports,'hostiles':hostiles}
(HERE/'deletion_probe.json').write_text(json.dumps(res,indent=2)+'\n');print('bad deleted banks',sum(z['disconnected_at22']>0 for z in reports),'bad fibres',len(hostiles));print(hostiles[:3])
