"""Independent exact budget-feasibility audit of all endpoint observer fibres."""
from collections import Counter,defaultdict
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import argparse,json
HERE=Path(__file__).resolve().parent
CHECKS=0

def need(ok,msg):
 global CHECKS
 CHECKS+=1
 if not ok:raise RuntimeError(msg)

def costs(rows,key):
 paths=[tuple(r['path']) for r in rows];values=[tuple(r['path']) if key=='path' else r[key] for r in rows]
 @lru_cache(None)
 def feasible(ids,b):
  if len({values[i] for i in ids})==1:return True
  if not b:return False
  for j in range(1,7):
   bins=defaultdict(list)
   for i in ids:bins[paths[i][j]].append(i)
   if len(bins)>1 and all(feasible(tuple(v),b-1) for v in bins.values()):return True
  return False
 root=tuple(range(len(rows)));adaptive=next(b for b in range(7) if feasible(root,b))
 for b in range(7):
  found=False
  for clocks in combinations(range(1,7),b):
   bins=defaultdict(set)
   for p,v in zip(paths,values):bins[tuple(p[j] for j in clocks)].add(v)
   if all(len(v)==1 for v in bins.values()):found=True;break
  if found:break
 return dict(adaptive=adaptive,nonadaptive=b)

def main():
 p=argparse.ArgumentParser();p.add_argument('--packet',type=Path,default=HERE.parent/'session44-actual-guards');p.add_argument('--write',action='store_true');args=p.parse_args();P=args.packet
 files=['ENDPOINT_OBSERVERS.md','endpoint_observers.py','endpoint.normal.json','endpoint.optimized.json','ENDPOINT_MANIFEST.json','path_records.json']
 hashes={f:sha256((P/f).read_bytes()).hexdigest() for f in files}
 need(hashes==json.loads((HERE/'ENDPOINT_SOURCE_HASHES.json').read_text()),'frozen endpoint snapshot')
 need((P/'endpoint.normal.json').read_bytes()==(P/'endpoint.optimized.json').read_bytes(),'producer endpoint normal optimized equality')
 R=json.loads((P/'endpoint.normal.json').read_text());data=json.loads((P/'path_records.json').read_text());out=[];total=0
 for sr in R['sources']:
  source=sr['source'];records=data[str(source)];bags=defaultdict(list)
  for r in records:bags[tuple((r['path'][0],r['path'][-1]))].append(r)
  need(len(bags)==len(sr['fibres']),'complete endpoint partition')
  hist={key:Counter() for key in ('phase','target','path')};fixed={key:Counter() for key in hist};improve={key:[] for key in hist}
  for row in sr['fibres']:
   pins=tuple(row['endpoints']);rows=bags[pins];need(len(rows)==row['actual_paths'],'exact endpoint mass')
   total+=1
   for key in hist:
    c=costs(rows,key);need(c==row['costs'][key],'endpoint exact adaptive lower and upper plus nonadaptive menu')
    hist[key][c['adaptive']]+=1;fixed[key][c['nonadaptive']]+=1
    if c['adaptive']<c['nonadaptive']:improve[key].append(dict(endpoints=pins,paths=len(rows),costs=c))
  need(all(max(h)==3 for h in hist.values()),'three-query sharp maximum separately for all objectives')
  need({key:{str(k):v for k,v in h.items()} for key,h in hist.items()}==sr['extra_adaptive_depth_distributions'],'complete adaptive distributions')
  need({key:{str(k):v for k,v in h.items()} for key,h in fixed.items()}==sr['extra_nonadaptive_depth_distributions'],'complete nonadaptive distributions')
  out.append(dict(source=source,endpoint_fibres=len(bags),adaptive_histograms=hist,nonadaptive_histograms=fixed,adaptive_improvements=improve))
 need(total==210,'all210 endpoint fibres')
 text=json.dumps(dict(status='PASS independent endpoint lower and upper certificates',checks=CHECKS,source_sha256=hashes,sources=out,total_endpoint_fibres=total),indent=2,sort_keys=True)+'\n'
 if args.write:(HERE/'endpoint_receipt.json').write_text(text)
 else:need((HERE/'endpoint_receipt.json').read_text()==text,'frozen endpoint audit receipt')
 print(json.dumps(dict(status='PASS',checks=CHECKS,endpoint_fibres=total,receipt_sha256=sha256(text.encode()).hexdigest()),sort_keys=True))

if __name__=='__main__':main()
