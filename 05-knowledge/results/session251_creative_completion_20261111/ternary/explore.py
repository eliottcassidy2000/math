#!/usr/bin/env python3
"""Exact bounded support of a stated fixed-core bank; no unrestricted census."""
import argparse,json
from pathlib import Path
from compiler import bank,signature,generate,hp_subsets,construct
p=argparse.ArgumentParser();p.add_argument('--bound',type=int,default=65535);args=p.parse_args()
cores=[(0,),(2,4,1),(2,4,9,3),(2,4,9,19,7),(2,4,25,3,11)]
# Additional positive source grammar: two ordered triangles.
cores.append((2|56,4|56,1|56,16,32,8))
sigs={}
for rows in cores:
 for O in range(1<<len(rows)):
  s=signature(rows,O);sigs.setdefault(s,(rows,O))
covered={};counts={};single=set()
for sig,(r,O) in sorted(sigs.items()):
 for val,a,b in generate(sig,args.bound):
  covered.setdefault(val,(r,O,a,b,sig));
  if len(r)<=5:single.add(val)
  counts[val]=min(counts.get(val,100),sig[0])
misses=[h for h in range(1,args.bound+1,2) if h not in covered]
w={}
for h in (173,335,663,1383,2373,2457):
 if h in covered:
  r,O,a,b,s=covered[h];child=construct(r,O,a,b)
  if hp_subsets(child)[-1]!=h:raise AssertionError(h)
  w[str(h)]={'core':r,'O':O,'A':a,'B':b,'signature':s,'child':child}
report={'bound':args.bound,'source_H_at_most':9,'bank_description':'Five cores: K1, C3, the strong order-4 core, two strong order-5 H9 cores; plus the gapless ordered pair of triangles. Not all H<=9 sources.','core_count':len(cores),'single_core_first_misses':[h for h in range(1,args.bound+1,2) if h not in single][:100],'single_core_missing_count':sum(h not in single for h in range(1,args.bound+1,2)),'signature_count':len(sigs),'missing_count':len(misses),'first_misses':misses[:100],'targets_requiring_parent9_in_this_bank':[h for h in sorted(counts) if counts[h]==9],'parent_histogram':{str(d):list(counts.values()).count(d) for d in sorted(set(counts.values()))},'witnesses':w}
Path(__file__).with_name('explore-results.json').write_text(json.dumps(report,sort_keys=True,indent=2)+'\n');print(json.dumps(report,sort_keys=True,indent=2))
