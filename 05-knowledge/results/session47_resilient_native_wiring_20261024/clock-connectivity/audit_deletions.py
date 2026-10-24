#!/usr/bin/env python3
"""Independent direct27-combinations/Kruskal audit of deletion fibres."""
import argparse,json,math
from pathlib import Path
from itertools import combinations
from collections import defaultdict,Counter
HERE=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--output',type=Path);a=p.parse_args()
rows=[tuple(map(int,z.split()))for z in (HERE/'input.txt').read_text().splitlines()];src=[r[0]for r in rows];x=[r[2]for r in rows];y=[r[3]for r in rows];w=[r[4]for r in rows];frozen=json.loads((HERE/'deletion_probe.json').read_text());checks=0
def need(ok,msg):
 global checks
 checks+=1
 if not ok:raise RuntimeError(msg)
def work(m):return sum(w[i]for i in range(28)if m>>i&1)
report=[];bad=[]
for h in range(28):
 g=defaultdict(list);universe=[i for i in range(28)if i!=h];total=0
 for ii in combinations(universe,4):
  yy=0
  for i in ii:yy^=y[i]
  g[sum(x[i]for i in ii)%128,yy].append(sum(1<<i for i in ii));total+=1
 need(total==math.comb(27,4)==17550,'complete deleted inventory enumeration')
 hist=Counter();components22=Counter()
 for q,states in g.items():
  par=list(range(len(states)))
  def root(j):
   while par[j]!=j:par[j]=par[par[j]];j=par[j]
   return j
  edges=sorted((work(m^n),i,j)for i,m in enumerate(states)for j,n in enumerate(states[:i]));last=0;at22=None
  for c,i,j in edges:
   if c>22 and at22 is None:
    out=defaultdict(list)
    for t,m in enumerate(states):out[root(t)].append(m)
    at22=sorted(sorted(z)for z in out.values())
   u,v=root(i),root(j)
   if u!=v:par[u]=v;last=c
  if at22 is None:at22=[sorted(states)]
  hist[last]+=1;components22[len(at22)]+=1
  if len(at22)>1:bad.append({'deleted_source':src[h],'charge':list(q),'states':sorted(states),'components':at22,'bottleneck':last})
 orig=frozen['reports'][h]
 need(orig['deleted_source']==src[h] and orig['nonempty_fibres']==len(g),'independent nonempty-fibre counts')
 need({int(k):v for k,v in orig['component_counts'].items()}==dict(components22),'independent complete22 component histograms')
 report.append({'deleted_source':src[h],'states':total,'fibre_bottleneck_histogram':dict(sorted(hist.items())),'maximum_bottleneck':max(hist)})
expected=sorted((z['deleted_source'],tuple(z['charge']),tuple(sorted(z['states'])),tuple(sorted(tuple(sorted(c))for c in z['components'])))for z in frozen['hostiles'])
actual=sorted((z['deleted_source'],tuple(z['charge']),tuple(z['states']),tuple(tuple(c)for c in z['components']))for z in bad)
need(actual==expected,'all exact hostile states and components match independent algorithm')
need(len(bad)==8 and all(z['bottleneck']==23 for z in bad),'all hostile fibres sharp23')
res={'status':'PASS','checks':checks,'directly_enumerated_states':491400,'deleted_bank_count':28,'sharp_uniform_margin4_deleted_bank_clock':23,'bad_at22_count':8,'bad_deleted_sources':sorted({z['deleted_source']for z in bad}),'independent_kruskal_reports':report}
s=json.dumps(res,indent=2,sort_keys=True)+'\n'
if a.output:a.output.write_text(s)
else:print(s,end='')
