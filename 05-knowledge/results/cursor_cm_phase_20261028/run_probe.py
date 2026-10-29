"""Exact supplied-run contexts: C3 quotient, two auxiliary singleton-run children."""
from itertools import permutations,product
from collections import defaultdict
from pathlib import Path
import json
import cursor_probe as c
HERE=Path(__file__).resolve().parent
o=c.o
def all_tournaments(n):
 for mask in range(1<<(n*(n-1)//2)):
  t=[0]*n;k=0
  for i in range(n):
   for j in range(i+1,n):
    a,b=(i,j) if mask>>k&1 else (j,i);t[a]|=1<<b;k+=1
  yield tuple(t)
def representatives(n):
 buckets=defaultdict(list)
 for t in all_tournaments(n):
  bucket=buckets[o.colours(t)[1]]
  if not any(o.isomorphisms(t,u) for u in bucket):bucket.append(t)
 return [t for bucket in buckets.values() for t in bucket]
def covers(t):
 n=len(t);out=set()
 for p in permutations(range(n)):
  for cuts in range(1<<(n-1)):
   runs=[];run=[p[0]];okay=True
   for j in range(n-1):
    if cuts>>j&1:runs.append(tuple(run));run=[p[j+1]]
    elif c.c.edge(t,p[j],p[j+1]):run.append(p[j+1])
    else:okay=False;break
   if okay:out.add(tuple(sorted(runs+[tuple(run)])))
 return sorted(out)
def act_cover(cover,g):return tuple(tuple(g[x] for x in run) for run in cover)
def trans(n):return tuple(sum(1<<j for j in range(i+1,n)) for i in range(n))
def context(a,run_order,b=None,d=None):
 k=len(run_order);m=k-1
 if k==1:
  path=run_order[0];return c.c.tau(c.c.relabel(a,path))
 b=trans(m) if b is None else b;d=trans(m) if d is None else d
 children=(a,b,d);owners=[(i,j) for i,t in enumerate(children) for j in range(len(t))]
 source=tuple(sum((c.c.edge(children[i],x,y) if i==j else ((j-i)%3==1))<<v
                  for v,(j,y) in enumerate(owners) if (i,x)!=(j,y)) for i,x in owners)
 path=[]
 for j,run in enumerate(run_order):
  path.extend(run)
  if j<m:path.extend((len(a)+j,len(a)+m+j))
 if not all(c.c.edge(source,x,y) for x,y in zip(path,path[1:])):raise RuntimeError('illegal path')
 return c.c.tau(c.c.relabel(source,path))
def output_partition(outputs):
 buckets=defaultdict(list);classes=[]
 for i,t in enumerate(outputs):
  key=o.colours(t)[1];slot=None
  for j in buckets[key]:
   if o.isomorphisms(t,outputs[classes[j][0]]):slot=j;break
  if slot is None:slot=len(classes);classes.append([]);buckets[key].append(slot)
  classes[slot].append(i)
 return classes
def main():
 rows=[];bad=[];controls=[];total=0
 for n in range(2,6):
  sources=representatives(n);tested=0;distinct=0;failed=0
  for a in sources:
   aut=o.isomorphisms(a,a,True);seen=set()
   for F in covers(a):
    # One representative of each forest under the whole source group.
    key=min(tuple(sorted(act_cover(F,g))) for g in aut)
    if key in seen:continue
    seen.add(key)
    G=[g for g in aut if tuple(sorted(act_cover(F,g)))==F]
    orders=sorted({min(act_cover(P,g) for g in G) for P in permutations(F)})
    if len(orders)<2:continue
    outs=[context(a,P) for P in orders];parts=output_partition(outs)
    tested+=1;distinct+=len(orders);total+=1
    if any(len(b)>1 for b in parts):
     failed+=1
     for block in parts:
      if len(block)>1 and len(bad)<30:
       bad.append({'source':a,'forest':F,'orders':[orders[x] for x in block],
                   'group':len(G),'output':outs[block[0]],'orbit_count':len(orders),'observed_cells':len(parts)})
    if len(controls)<10:controls.append({'source':a,'forest':F,'group':len(G),'orbit_count':len(orders),'observed_cells':len(parts)})
  row={'n':n,'source_classes':len(sources),'forests_tested':tested,'order_orbits':distinct,'one_context_failures':failed}
  rows.append(row);print(json.dumps(row),flush=True)
 result={'status':'FINITE-EXACT','menu':'C3 walk A,B,C repeated, transitive auxiliary singleton runs','rows':rows,'first_failures':bad,'controls':controls}
 (HERE/'run_probe.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
