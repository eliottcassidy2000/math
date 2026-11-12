"""New genuine clock births and a complete narrowly scoped continuation test."""
from pathlib import Path
from itertools import product
import json
from explore import val
from tensor import kernel
from cones import twinquot,maps
from normalforms import need
HERE=Path(__file__).resolve().parent
CASES=[{'rows':[158,156,8,16,36,335,159,60,223],'O':65,'ports':[1,2,8]}, {'rows':[670,668,8,528,36,847,671,60,735,148],'O':65,'ports':[1,2,8,9]}]

def strong(r):
 n=len(r)
 for root in range(n):
  reached={root};todo=[root]
  while todo:
   v=todo.pop()
   for w in range(n):
    if r[v]>>w&1 and w not in reached:reached.add(w);todo.append(w)
  if len(reached)!=n:return False
 return True

def run():
 cert=[]
 for item in CASES:
  r=item['rows'];O=item['O'];P=item['ports'];n=len(r)
  need(strong(r),'strong birth')
  q,o,cc=twinquot(r,O);need(len(cc)==n,'genuine twin-reduced core')
  need(len(maps(q,o)[1])==1,'rigid marked quotient')
  K=kernel(r,O,P);need(K=={(1,)*len(P):1},'full tensor')
  tests=[]
  for sizes in [tuple([1]*len(P)),tuple([2]*len(P))]+[tuple(3 if i==j else 1 for i in range(len(P))) for j in range(len(P))]:
   d=val(r,O,dict(zip(P,sizes)));need(d==1,'literal length control');tests.append({'sizes':sizes,'delta':d})
  cert.append({**item,'n':n,'kernel':[[list(c),v] for c,v in K.items()],'literal_controls':tests,'series':f'z^{n}/(1-z)^{len(P)}'})
 # Exhaustive next transitive clock, with the six old exterior vertices fixed.
 r=CASES[-1]['rows'];O=65;P=CASES[-1]['ports'];n=len(r);ext=[v for v in range(n) if v not in P]
 order=sorted(P,key=lambda v:-sum(bool(r[v]>>u&1) for u in P));counts={'templates':0,'strong':0,'twin_reduced':0,'unit':0,'all_single_double_screens':0};survivors=[]
 for gap in range(len(P)+1):
  for sig in range(1<<len(ext)):
   counts['templates']+=1
   cut=sum(1<<v for v in order[gap:])+sum(1<<v for i,v in enumerate(ext) if sig>>i&1)
   if cut in (0,(1<<n)-1):continue
   counts['strong']+=1
   rr=[r[v]|((1-(cut>>v&1))<<n) for v in range(n)]+[cut]
   q,o,cc=twinquot(rr,O)
   if len(cc)!=n+1:continue
   counts['twin_reduced']+=1
   if val(rr,O)!=1:continue
   counts['unit']+=1
   pp=P+[n]
   if any(val(rr,O,{v:2})!=1 for v in pp):continue
   counts['all_single_double_screens']+=1;survivors.append({'rows':rr,'O':O,'ports':pp})
 need(not survivors,'unexpected continuation; revise theorem')
 out={'scope':'Two exact jointly stable strong twin-reduced charts, plus all320 transitive-chain one-new-clock templates from the displayed ten-core with six exterior vertices and marks fixed. Does not exclude other operations or other cores.','births':cert,'continuation_counts':counts,'five_clock_survivors':survivors,'first_three_clock_order':{'order':9,'lower_bound_source':'inherited complete strong marked stable-port bank through8; at most two canonical clocks'},'four_clock_order10_is_not_claimed_minimal':True}
 print(json.dumps(out,sort_keys=True,indent=2));return out
if __name__=='__main__':
 out=run();(HERE/('clocks.optimized.json' if not __debug__ else 'clocks.json')).write_text(json.dumps(out,sort_keys=True,indent=2)+'\n')
