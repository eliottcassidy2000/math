#!/usr/bin/env python3
"""Independent complete-graph minimax, event-count and total-clock optimization."""
import argparse,hashlib,heapq,json
from collections import deque
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--packet',type=Path,required=True);p.add_argument('--routes',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
checks=0
def require(ok,why):
 global checks
 checks+=1
 if not ok:raise RuntimeError(why)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
r=json.loads(a.routes.read_text());p=json.loads(a.packet.read_text());slots={z['source']:z for z in p['source_options']}
require(r['status']=='PASS'and r['input_sha256']['packet']==sha(a.packet),'complete physical-state audit tied to unchanged packet')
reports=[]
for z in r['routes']:
 states=z['accepted_states'];start,end=z['advertised_path'][0],z['advertised_path'][-1];sources=z['sources'];costs=[sum(x!=y for x,y in zip(*[o['physical_path']for o in slots[s]['options']]))for s in sources];graph={x:[]for x in states}
 require(all(c==slots[s]['cost']for s,c in zip(sources,costs)),'weights from actual physical clock paths')
 for x in states:
  for y in states:
   if x==y:continue
   support=(x^y).bit_count();clock=sum(c for j,c in enumerate(costs)if(x^y)>>j&1);graph[x].append((y,support,clock))
 def minimax(index):
  d={start:0};h=[(0,start)]
  while h:
   cost,x=heapq.heappop(h)
   if d[x]!=cost:continue
   for edge in graph[x]:
    y=edge[0];new=max(cost,edge[index])
    if y not in d or new<d[y]:d[y]=new;heapq.heappush(h,(new,y))
  return d[end]
 def optimize(index,budget):
  levels={start:0};q=deque([start])
  while q:
   x=q.popleft()
   for edge in graph[x]:
    y=edge[0]
    if edge[index]<=budget and y not in levels:levels[y]=levels[x]+1;q.append(y)
  d={start:(0,0)};parent={start:None};h=[(0,0,start)]
  while h:
   work,count,x=heapq.heappop(h)
   if d[x]!=(work,count):continue
   for edge in graph[x]:
    y,_,clock=edge
    if edge[index]>budget:continue
    new=(work+clock,count+1)
    if y not in d or new<d[y]:d[y]=new;parent[y]=x;heapq.heappush(h,(*new,y))
  path=[end]
  while path[-1]!=start:path.append(parent[path[-1]])
  path.reverse();weights=[]
  for x,y in zip(path,path[1:]):weights.append(next((s,c)for t,s,c in graph[x]if t==y))
  require(sum(c for s,c in weights)==d[end][0]and len(path)-1==d[end][1],'minimum-work path exact accounting')
  require(d[end][1]==levels[end],'minimum-work path also attains minimum event count in this menu')
  return {'threshold':budget,'minimum_events':levels[end],'minimum_total_clock_work':d[end][0],'minimum_work_path':path,'step_slot_supports':[s for s,c in weights],'step_clock_costs':[c for s,c in weights]}
 bslot=minimax(1);bclock=minimax(2)
 require(bslot==z['exact_minimax_slot_budget']and bclock==z['exact_minimax_ordinary_clock_budget'],'independent weighted minimax Dijkstra reproduces exhaustive threshold scan')
 result={'name':z['name'],'at_optimal_slot_budget':optimize(1,bslot),'at_optimal_clock_budget':optimize(2,bclock)}
 reports.append(result)
expected={'buffer_3015':(4,50),'buffer_1394':(3,40),'buffer_1394_2905':(4,48),'eighteen_buffer_1394':(4,55),'eighteen37_buffer_2071':(4,49)}
for z in reports:
 if z['name']in expected:
  o=z['at_optimal_slot_budget'];require((o['minimum_events'],o['minimum_total_clock_work'])==expected[z['name']],'independently checked finite optimization summary')
cost37=next(z for z in reports if z['name']=='eighteen37_buffer_2071')['at_optimal_clock_budget'];require((cost37['threshold'],cost37['minimum_events'],cost37['minimum_total_clock_work'])==(12,5,57),'exact cost37 Pareto boundary: lower clock bottleneck needs five events and work57')
out={'status':'PASS','checks':checks,'input_sha256':{'packet':sha(a.packet),'audited_complete_routes':sha(a.routes)},'routes':reports,'scope':'Exact finite optimization in each separately authenticated complete pinned legal fibre. Full pairwise graph, actual owner-clock distances. Does not optimize over other buffer sets or undeclared sources; strict temporal nonreturn is not imposed on these reversible routes.'}
a.out.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps({'status':'PASS','checks':checks,'routes':reports},sort_keys=True))
