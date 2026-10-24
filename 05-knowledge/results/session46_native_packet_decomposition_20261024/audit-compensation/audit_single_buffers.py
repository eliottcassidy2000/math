#!/usr/bin/env python3
"""Independent exhaustive boundary: every single supplied buffer and both starting bits."""
import argparse,json,hashlib
from collections import deque
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,required=True);p.add_argument('--foundation',type=Path,required=True);p.add_argument('--packet',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
checks=0
def require(ok,why):
 global checks
 checks+=1
 if not ok:raise RuntimeError(why)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
Epath=a.repo/'05-knowledge/results/metagraph_order8_electrical_20261023/receipt.json';E=json.loads(Epath.read_text())['order8'];Y=E['radical_lifts_mod4'];X=E['higher_layers'][-1]['lifts'][0];chars=[[64*y for y in q]for q in Y[:4]]+[X];U,V=1965,2234;N=3528
F=json.loads(a.foundation.read_text());packet=json.loads(a.packet.read_text());base=packet['base_source_order'];sources={r['source']:r for r in F['sources']}
require(F['status']=='PASS'and F['input_sha256']['electrical']==sha(Epath),'independent foundation uses same original harmonic vectors')
require(len(base)==16 and all(8-len(sources[s]['fixed_clocks'])==2 for s in base),'sixteen original sources all cost two')
buffers=sorted(s for s,r in sources.items()if s not in base);cost2=[s for s in buffers if 8-len(sources[s]['fixed_clocks'])==2]
require(cost2==[1394,1921,2905,3369],'all available external cost-two sources, not a chosen list')
def current(ss,mask,chi):
 out=[0]*N
 for j,s in enumerate(ss):
  t=V if mask>>j&1 else U;z=2*(chi[s]-chi[t]);out[s]+=z;out[t]-=z
 return out
def reach(states,start,end,budget,costs=None):
 parents={start:None};q=deque([start])
 while q:
  old=q.popleft()
  for new in states:
   weight=(old^new).bit_count()if costs is None else sum(c for j,c in enumerate(costs)if(old^new)>>j&1)
   if new not in parents and weight<=budget:parents[new]=old;q.append(new)
 if end not in parents:return sorted(parents),None
 path=[end]
 while path[-1]!=start:path.append(parents[path[-1]])
 return sorted(parents),path[::-1]
rows=[]
for b in buffers:
 ss=base+[b];contrib=[tuple(-2*(chi[s]-chi[V])%256 for chi in chars)for s in ss];sums=[(0,)*5]*(1<<17)
 for mask in range(1,1<<17):
  low=mask&-mask;j=low.bit_length()-1;sums[mask]=tuple((x+y)%256 for x,y in zip(sums[mask^low],contrib[j]))
 for initial in(0,1):
  start=65280+(initial<<16);end=255+(initial<<16);baseline=[current(ss,start,chi)for chi in chars];target=tuple(row[V]%256 for row in baseline)
  legal=[mask for mask in range(1<<17)if mask.bit_count()==start.bit_count()and sums[mask]==target]
  require(start in legal and end in legal,'both original endpoints present')
  for mask in legal:
   require(all(all((z-old)%256==0 for z,old in zip(current(ss,mask,chi),baseline[j]))for j,chi in enumerate(chars)),'every accepted state full nodewise physical-row current')
   ports=set()
   for j,s in enumerate(ss):
    t=V if mask>>j&1 else U;choice=next(z for z in sources[s]['complete_pinned_domain']if z['target']==t);T=set(choice['closed_ports'])
    require(len(T)==4 and not ports.intersection(T),'every accepted state genuine non-reused resource');ports.update(T)
  six,_=reach(legal,start,end,6)
  for threshold in range(2,18,2):
   _,path=reach(legal,start,end,threshold)
   if path is not None:slot_min=threshold;slot_path=path;break
  costs=[8-len(sources[s]['fixed_clocks'])for s in ss]
  for threshold in range(1,1+sum(costs)):
   _,path=reach(legal,start,end,threshold,costs)
   if path is not None:clock_min=threshold;clock_path=path;break
  require(b not in cost2 or(slot_min==8 and clock_min==16),'each cost-two buffer has exact eight-slot and sixteen-clock bottleneck')
  rows.append({'buffer':b,'buffer_clock_cost':costs[-1],'initial':initial,'legal_states':legal,'six_move_start_component':six,'slot_optimal_path':slot_path,'clock_optimal_path':clock_path,'exact_slot_bottleneck':slot_min,'exact_ordinary_clock_bottleneck':clock_min})
require(len(rows)==24 and min(r['exact_ordinary_clock_bottleneck']for r in rows)==13,'all twelve single buffers and both initial bits: sharp global clock minimum thirteen')
report={'status':'PASS','checks':checks,'input_sha256':{'electrical':sha(Epath),'foundation':sha(a.foundation),'packet':sha(a.packet)},'cases':rows,'scope':'All twelve external sources from the authenticated28 supply, both initial bits; complete fixed-clock U/V menus and all-dyadic physical current via independently certified bank. No assertion for additional source supply.'}
a.out.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');print(json.dumps({'status':'PASS','checks':checks,'cases':[[r['buffer'],r['initial'],len(r['legal_states']),len(r['six_move_start_component']),r['exact_slot_bottleneck']]for r in rows]},sort_keys=True))
