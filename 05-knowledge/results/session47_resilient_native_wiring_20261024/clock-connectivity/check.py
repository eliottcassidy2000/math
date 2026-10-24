#!/usr/bin/env python3
"""Optimization-safe independent theorem/receipt and route audit."""
import argparse,hashlib,json,math
from pathlib import Path
from collections import Counter,deque
from itertools import combinations
HERE=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,default=Path('session35-publication/clean-repo'));p.add_argument('--output',type=Path);args=p.parse_args()
RESULTS=args.repo/'05-knowledge/results';domain=RESULTS/'native_packet_decomposition_20261023/native/domains.json'
rows=sorted(json.loads(domain.read_text())['sources'],key=lambda r:r['source']);n=len(rows);checks=0
def need(ok,msg):
 global checks
 checks+=1
 if not ok:raise RuntimeError(msg)
need(n==28,'full28 source supply')
src=[r['source']for r in rows];X=[r['X']for r in rows];Y=[sum(b<<j for j,b in enumerate(r['code'][1:]))for r in rows];W=[sum(a!=b for a,b in zip(r['old_before'],r['old_after']))for r in rows]
need(len(set((r['code'][0],y)for r,y in zip(rows,Y)))==28,'all source charge labels distinct')
expected='\n'.join(f'{s} {x} {r["code"][0]} {y} {w}'for s,x,r,y,w in zip(src,X,rows,Y,W))+'\n'
need((HERE/'input.txt').read_text()==expected,'enumerators use exact authenticated native labels and costs')
need(Counter(W)=={2:20,3:5,4:3},'native clock distribution')
def ids(m):return [i for i in range(n)if m>>i&1]
def info(m):
 ii=ids(m);yy=0
 for i in ii:yy^=Y[i]
 return len(ii),sum(X[i]for i in ii),yy,sum(W[i]for i in ii)
def mask(ss):return sum(1<<src.index(s)for s in ss)
def read(name):return [json.loads(s)for s in (HERE/name).read_text().splitlines()]
producer=read('budget22-support8.jsonl');auditor=read('independent-budget22-support8.jsonl');plain=read('budget22.jsonl');depth=read('independent-depth-budget22-support8.jsonl');need(len(depth)==15,'complete depth replay');need(len(producer)==len(auditor)==len(plain)==15,'all fifteen half-margin runs completed')
summaries=[]
for k,(a,b,c)in enumerate(zip(producer,auditor,plain)):
 need(a['margin']==b['margin']==c['margin']==k,'margin coverage')
 need(a['states']==b['states']==c['states']==math.comb(28,k),'exact exhaustive subset counts')
 for key in ('budget','fibres','disconnected_fibres','component_counts'):need(a[key]==b[key]==c[key],'two independent graph algorithms and relaxed run agree: '+key)
 need(a['support']==b['support']==8 and a['budget']==22 and a['disconnected_fibres']==0,'joint8slot22clock connectivity')
 need(a['component_counts']=={'1':a['fibres']},'every nonempty fibre connected')
 need(a['fibres']==([1,28,351,1646][k]if k<4 else 2048),'all actual nonempty charge counts')
 need(all(depth[k][z]==a[z]for z in ('margin','budget','support','states','fibres','disconnected_fibres','component_counts')),'depth replay matches complete components')
 need(depth[k]['maximum_root_depth']==[0,0,1,1,2,2,2,3,3,3,3,3,3,4,4][k],'exact recorded maximum root depths')
 summaries.append({kk:a[kk]for kk in ('margin','states','fibres','component_counts')})
need(sum(z['states']for z in summaries)==154276028,'lower-half exhaustive state total')
totalY=0
for y in Y:totalY^=y
need(sum(X)==3733 and sum(X)%128==21 and totalY==15,'complement affine charge centre')
# Independent complete 4-subset reconstruction of the sharp lower witness.
A={1532,1691,1923,1974};B={1841,2176,2837,3360};a,b=mask(A),mask(B);found=[]
for ii in combinations(range(n),4):
 yy=0
 for i in ii:yy^=Y[i]
 if sum(X[i]for i in ii)%128==123 and yy==11:found.append(sum(1<<i for i in ii))
need(set(found)=={a,b},'entire full28 margin4 charge123,11 fibre has exactly two states')
need(a&b==0 and info(a)[1:3]==(507,11) and info(b)[1:3]==(635,11),'lower witness exact integer sums and parity')
need(info(a^b)[0]==8 and info(a^b)[3]==22,'joint sharp lower distances')
cycle_details=[]
for i in ids(a^b):
 before,after=rows[i]['old_before'],rows[i]['old_after'];f=dict(zip(before,after));seen=set();cycles=[]
 for v in range(8):
  if v in seen:continue
  cy=[];u=v
  while u not in seen:seen.add(u);cy.append(u);u=f[u]
  if len(cy)>1:cycles.append(cy)
 need(len(cycles)==1 and len(cycles[0])==W[i],'literal single-cycle witness costs')
 cycle_details.append({'source':src[i],'cycle':cycles[0],'clock_cost':W[i]})
need(Counter(z['clock_cost']for z in cycle_details)=={2:4,3:2,4:2},'2/3/4 cycle mechanism')
# Direct combinations, independent of the route producer's 14+14 MITM.
route=json.loads((HERE/'route18_37.json').read_text());full=[];exact=[]
for ii in combinations(range(n),9):
 xx=sum(X[i]for i in ii)
 if xx%128!=9:continue
 yy=0
 for i in ii:yy^=Y[i]
 if yy!=7:continue
 m=sum(1<<i for i in ii);full.append(m)
 if xx==1161:exact.append(m)
need(len(full)==3337 and len(exact)==1355,'independent complete full28 circuit-context fibre counts')
start=mask({1477,1549,1691,1921,1974,2837,2905,3000,3369});end=mask({1394,1425,1457,1480,2176,2185,2581,2835,2945});core=start|end;outside=((1<<28)-1)^core
restricted=[m for m in full if not m&outside]
need(set(restricted)=={start,end},'complete restricted18-slot fibre remains primitive')
need(info(start^end)[3]==37,'restricted18slot exact clock37')
three=sum(1<<i for i,v in enumerate(W)if v==3);four=sum(1<<i for i,v in enumerate(W)if v==4)
def work(m):return 2*m.bit_count()+(m&three).bit_count()+2*(m&four).bit_count()
def component(states,budget):
 left=set(states);left.remove(start);todo=deque([start]);out=[start]
 while todo:
  m=todo.popleft();found=[z for z in left if work(m^z)<=budget]
  for z in found:left.remove(z);todo.append(z);out.append(z)
 return sorted(out)
route_summary={}
for key,states,minimum,total in [('dyadic',full,12,57),('exact_X',exact,13,49)]:
 rr=route[key]['route'];path=[z['mask']for z in rr['states']]
 need(path[0]==start and path[-1]==end and all(m in states for m in path),'verified endpoint route all arrivals in full authenticated fibre')
 need(all(not m&outside for m in (path[0],path[-1])),'all ten context slots restored')
 for m,z in zip(path,rr['states']):
  need(z['V_sources']==[src[i]for i in ids(m)] and z['integer_X_sum']==info(m)[1],'exact arrival sources/integer source sums')
  need(z['U_integer_current_change']==2*(info(m)[1]-1161),'actual target-U integer current at fixed margin')
 for (m,z),ev in zip(zip(path,path[1:]),rr['events']):
  need(ev['changed_sources']==[src[i]for i in ids(m^z)] and ev['ordinary_owner_clock_work']==work(m^z) and ev['changed_slots']==(m^z).bit_count(),'actual event source and owner-clock work')
 need(max(work(m^z)for m,z in zip(path,path[1:]))==minimum and sum(work(m^z)for m,z in zip(path,path[1:]))==total,'reported route peak and total')
 cut=component(states,minimum-1);need(end not in cut,'complete lower-budget component separates endpoints')
 route_summary[key]={'states':len(states),'sharp_endpoint_clock_bottleneck':minimum,'verified_route_events':len(path)-1,'route_total_clock_work':total,'lower_budget':minimum-1,'lower_component':cut,'lower_component_size':len(cut)}
inputs=[domain,RESULTS/'native_packet_decomposition_20261023/native/check.py',RESULTS/'native_packet_decomposition_20261023/native/check_receipt.json',RESULTS/'metagraph_order8_20261023/merged_nodes.tsv',RESULTS/'metagraph_order8_20261023/owner_map.u16.gz',RESULTS/'metagraph_order8_electrical_20261023/receipt.json']
result={'status':'PASS','checks':checks,'scope':'FINITE-EXACT: complete28-source fixed-margin full-paired-dyadic fibres, joint support8 clock22; no all-tournament-order theorem','directly_enumerated_states':154276028,'complement_covered_states':268435456,'all_margins_nonempty_fibres':47060,'complement_charge':[21,15],'dense_root_depths':[z['maximum_root_depth']for z in depth],'dense_route_event_upper_bound':8,'dense_route_clock_upper_bound':176,'half_margin_census':summaries,'witness22':{'charge':[123,11],'states':[sorted(A),sorted(B)],'masks':[a,b],'integer_X_sums':[507,635],'target_U_integer_current_difference':256,'cycles':cycle_details},'route18_37_audit':route_summary,'input_sha256':{q.relative_to(args.repo).as_posix():hashlib.sha256(q.read_bytes()).hexdigest()for q in inputs}}
out=json.dumps(result,indent=2,sort_keys=True)+'\n'
if args.output:args.output.write_text(out)
else:print(out,end='')
