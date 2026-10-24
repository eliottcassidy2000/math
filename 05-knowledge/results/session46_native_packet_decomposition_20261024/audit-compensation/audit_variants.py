#!/usr/bin/env python3
"""Independent native and full-image validation of event/work Pareto route variants."""
import argparse,csv,gzip,hashlib,json
from collections import Counter
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,required=True);p.add_argument('--packet',type=Path,required=True);p.add_argument('--variants',type=Path,required=True);p.add_argument('--images',type=Path,required=True);p.add_argument('--base-images',type=Path,nargs=3,required=True);p.add_argument('--routes',type=Path,required=True);p.add_argument('--optimization',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
checks=0
def require(ok,why):
 global checks
 checks+=1
 if not ok:raise RuntimeError(why)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
R=a.repo/'05-knowledge/results';edgepath=R/'metagraph_order8_20261023/merged_edges.tsv.gz';electrical=R/'metagraph_order8_electrical_20261023/receipt.json';paths={'packet':a.packet,'variants':a.variants,'variant_images':a.images,'audited_routes':a.routes,'optimization':a.optimization,'edges':edgepath,'electrical':electrical}
paths.update({f'base_images_{j}':f for j,f in enumerate(a.base_images)})
hashes={k:sha(v)for k,v in paths.items()};prior=json.loads(a.routes.read_text());op=json.loads(a.optimization.read_text());packet=json.loads(a.packet.read_text());variants=json.loads(a.variants.read_text());slots={z['source']:z for z in packet['source_options']}
require(prior['status']=='PASS'and all(prior['input_sha256'][k]==hashes[k]for k in('packet','edges','electrical')),'prior full native/current audit exactly matches originals')
require(op['status']=='PASS'and op['input_sha256']['audited_complete_routes']==hashes['audited_routes']and op['input_sha256']['packet']==hashes['packet'],'exact complete-graph optimization tied to same packet and audit')
known=set()
for key,f in zip(('images','extra_images','extra37_images'),a.base_images):
 require(prior['input_sha256'][key]==sha(f),'original image vectors unchanged since independent full multiplication')
 known.update(tuple(sorted(map(tuple,z['charge'])))for z in json.loads(f.read_text())['solutions'])
edges=[tuple(int(r[k])for k in('source','target','blue','black'))for r in csv.DictReader(gzip.open(edgepath,'rt'),delimiter='\t')];N=3528;U,V=1965,2234
E=json.loads(electrical.read_text())['order8'];X=E['higher_layers'][-1]['lifts'][0];chars=[[64*y for y in yy]for yy in E['radical_lifts_mod4'][:4]]+[X]
def lap(z):
 out=[0]*N
 for i,j,b,k in edges:v=(b+k)*(z[i]-z[j]);out[i]+=v;out[j]-=v
 return out
images=json.loads(a.images.read_text());require(images['input_sha256']==hashes['edges']and images['modulus']==256,'new images pinned to original full graph')
for sol in images['solutions']:
 z=sol['potential'];q=dict(sol['charge']);w=dict(sol['integer_correction_nonzero']);require(len(z)==N and z[-1]==0 and all(type(v)is int and 0<=v<256 for v in z),'literal full grounded image vector')
 require(sum(q.values())==sum(w.values())==0,'zero sums in exact image equation')
 for i,v in enumerate(lap(z)):require(v==q.get(i,0)+256*w.get(i,0),'every physical row of new exact integer image equation')
 known.add(tuple(sorted(q.items())))
full=(1<<21)-1;chords=[(i,j)for i in range(8)for j in range(i+2,8)];idx={e:k for k,e in enumerate(chords)}
def resource(w):
 rev=sum(((w>>idx[7-j,7-i])&1)<<k for k,(i,j)in enumerate(chords));return{w,w^full,rev,rev^full}
def current(ss,mask,chi):
 out=[0]*N
 for j,s in enumerate(ss):
  t=slots[s]['options'][mask>>j&1]['target'];v=2*(chi[s]-chi[t]);out[s]+=v;out[t]-=v
 return out
reports=[]
for route in variants['routes']:
 ss=route['sources'];old=next(z for z in prior['routes']if z['sources']==ss and not z['name'].startswith('no_buffer'));optimization=next(z for z in op['routes']if z['name']==old['name']);path=route['path'];start,end=path[0],path[-1]
 require((start,end)==(old['advertised_path'][0],old['advertised_path'][-1]),'exact same target endpoints and restored buffers')
 require(all(z in old['accepted_states']for z in path),'all variant states inside independently complete legal fibre')
 baseline=[current(ss,start,chi)for chi in chars];initial_targets=Counter(slots[s]['options'][start>>j&1]['target']for j,s in enumerate(ss));seen=set();previous=set();reentries=[]
 for index,mask in enumerate(path):
  actualports=set();resources=set()
  for j,s in enumerate(ss):
   opt=slots[s]['options'][mask>>j&1];T=resource(opt['mask']);require(len(T)==4 and not actualports.intersection(T),'no simultaneous raw-port or normalized-resource reuse');actualports.update(T);resources.add(min(T))
   require(all(opt['physical_path'][clock]==owner for clock,owner in slots[s]['fixed_physical_clock_owners']),'every actual physical owner clock retained')
  require(Counter(slots[s]['options'][mask>>j&1]['target']for j,s in enumerate(ss))==initial_targets,'every target margin at every intermediate')
  for k,chi in enumerate(chars):require(all((x-y)%256==0 for x,y in zip(current(ss,mask,chi),baseline[k])),'every original physical nodewise current at every intermediate')
  if index:reentries.append({'to_mask':mask,'reentered_resources':sorted((resources-previous)&seen)})
  seen.update(resources);previous=resources
 steps=[]
 for x,y in zip(path,path[1:]):
  changed=[j for j in range(len(ss))if(x^y)>>j&1];q=tuple(sorted((ss[j],2*((y>>j&1)-(x>>j&1)))for j in changed));require(q in known,'each new performed step has independently multiplied exact source image vector')
  clock=sum(sum(v!=w for v,w in zip(slots[ss[j]]['options'][x>>j&1]['physical_path'],slots[ss[j]]['options'][y>>j&1]['physical_path']))for j in changed);steps.append({'slot_support':len(changed),'owner_clock_changes':clock,'source_charge':q})
 require([s['slot_support']for s in steps]==route['step_sizes']and max(s['slot_support']for s in steps)==route['budget'],'advertised step support accounting')
 metric='at_optimal_clock_budget'if route['name'].endswith('min_clock')else'at_optimal_slot_budget';goal=optimization[metric];work=sum(s['owner_clock_changes']for s in steps);events=len(steps)
 require(events==goal['minimum_events']and work==goal['minimum_total_clock_work'],'route attains independently proven minimum event count and minimum clock work at declared optimal threshold')
 if metric=='at_optimal_clock_budget':require(max(s['owner_clock_changes']for s in steps)==goal['threshold'],'route attains exact optimal clock bottleneck')
 else:require(max(s['slot_support']for s in steps)==goal['threshold'],'route attains exact optimal slot bottleneck')
 reports.append({'name':route['name'],'path':path,'steps':steps,'maximum_slots':max(s['slot_support']for s in steps),'maximum_owner_clocks':max(s['owner_clock_changes']for s in steps),'total_owner_clock_work':work,'events':events,'resource_reentry_count':sum(len(x['reentered_resources'])for x in reentries),'resource_reentry_ledger':reentries,'optimization_metric':metric})
report={'status':'PASS','checks':checks,'input_sha256':hashes,'new_exact_image_vectors':len(images['solutions']),'verified_total_charge_vectors_available':len(known),'routes':reports,'scope':'Three new variants of already exhaustively audited native fibres. Every intermediate full physical current/margin/pin/resource and every new exact image vector checked. Optimizations are menu-specific and reversible; none imposes strict temporal nonreturn.'}
a.out.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');print(json.dumps({'status':'PASS','checks':checks,'new_images':len(images['solutions']),'routes':[{k:r[k]for k in('name','events','maximum_slots','maximum_owner_clocks','total_owner_clock_work','resource_reentry_count')}for r in reports]},sort_keys=True))
