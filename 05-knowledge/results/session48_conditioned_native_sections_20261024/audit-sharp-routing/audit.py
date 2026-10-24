#!/usr/bin/env python3
"""Independent blocked-owner-4 route audit. Explicit checks survive Python -O."""
from pathlib import Path
from itertools import combinations
from collections import defaultdict, Counter
from hashlib import sha256
import json, subprocess, sys, argparse
P=Path(__file__).resolve().parent
I=P/'inputs'
CHECKS=0

def need(x,why):
 global CHECKS
 CHECKS+=1
 if not x:raise RuntimeError(why)
def read(name):return json.loads((I/name).read_text())
def add(a,b):return (((a//16+b//16)%128)*16)|((a^b)&15)
def neg(a):return ((-(a//16))%128)*16|(a&15)
def charge(ss):
 q=0
 for s in ss:q=add(q,L[s])
 return q

def sources():
 for row in json.loads((P/'INPUT_PROVENANCE.json').read_text()):
  b=(P/row['path']).read_bytes();need(len(b)==row['bytes'] and sha256(b).hexdigest()==row['sha256'],'pinned source hash '+row['path'])
 old=read('native_packet_decomposition_20261023/native/domains.json')['sources']
 supply=read('compatible_charged_helpers_20261023/supply/selected_supply.json')['helpers']
 repair=read('guarded_component_completion_20261024/sections/native_2125_receipt.json')
 ee=read('metagraph_order8_electrical_20261023/receipt.json')['order8'];X=ee['higher_layers'][-1]['lifts'][0];Y=ee['radical_lifts_mod4']
 paths={r['source']:sorted(r['exact_all_shared_clocks_target_domain'],key=lambda t:t['target']) for r in old}
 for r in supply:paths[r['source']]=[r['selected']['before'],r['selected']['after']]
 paths[2125]=repair['records']
 core=sorted(r['source'] for r in old);free=sorted(set(paths)-set(core)-{368,776,3139})
 costs={s:sum(a!=b for a,b in zip(*(r['path'] for r in pp))) for s,pp in paths.items()}
 labels={s:(X[s]%128)*16+sum((Y[j][s]%2)<<(j-1) for j in range(1,5)) for s in paths}
 need(len(core)==28 and len(paths)==82 and len(free)==51,'authentic declared populations')
 native=json.loads((P/'native_certificate.json').read_text())
 need(core==native['cores'] and free==native['free_helpers'],'native auditor source identities')
 need(costs=={int(s):v for s,v in native['costs'].items()} and labels=={int(s):v for s,v in native['labels'].items()},'labels and weights independently derived from literal original records')
 oldlabels=[tuple(map(int,line.split())) for line in (I/'sixslot570_20261024/labels.txt').read_text().splitlines()]
 need(len(oldlabels)==28,'full inherited core label table')
 for s,(b,x) in zip(core,oldlabels):need(labels[s]==16*x+(b>>1)+(((x%2)^(b&1))<<3),'explicit additive old Y0..3 to Y1..4 adapter')
 bridge=read('sixslot570_20261024/bridges.json');need(bridge['sources']==core,'samecharge theorem same source identities')
 def maskcharge(m):
  q=0
  for j,s in enumerate(core):
   if m>>j&1:q=add(q,labels[s])
  return q
 def maskwork(m):return sum(costs[s] for j,s in enumerate(core) if m>>j&1)
 need(len(bridge['moves'])==220,'all inherited completion edges')
 for a,b in bridge['moves']:
  need(not a&b and a.bit_count()==b.bit_count()==4 and maskcharge(a)==maskcharge(b) and maskwork(a|b)<=22,'inherited eight-slot edge physical labels and exact work')
 small=set()
 for k in (2,3):
  bank=defaultdict(list)
  for sub in combinations(range(28),k):
   m=sum(1<<j for j in sub);bank[maskcharge(m)].append(m)
  for values in bank.values():
   for a,b in combinations(values,2):
    u,v=a&~b,b&~a;small.add(tuple(sorted((u,v))))
 need(len(small)==1944 and max(maskwork(a|b) for a,b in small)==20,'complete inherited four/six-slot types actual cost20')
 # Frozen census is an inherited theorem receipt, not a new rerun of 2^28 states.
 census=[]
 for budget in (6,8):
  a=[json.loads(t) for t in (I/f'sixslot570_20261024/fibres{budget}.jsonl').read_text().splitlines()]
  b=[json.loads(t) for t in (I/f'sixslot570_20261024/audit_fibres{budget}.jsonl').read_text().splitlines()]
  need(len(a)==len(b)==(15 if budget==6 else 2),'inherited complete grade range')
  for u,v in zip(a,b):
   need(u['rows']==v['rows'],'inherited independently enumerated census rows agree')
   need(sum(r[2]>1 for r in u['rows'])==({4:197,5:11}.get(u['k'],0) if budget==6 else 0),'inherited exceptional grades')
   census.append({k:u[k] for k in ('k','budget','states','disconnected')})
 return core,free,costs,labels,paths,census

C,H,W,L,PATHS,CENSUS=sources()

def cover(k,max_expensive):
 cheap=[s for s in C if W[s]==2];expensive=sorted(set(C)-set(cheap));need(len(cheap)==20,'twenty cheap records')
 rows=[0]*2048;stages=[];charts=0;entries=0
 for e in range(max_expensive+1):
  for E in combinations(expensive,e):
   # Chart fixes E and k-e-3 cheap sources. Every pair changes <=6 cheap cells.
   for F in combinations(cheap,k-e-3):
    fixed=E+F;remaining=sorted(set(cheap)-set(F));base=charge(fixed)
    image={add(base,charge(T)) for T in combinations(remaining,3)}
    bits=0
    for q in image:bits|=1<<q
    for q in image:rows[q]|=bits
    charts+=1;entries+=len(image)
  stages.append({'fixed_expensive_up_to':e,'uncovered_ordered_charge_pairs':sum(2048-r.bit_count() for r in rows),'charts':charts,'charge_image_entries':entries})
 need(all(r==(1<<2048)-1 for r in rows),'all ordered charge pairs covered at grade '+str(k))
 return {'grade':k,'stages':stages,'covered_ordered_pairs':2048**2,'max_core_support':6,'max_core_clock':12,'rows_sha256':sha256(b''.join(r.to_bytes(256,'little') for r in rows)).hexdigest()}

def work(mask,ids):return sum(W[s] for j,s in enumerate(ids) if mask>>j&1)
def complete_low():
 output=[];hostiles=[];ids=C+H;index={s:i for i,s in enumerate(ids)}
 for k in (0,1):
  for sign in (1,-1):
   bank=defaultdict(list)
   for c in combinations(C,k):
    for h in combinations(H,2):
     q=add(charge(c) if sign==1 else neg(charge(c)),charge(h));m=sum(1<<index[s] for s in c+h)
     bank[q].append((m,c,h))
   hist=Counter();trees_digest=sha256()
   for q,states in sorted(bank.items()):
    n=len(states);dist=[999]*n;dist[0]=0;prev=[-1]*n;seen=set();tree=[]
    for z in range(n):
     a=min((i for i in range(n) if i not in seen),key=lambda i:(dist[i],i));seen.add(a)
     if prev[a]>=0:tree.append((prev[a],a,dist[a]))
     for b in range(n):
      if b in seen:continue
      w=work(states[a][0]^states[b][0],ids)
      if w<dist[b]:dist[b]=w;prev[b]=a
    bottleneck=max(dist);hist[bottleneck]+=1;trees_digest.update(json.dumps([q,tree],separators=(',',':')).encode())
    need(bottleneck<=33,'complete low fibre has 33-clock spanning tree')
    if bottleneck==33:
     # Complete adjacency matrix, not a proposed alphabet or sheet quotient.
     matrix=[[work(a[0]^b[0],ids) for b in states] for a in states]
     unused=set(range(n));components=[]
     while unused:
      todo=[min(unused)];unused.remove(todo[0])
      for a in todo:
       hits=[b for b in unused if matrix[a][b]<=32]
       unused.difference_update(hits);todo.extend(hits)
      components.append(sorted(todo))
     need(len(components)==2,'two complete low-grade components at 32')
     minimum=min(matrix[a][b] for a in components[0] for b in components[1]);need(minimum==33,'true full-graph minimum cut equals 33')
     physical=[]
     for m,c,h in states:
      v=sorted(c if sign==1 else set(C)-set(c));helpers=sorted(h+(776,3139));selected=set(v+helpers)
      need(368 not in selected and 776 in selected and 3139 in selected,'external owner4 literal phase')
      pathrows={str(s):PATHS[s][int(s in selected)]['path'] for s in sorted(PATHS)}
      clockvals=[pathrows['368'][1], [8,9,4,10,11,12,13,3][pathrows['776'][4]], [14,15,16,17,18,3,4,19][pathrows['3139'][7]]]
      need(clockvals==[4,4,4],'actual shared owner value is4, not abstract bit phase')
      physical.append({'core_V':v,'helper_V':helpers,'physical_total_charge':charge(v+helpers),'owner_values':clockvals,'full_paths':pathrows})
     need(len(set(r['physical_total_charge'] for r in physical))==1,'all full states in same physical charged fibre')
     need(all((a[0]^b[0]).bit_count()<=6 for a in states for b in states),'support cap cannot cause lowercut')
     hostiles.append({'virtual_grade':k,'sign':sign,'virtual_charge':q,'physical_core_grade':k if sign==1 else 28-k,'states':[{'core_virtual':c,'helper_free':h} for m,c,h in states],'physical_states':physical,'complete_clock_matrix':matrix,'components_at32':components,'minimum_crossing_clock':minimum,'spanning_tree_at33':tree})
   output.append({'grade':k,'sign':sign,'states':sum(map(len,bank.values())),'fibres':len(bank),'minimum_fibre_size':min(map(len,bank.values())),'bottleneck_histogram':dict(sorted(hist.items())),'spanning_trees_sha256':trees_digest.hexdigest()})
 need([(x['sign'],x['virtual_charge'],len(x['states'])) for x in hostiles]==[(1,1156,9),(-1,700,13)],'both independent sharp hostiles')
 return {'summary':output,'sharp_hostiles':hostiles}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');ap.add_argument('--skip-cover',action='store_true');args=ap.parse_args()
 txt='28 51\n'+'\n'.join(f'{s} {L[s]} {W[s]}' for s in C+H)+'\n';(P/'routing_input.txt').write_text(txt)
 subprocess.run(['g++','-std=c++17','-O2',str(P/'section_graph.cpp'),'-o',str(P/'section_graph')],check=True)
 section=json.loads(subprocess.check_output([str(P/'section_graph'),str(P/'routing_input.txt')],text=True))
 result={'status':'PASS','scope':'native82 fixed owner4, exact separate core/helper counts, original full nodewise dyadic bank, primary slot-clock occurrences','inherited_core_sheet_theorem':'samecharge (8,22); c09b sixslot570_20261024; full census retained, not rerun','inherited_census':CENSUS,'interior_covers':[cover(5,1),cover(6,0)] if not args.skip_cover else [],'boundary_section_graphs':section,'complete_small_fibres':complete_low(),'checks':CHECKS}
 result=json.loads(json.dumps(result))
 if args.write:(P/'audit_receipt.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
 else:need(result==json.loads((P/'audit_receipt.json').read_text()),'frozen audit receipt equality')
 print(json.dumps({'status':'PASS','checks':CHECKS,'interior_ordered_pairs':2048**2,'clock_bound':33,'support_bound':8,'sharp_physical_core_grades':[1,27]},sort_keys=True))
if __name__=='__main__':main()
