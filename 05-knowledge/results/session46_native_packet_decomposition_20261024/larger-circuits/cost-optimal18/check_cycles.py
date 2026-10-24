#!/usr/bin/env python3
"""Literal relative-owner cycles and complete weighted menu polynomial."""
from pathlib import Path
from collections import Counter
from itertools import combinations
from math import comb
import json,csv,argparse,hashlib
H=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,required=True);a=p.parse_args()
P=json.loads((H/'pool.json').read_text());C=json.loads((H/'candidates.json').read_text())['found']['9'];O={r['source']:r for r in P['options']};L={int(s):tuple(q) for s,q in P['labels'].items()};checks=0
def need(ok,why):
 global checks
 checks+=1
 if not ok:raise RuntimeError(why)
records=[];costs=[];ids=C['plus']+C['minus'];signs=[1]*9+[-1]*9
for s,sign in zip(ids,signs):
 r=O[s];old=r['before' if sign==1 else 'after']['path'];new=r['after' if sign==1 else 'before']['path'];perm=dict(zip(old,new));seen=set();cycles=[]
 for owner in range(8):
  if owner in seen:continue
  cycle=[];v=owner
  while v not in seen:seen.add(v);cycle.append(v);v=perm[v]
  need(v==owner,'literal relative-owner map is a permutation')
  if len(cycle)>1:cycles.append(cycle)
 cost=sum(len(c) for c in cycles);need(cost==sum(x!=y for x,y in zip(old,new))==r['cost'],'cycle support equals changed owner clocks')
 need([len(c) for c in cycles]==([3] if s==2835 else [2]),'seventeen transpositions and one three-cycle')
 costs.append(cost);records.append({'source':s,'old_path':old,'new_path':new,'old_owner_to_new_owner_at_same_clock':perm,'nontrivial_cycles':cycles,'clock_cost':cost})
need(len(records)==18 and costs.count(2)==17 and costs.count(3)==1,'exact cycle census')
r=next(r for r in records if r['source']==2835);need(r['nontrivial_cycles']==[[0,4,6]],'exact relative owner cycle(0 4 6)')
third=[r['old_owner_to_new_owner_at_same_clock'][v] for v in r['new_path']]
need(third==[2,1,7,5,4,0,6,3],'third ambient cycle iterate')
name='05-knowledge/results/metagraph_order8_20261023/merged_nodes.tsv';raw=(a.repo/name).read_bytes();pins=json.loads((H/'SOURCE_PINS.json').read_text());need(hashlib.sha256(raw).hexdigest()==pins['originals'][name],'original matrix source pin')
nodes=list(csv.DictReader(raw.decode().splitlines(),delimiter='\t'));code=int(nodes[2835]['canon_a']);B=[[0]*8 for _ in range(8)]
for bit,(i,j) in enumerate(combinations(range(8),2)):B[i][j]=(code>>(27-bit))&1;B[j][i]=1-B[i][j]
bad=[[u,v] for u,v in zip(third,third[1:]) if not B[u][v]];need(bool(bad),'ambient order3 permutation is not a three-state native action')
need(all(B[u][v] for p in [r['old_path'],r['new_path']] for u,v in zip(p,p[1:])),'both selected source2835 paths are native')
# Complete bivariate polynomial by direct assignment and full-Q condition.
n=18;count=[0]*(1<<n);mass=[0]*(1<<n);parity=[0]*(1<<n);dist=[0]*(1<<n);slots=[0]*(1<<n);rawpoly=Counter();guardpoly=Counter();levels={m:Counter() for m in [2,4,8,16,32,64,128,256]};odd_counts=Counter();odd_slot_counts=Counter();odd_slot=ids.index(2835)
for mask in range(1,1<<n):
 low=mask&-mask;j=low.bit_length()-1;prev=mask^low
 count[mask]=count[prev]+signs[j];mass[mask]=(mass[prev]+signs[j]*L[ids[j]][1])%128;parity[mask]=parity[prev]^L[ids[j]][0];dist[mask]=dist[prev]+costs[j];slots[mask]=slots[prev]+1
for mask in range(1<<n):
 rawpoly[slots[mask],dist[mask]]+=1
 if count[mask]==mass[mask]==parity[mask]==0:guardpoly[slots[mask],dist[mask]]+=1
 if count[mask]==0:
  for m,poly in levels.items():
   if m==2 or (parity[mask]==0 and mass[mask]%(m//2)==0):
    poly[slots[mask],dist[mask]]+=1;odd_counts[m]+=dist[mask]%2;odd_slot_counts[m]+=(mask>>odd_slot)&1
expected=Counter()
for a2 in range(18):
 for a3 in range(2):expected[a2+a3,2*a2+3*a3]+=comb(17,a2)
need(rawpoly==expected,'raw polynomial is(1+z*t^2)^17*(1+z*t^3)')
need(guardpoly==Counter({(0,0):1,(18,37):1}),'full guard polynomial is1+z^18*t^37')
tower=[]
for m,poly in levels.items():
 need(all(poly.get((18-s,37-c),0)==v for (s,c),v in poly.items()),'complement reciprocity at modulus'+str(m))
 need(2*odd_counts[m]==sum(poly.values()) and odd_slot_counts[m]==odd_counts[m],'odd clock parity is exactly the source2835 bit, with half the states')
 need(sum((-1)**c*v for (s,c),v in poly.items())==0,'clock polynomial divisible by1+t')
 tower.append({'modulus':m,'accepted_states':sum(poly.values()),'odd_clock_cost_states':odd_counts[m],'even_clock_cost_states':sum(poly.values())-odd_counts[m],'source2835_flipped_states':odd_slot_counts[m],'reciprocity':'P(z,t)=z^18*t^37*P(1/z,1/t)','coefficients':[{'changed_slots':s,'changed_clocks':c,'coefficient':v} for (s,c),v in sorted(poly.items())]})
need([r['accepted_states'] for r in tower]==[48620,1526,754,350,152,72,34,2],'exact full-bank state tower')
print(json.dumps({'status':'PASS','checks':checks,'dyadic_polynomial_tower':tower,'relative_cycles':records,'source2835_third_ambient_iterate':third,'third_iterate_forbidden_arcs':bad,'raw_bivariate_polynomial':'(1+z*t^2)^17*(1+z*t^3)','raw_coefficients':[{'changed_slots':s,'changed_clocks':c,'coefficient':v} for (s,c),v in sorted(rawpoly.items())],'full_guard_polynomial':'1+z^18*t^37','full_guard_coefficients':[{'changed_slots':s,'changed_clocks':c,'coefficient':v} for (s,c),v in sorted(guardpoly.items())],'scope':'Relative permutations on each literal source owner set, not automorphisms or a global common-host permutation action.'},indent=2,sort_keys=True))
