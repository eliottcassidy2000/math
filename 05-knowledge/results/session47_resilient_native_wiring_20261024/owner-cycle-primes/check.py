#!/usr/bin/env python3
"""Exact cycle/sign diagnostics. No native resource/authentication replay claim."""
import argparse,collections,hashlib,json,math,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,default=Path('session35-publication/clean-repo'));p.add_argument('--output',type=Path);a=p.parse_args()
r=a.repo
paths={'core':r/'05-knowledge/results/session46_native_packet_decomposition_20261024/pinned-domains/normal.json','helpers':r/'05-knowledge/results/compatible_charged_helpers_20261023/supply/selected_supply.json'}
def raw(path):
 if path.exists():return path.read_bytes()
 return subprocess.check_output(["git","-C",str(r),"show","HEAD:"+str(path.relative_to(r))])
checks=0
def check(x):
 global checks
 checks+=1
 if not x:raise ValueError('failed exact control '+str(checks))
def cycles(p,q):
 check(sorted(p)==list(range(8)) and sorted(q)==list(range(8)))
 f=dict(zip(p,q));seen=set();out=[]
 for i in range(8):
  if i in seen:continue
  c=[];j=i
  while j not in seen:seen.add(j);c.append(j);j=f[j]
  if len(c)>1:out.append(c)
 w=sum(len(c) for c in out)
 check(w==sum(x!=y for x,y in zip(p,q)))
 inv=sum(f[i]>f[j] for i in range(8) for j in range(i+1,8))
 check(inv%2==(w-len(out))%2)
 return out,w,inv%2
core=[]
for d in json.loads(raw(paths['core']))['domains']:
 ps=sorted((x for x in d['complete_pinned_paths'] if x['target_in_original_pair']),key=lambda x:x['target'])
 check(len(ps)==2)
 c,w,s=cycles(ps[0]['path'],ps[1]['path']);check(len(c)==1)
 core.append({'source':d['source'],'cycles':c,'clock_cost':w,'sign_bit':s})
helper=[]
for d in json.loads(raw(paths['helpers']))['helpers']:
 c,w,s=cycles(d['selected']['before']['path'],d['selected']['after']['path'])
 check(w==d['selected']['cost'])
 helper.append({'source':d['source'],'cycles':c,'clock_cost':w,'sign_bit':s})
# Every possible even changed-slot set has sum cycle count even: all core have one.
# Do not enumerate 2^28 subsets to establish this elementary implication.
check(len(core)==28 and all(len(x['cycles'])==1 for x in core))
A={1532,1691,1923,1974};B={1841,2176,2837,3360}
witness=[x for x in core if x['source'] in A|B]
check(len(witness)==8 and sum(x['clock_cost'] for x in witness)==22)
check(sum(x['sign_bit'] for x in witness)%2==0)
check(collections.Counter(x['clock_cost'] for x in witness)=={2:4,3:2,4:2})
# Basis swap: (a,b,c,d) -> (d,b,c,a) sends packed 1011 (11) to itself;
# instead explicit invertible XOR row map sends 1011 -> 0001.
# output bits (x0, x1+x0, x2, x3+x0), self-inverse.
def gauge(v):return (v&1)|((((v>>1)&1)^(v&1))<<1)|(v&4)|((((v>>3)&1)^(v&1))<<3)
check(gauge(11)==1)
for x in range(16):
 check(gauge(gauge(x))==x)
 for y in range(16):check(gauge(x^y)==gauge(x)^gauge(y))
# All images of generators in odd cyclic groups must be zero.
for odd in (3,11,33):
 for order in (2,128):
  check([v for v in range(odd) if order*v%odd==0]==[0])
# Explicit multiple-cycle helper breaks per-slot sign=(-1)^(w-1).
h=next(x for x in helper if x['source']==650)
check(sorted(map(len,h['cycles']))==[2,4])
check(h['sign_bit']!=(h['clock_cost']-1)%2)
result={'status':'PASS','checks':checks,'scope':'Cycle diagnostics from previously authenticated records; not a new native path/resource replay. No prime connection inferred.','input_sha256':{k:hashlib.sha256(raw(v)).hexdigest() for k,v in paths.items()},'core':core,'helpers':helper,'witness22':witness,'helper_cycle_histogram':dict(sorted(collections.Counter(','.join(map(str,sorted(map(len,x['cycles'])))) for x in helper).items())),'gauge':{'old_parity_label':11,'new_parity_label':1,'invertible':True},'odd_character_targets':[3,11,33]}
s=json.dumps(result,indent=2,sort_keys=True)+'\n'
if a.output:a.output.write_text(s)
else:print(s,end='')
