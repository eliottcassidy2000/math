#!/usr/bin/env python3
import argparse,itertools,json,hashlib,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--output',type=Path);p.add_argument('--repo',type=Path,default=Path('session35-publication/clean-repo'));args=p.parse_args();checks=0

def req(x):
 global checks
 checks+=1
 if not x:raise ValueError('control '+str(checks))
def compile_paths(T,required,forbidden,pins={}):
 n=len(T);succ={};pred={}
 if required&forbidden:return set()
 for u,v in required:
  i,j=(u,v) if T[u][v] else (v,u)
  if i in succ and succ[i]!=j:return set()
  if j in pred and pred[j]!=i:return set()
  succ[i]=j;pred[j]=i
 blocks=[];seen=set()
 for x in range(n):
  if x in pred:continue
  b=[];v=x
  while v not in seen:
   seen.add(v);b.append(v)
   if v not in succ:break
   v=succ[v]
  blocks.append(b)
 if len(seen)!=n:return set()
 starts={}
 for k,b in enumerate(blocks):
  ss={pins[v]-j for j,v in enumerate(b) if v in pins}
  if len(ss)>1:return set()
  if ss:starts[k]=next(iter(ss))
 out=set()
 for order in itertools.permutations(range(len(blocks))):
  t=0;ok=True
  for k in order:
   if k in starts and starts[k]!=t:ok=False;break
   t+=len(blocks[k])
  if not ok:continue
  for k,l in zip(order,order[1:]):
   u,v=blocks[k][-1],blocks[l][0]
   if not T[u][v] or tuple(sorted((u,v))) in forbidden:ok=False;break
  if ok:out.add(tuple(v for k in order for v in blocks[k]))
 return out
n=4;pairs=list(itertools.combinations(range(n),2));perms=list(itertools.permutations(range(n)));cases=0;pin_cases=0
for code in range(1<<len(pairs)):
 T=[[0]*n for _ in range(n)]
 for k,(u,v) in enumerate(pairs):T[u][v]=(code>>k)&1;T[v][u]=1-T[u][v]
 legal=[P for P in perms if all(T[u][v] for u,v in zip(P,P[1:]))]
 masks={P:{tuple(sorted((u,v))) for u,v in zip(P,P[1:])} for P in legal}
 for flags in itertools.product(range(3),repeat=len(pairs)):
  positive={e for e,x in zip(pairs,flags) if x==1};negative={e for e,x in zip(pairs,flags) if x==2}
  direct={P for P in legal if positive<=masks[P] and not negative&masks[P]}
  compiled=compile_paths(T,positive,negative);req(compiled==direct);cases+=1
  # Every single-global-clock pin is separately tested on one source in eight;
  # no inference that local cell clocks are globally aligned.
  if code%8==0:
   for v in range(n):
    for c in range(n):
     pin={v:c};got=compile_paths(T,positive,negative,pin)
     req(got=={P for P in direct if P[c]==v});pin_cases+=1
  if code%8==0:
   converse=list(map(list,zip(*T)))
   req(compile_paths(converse,positive,negative)=={P[::-1] for P in direct})
# Incoming actual shared-owner pair, independent of unspecified cross arcs.
coreU=[5,7,2,3,1,4,6,0];coreV=[5,7,6,3,1,4,2,0];rho=[8,9,10,2,6,11,12,13]
helperU=[rho[i] for i in [7,4,6,5,3,0,2,1]];helperV=[rho[i] for i in [7,3,6,1,5,4,0,2]]
rel='05-knowledge/results/compatible_charged_helpers_20261023/supply/shared_owner.json'
source=args.repo/rel
raw=source.read_bytes() if source.exists() else subprocess.check_output(['git','-C',str(args.repo),'show','HEAD:'+rel])
data=json.loads(raw)
req(data['two_cell_records']['176'][0]['global_path']==coreU)
req(data['two_cell_records']['176'][1]['global_path']==coreV)
req(data['two_cell_records']['368'][0]['global_path']==helperU)
req(data['two_cell_records']['368'][1]['global_path']==helperV)
req(data['embeddings']['368']==rho)
for a,b,branch in [(coreU,helperV,2),(coreV,helperU,6)]:
 edges=set(zip(a,a[1:]))|set(zip(b,b[1:]));req({v for u,v in edges if u==branch}=={3,12})
# A semicomplete block carrier can still fail a separate global-clock demand.
T=[[0,1,1],[0,0,1],[0,0,0]]
req(compile_paths(T,set(),set())=={(0,1,2)})
req(compile_paths(T,set(),set(),{2:0})==set())
out={'source_sha256':hashlib.sha256(raw).hexdigest(),'status':'PASS','checks':checks,'all_four_owner_flag_instances':cases,'single_global_pin_instances':pin_cases,'native_host_obstructions':[{'core':'U','helper':'V','branch':2,'successors':[3,12]},{'core':'V','helper':'U','branch':6,'successors':[3,12]}],'scope':'Exact all-flag four-owner compiler controls, selected one-global-clock and converse controls, plus incoming local-host global-observer obstruction.'}
s=json.dumps(out,indent=2,sort_keys=True)+'\n'
if args.output:args.output.write_text(s)
else:print(s,end='')
