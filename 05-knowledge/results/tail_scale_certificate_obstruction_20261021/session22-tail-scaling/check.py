#!/usr/bin/env python3
import argparse,hashlib,json,importlib.util
from pathlib import Path
base=None
# Explicit base path keeps staging and published layouts equally reproducible.
CHECKS=0
def ck(x,msg):
 global CHECKS
 CHECKS+=1
 if not x:raise ValueError(msg)
def prefix(s,h,V):
 u=s;t=0;path=[]
 for j in range(h):
  ck(u not in V,'literal prefix remains unfinished')
  v,c=base.edge(u);path.append([u,c]);u=v;t+=c
 ck(u not in V,'literal final cutoff remains unfinished')
 return u,t

def main():
 global base
 p=argparse.ArgumentParser();p.add_argument('--base-dir',type=Path,required=True);p.add_argument('--repo-root',type=Path,required=True);p.add_argument('--base-receipt',type=Path,required=True);p.add_argument('--output',type=Path);p.add_argument('--verify',type=Path);a=p.parse_args()
 spec=importlib.util.spec_from_file_location('routefree_base',a.base_dir/'check.py');base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
 raw=(a.repo_root/'05-knowledge/results/session18_native_annulus/frontier.json').read_bytes().replace(b'\r\n',b'\n');ck(hashlib.sha256(raw).hexdigest()==base.FRONTIER_SHA,'frontier pin');V=set(json.loads(raw)['frontier'])
 br=a.base_receipt.read_bytes();ck(hashlib.sha256(br).hexdigest()=='e5ec575906451e4a6ee293a98ed85dc3d90a2a8b18b878a921df151166d9032a','base all-horizon receipt pin')
 b=next(c for c in json.loads(br)['cases']if c['N']==4096)
 N=64*b['N'];h=b['horizon']+6;r=base.PrefixRoutes(V);rows=[];counts=[0]*12;scanned=0
 for s in range(128*N+1,256*N,2):
  if not base.primitive(s):continue
  scanned+=1
  if r.get(s,h)is not None:continue
  v=base.charge(64*N,s);ck(base.Q(v)==1,'upper-band unit pending charge');end,t=prefix(s,h,V)
  counts[base.INDEX[v]]+=1;rows.append([s,s%15,s%35,*v,end,t])
  if len(rows)==4096:break
 ck(len(rows)==4096,'explicit pending subset completed')
 D=base.zonotope(counts);ck(D>64*b['D_new'],'refute linear one-scale D propagation')
 ck(len(rows)>8*b['pending_nonzero'],'refute square-root pending-count propagation')
 # Exact source lift and its missing-density coordinate.
 images=set()
 for s in range(1,513,2):
  if not base.primitive(s):continue
  u,t=base.edge(s);aa=2 if s%8==1 else 1;c=7*(2**aa+3);z=64*s+c
  ck(base.edge(z)==(u,t+6),'native inverse-address lift: same parent, clock plus six')
  images.add(z%512)
 ck(images=={113,227,483},'exact image residue classes')
 # One lost-magnitude hostile with unchanged source charge, guards and edge clock.
 x,y=32881,247921;u,t=base.edge(x);v,tp=base.edge(y)
 ck(y-x==840*2**t and v-u==2520,'literal translation cylinder')
 ck((u,v,t,tp)==(385,2905,8,8),'first-edge clock and endpoints')
 ck(x%15==y%15 and x%35==y%35,'same joint guards')
 ck(base.charge(64*4096,x)==base.charge(64*4096,y)==(-1,-1),'same nonzero source charge')
 ck(x not in V and y not in V and u in V and v not in V,'different next first-arrival predicate')
 translations=[]
 for s in (3,7,27,97):
  for length in (1,2,4,8):
   states=[s];clocks=[]
   for j in range(length):
    z,c=base.edge(states[-1]);states.append(z);clocks.append(c)
   T=sum(clocks);shift=840*2**T;z=s+shift;cum=0
   for j,c in enumerate(clocks):
    ck(z-states[j]==840*3**j*2**(T-cum),'exact translated prefix state')
    zz,cc=base.edge(z);ck(cc==c,'translated full native clock');ck(z%105==states[j]%105,'translated joint source guard');cum+=c;z=zz
   ck(z-states[-1]==840*3**length,'translated final endpoint');translations.append([s,length,T,shift])
 # Exact six-antipodal-axis compression and conditional cone controls.
 axis_checks=0
 for seed in range(100):
  cc=[((seed+3)*(i+5)+i*i)%7 for i in range(12)]
  axes=[cc[i]+cc[i+6]for i in range(6)];rep=axes+[0]*6
  ck(base.zonotope(cc)==base.zonotope(rep),'signed-zonotope antipodal axis compression')
  Ts=[(seed%11-5,seed%7-3),(2,-1),(-3,4)]
  ck(base.linear(cc,Ts)==base.linear(rep,Ts),'terminal support antipodal axis compression')
  small=[8*z-(seed%(z+1)if z else 0)for z in axes]
  ck(base.zonotope(small+[0]*6)<=64*base.zonotope(rep),'conditional scale cone zonotope inclusion')
  ell=[max(abs(base.beta(t,base.VECTORS[i]))for t in Ts)for i in range(6)]
  ellp=[8*z-(seed%(z+1)if z else 0)for z in ell]
  ck(sum(x*y for x,y in zip(small,ellp))<=64*sum(x*y for x,y in zip(axes,ell)),'conditional linear support cone')
  axis_checks+=1
 base_axes=[b['pending_new_counts'][i]+b['pending_new_counts'][i+6]for i in range(6)]
 subset_axes=[counts[i]+counts[i+6]for i in range(6)]
 ck(any(subset_axes[i]>8*base_axes[i]for i in range(6)),'actual subset violates naive count cone')
 result={'status':'PROVED native scale-image and translation-cylinder obstruction; FINITE-EXACT subset refutes naive D contraction; uniform153N OPEN','base_receipt_sha256':hashlib.sha256(br).hexdigest(),'N':N,'horizon':h,'base_N':4096,'base_horizon':105,'base_pending':b['pending_nonzero'],'base_D':b['D_new'],'candidate_D_limit':64*b['D_new'],'pending_subset_size':len(rows),'primitive_candidates_scanned':scanned,'pending_subset_D':D,'pending_direction_counts':counts,'base_axis_totals':base_axes,'pending_subset_axis_totals':subset_axes,'axis_cone_controls':axis_checks,'pending_rows_schema':['source','mod15','mod35','charge_a','charge_b','state_after_h','actual_clock_after_h'],'pending_rows':rows,'image_residues_mod512':sorted(images),'guard_clock_hostile':[x,y,u,v,t,tp,-1,-1,x%15,x%35],'translation_controls':translations,'checks':CHECKS}
 data=(json.dumps(result,sort_keys=True,separators=(',',':'))+'\n').encode()
 if a.output:a.output.write_bytes(data)
 if a.verify:ck(a.verify.read_bytes()==data,'immutable replay')
 print(json.dumps({k:v for k,v in result.items()if k not in ('pending_rows','translation_controls')},sort_keys=True));print(json.dumps({'checks':CHECKS,'sha256':hashlib.sha256(data).hexdigest()}))
if __name__=='__main__':main()
