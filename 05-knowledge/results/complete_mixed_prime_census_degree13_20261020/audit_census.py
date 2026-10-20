#!/usr/bin/env python3
import ast,json,sys,time,subprocess
from pathlib import Path
from collections import Counter
import extension_census as e
m=e.m
ROOT=Path(__file__).resolve().parent
# Only the explicit proposed permutation constructors are inherited from the incoming list.
src=ast.parse((ROOT/'incoming/odd_symmetry_mixed_prime_census_20261019.py').read_text())
names={'cyc','perm_from_dict','frob_mult','classes'}
ns={};exec(compile(ast.Module(body=[x for x in src.body if isinstance(x,ast.FunctionDef) and x.name in names],type_ignores=[]),'incoming constructors','exec'),ns)
ORIGINAL=e.relative_records

def fast_records(n,H,gens,tournament):
 pairs,index,orbits,base=m.pair_data(n,H);N=e.normalizer(n,H,gens)
 coord={x:i for i,x in enumerate(orbits)};mult=Counter()
 for p in N:
  r=m.prep(p,pairs,index);im=tuple(coord[m.move(x,r)] for x in orbits)
  shift=m.move(base,r,True)^base;flip=sum(1<<i for i,x in enumerate(orbits) if shift&x)
  m.need(shift==sum(o for i,o in enumerate(orbits) if flip>>i&1),'affine shift exact')
  signbits=sum(1<<i for i,x in enumerate(orbits) if m.sign(x,r)==-1)
  mult[im,flip,signbits]+=1
 tag=f'{n}_{len(H)}_{len(orbits)}_{"T" if tournament else "G"}'
 inp=ROOT/(tag+'.actions');out=ROOT/(tag+'.orbits')
 with inp.open('w') as f:
  f.write(f'{len(orbits)} {len(mult)} {len(N)}\n')
  for (im,flip,sgn),w in sorted(mult.items()):f.write(' '.join(map(str,(flip,sgn,w)+im))+'\n')
 subprocess.run([str(ROOT/'cube_orbits'),str(inp),str(out),'T' if tournament else 'G'],check=True)
 rec=[]
 for line in out.read_text().splitlines():
  s,stab,size=map(int,line.split());mask=(base if tournament else 0)^sum(o for i,o in enumerate(orbits) if s>>i&1)
  rec.append({'mask':mask,'relative_stabilizer_order':stab,'orbit_size':size})
 return rec,pairs,orbits,N
m.relative_records=fast_records

def allcases():
 for n in range(7,14):
  for name,gens in ns['classes'](n):yield n,name,gens

def run():
 mode=sys.argv[1] if len(sys.argv)>1 else 'new13'
 rows=[];existing=None
 if mode=='large13':cases=[(n,s,g) for n,s,g in allcases() if n==13 and s in ('C15 (3,5,1^5)','F21 (7,1^6)')]
 elif mode=='small13':cases=[(n,s,g) for n,s,g in allcases() if n==13 and s not in ('C15 (3,5,1^5)','F21 (7,1^6)')]
 elif mode=='all_lower':cases=[(n,s,g) for n,s,g in allcases() if n<=12]
 elif mode=='control':cases=[(n,s,g) for n,s,g in allcases() if n<=10]
 else:raise RuntimeError(mode)
 for n,name,gens in cases:
  t=time.monotonic()
  if mode=='control':
   H=m.closure(gens,n)
   for tour in (False,True):
    a=fast_records(n,H,gens,tour)[0];b=ORIGINAL(n,H,gens,tour)[0]
    m.need(a==b,'accelerator equals original full records')
  r=m.census_case(n,name,gens);rows.append(r)
  print(n,name,{k:r[k] for k in ('order','normalizer_order','pair_orbits','relative_tournament_orbits','tau_T','tau_E','difference','tournament_fusion_defect','even_graph_fusion_defect','sign_escape_relative_orbits','sign_escape_global_types')},'seconds',round(time.monotonic()-t,2),flush=True)
  (ROOT/(mode+'.json')).write_text(json.dumps({'status':'RUNNING','checks':m.CHECKS,'rows':rows},indent=2,sort_keys=True)+'\n')
 (ROOT/(mode+'.json')).write_text(json.dumps({'status':'PASS','checks':m.CHECKS,'rows':rows},indent=2,sort_keys=True)+'\n')
if __name__=='__main__':run()
