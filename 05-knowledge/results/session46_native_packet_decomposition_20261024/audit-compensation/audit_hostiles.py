#!/usr/bin/env python3
"""Negative controls: mutate only temporary copies, expect explicit failures even with -O."""
import argparse,copy,json,subprocess,sys,tempfile
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,required=True);p.add_argument('--packet',type=Path,required=True);p.add_argument('--images',type=Path,required=True);p.add_argument('--foundation',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
packet=json.loads(a.packet.read_text());images=json.loads(a.images.read_text());cases=[]
x=copy.deepcopy(packet);x['source_options'][0]['options'][0]['physical_path'][0]=x['source_options'][0]['options'][0]['physical_path'][1];cases.append(('altered_physical_owner_lift','packet',x,'exact retained owner lift'))
x=copy.deepcopy(packet);x['source_options'][0]['options'][0]['mask']=x['source_options'][0]['options'][1]['mask'];cases.append(('resource_renamed_as_extra_capacity','packet',x,'actual option recovered'))
x=copy.deepcopy(packet);x['routes'][0]['path'][1]^=1;cases.append(('illegal_intermediate_margin_or_current','packet',x,'every advertised intermediate literal state accepted'))
x=copy.deepcopy(images);x['solutions'][0]['potential'][0]=(x['solutions'][0]['potential'][0]+1)%256;cases.append(('corrupted_full_laplacian_image','images',x,'every original physical row exact integer image equation'))
results=[]
with tempfile.TemporaryDirectory(prefix='native-compensation-hostiles-')as tmp:
 d=Path(tmp)
 for name,target,content,expected in cases:
  f=d/(name+'.json');f.write_text(json.dumps(content));cmd=[sys.executable]
  if sys.flags.optimize:cmd+=['-O']
  cmd += [str(Path(__file__).with_name('audit_routes.py')),'--repo',str(a.repo),'--packet',str(f if target=='packet'else a.packet),'--foundation',str(a.foundation),'--images',str(f if target=='images'else a.images),'--out',str(d/'unused.json')]
  completed=subprocess.run(cmd,text=True,capture_output=True)
  if completed.returncode==0 or expected not in completed.stderr:raise RuntimeError('Hostile did not fail at expected guard: '+name+' '+completed.stderr)
  results.append({'control':name,'status':'REJECTED','expected_guard':expected})
report={'status':'PASS','controls':results,'scope':'Temporary corrupted copies only; no producer input modified. All guards are explicit and remain live under optimized Python.'}
a.out.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');print(json.dumps(report,sort_keys=True))
