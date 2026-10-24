#!/usr/bin/env python3
"""Physical-copy semantic mutations; normal and optimized checks must reject."""
import concurrent.futures as F
import json,shutil,subprocess,sys,tempfile
from pathlib import Path
HERE=Path(__file__).resolve().parent
CASES=[
 ('omit_one_stabilizer_slope','for slope in range(1,5):','for slope in range(1,4):','c.three_orbits()','forty line invariant degree three matrices'),
 ('erase_negative_matching_sign','((1,) if d==0 else (1,-1))','(1,)','c.seven_orbits()','five genuine q7 classes'),
 ('omit_second_q9_matching','for p in I.permutations(high_targets):','for p in list(I.permutations(high_targets))[:1]:','c.nine_orbits()','two q9 crossing matchings'),
 ('truncate_cyclic_universe','for bits in range(32)}','for bits in range(16)}','c.transitive_classes()','thirty two cyclic degree five connection sets'),
 ('count_hidden_Aut55_twice',"selected=q1+[r for r in q3 if r['error']]+q7+q9",'selected=q1+q3+q7+q9','c.run()','eighteen nonrigid regular eleven classes'),
 ('confuse_pointed_orbits_with_vertices',"pointed=sum(len(r['vertex_orbit_sizes']) for r in selected)",'pointed=sum(11 for r in selected)','c.run()','nonrigid pointed mass78 and error660')]
def test(case):
 name,old,new,command,label=case;results=[];manifest=json.loads((HERE/'MANIFEST.json').read_text())
 for opt in (False,True):
  with tempfile.TemporaryDirectory(prefix='s51-nonrigid-hostile-') as tmp:
   root=Path(tmp)
   for item in manifest['files']:
    target=root/item['path'];target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(HERE/item['path'],target)
   shutil.copyfile(HERE/'MANIFEST.json',root/'MANIFEST.json')
   path=root/'check.py';code=path.read_text()
   if code.count(old)!=1:raise ValueError('mutation anchor not unique: '+name)
   path.write_text(code.replace(old,new));flags=['-O'] if opt else []
   run=subprocess.run([sys.executable]+flags+['-c','import check as c; '+command],cwd=root,capture_output=True,text=True)
   if run.returncode==0 or 'ValueError: '+label not in run.stderr:raise ValueError('semantic mutation survived: '+name+'\n'+run.stderr)
   verify=subprocess.run([sys.executable]+flags+['verify_manifest.py'],cwd=root,capture_output=True,text=True)
   if verify.returncode==0 or 'manifest byte mismatch: check.py' not in verify.stderr:raise ValueError('mutation escaped pin: '+name+'\n'+verify.stderr)
   results.append({'optimized':opt,'semantic_rejection':label,'manifest_rejection':True})
 return {'case':name,'runs':results}
if __name__=='__main__':
 with F.ThreadPoolExecutor(max_workers=6) as pool:results=list(pool.map(test,CASES))
 print(json.dumps({'status':'PASS','hostiles':results},sort_keys=True,indent=2))
