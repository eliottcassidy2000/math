#!/usr/bin/env python3
import concurrent.futures as F,json,shutil,subprocess,sys,tempfile
from pathlib import Path
HERE=Path(__file__).resolve().parent
CASES=[
 ('erase_binomial_multiplicity','weight*comb(multiplicity,out)','weight','grouped score count matches literal subset recurrence'),
 ('reverse_residual_degree_ownership','remaining+[degree]*out+[degree-1]*(multiplicity-out)','remaining+[degree-1]*out+[degree]*(multiplicity-out)','grouped score count matches literal subset recurrence'),
 ('omit_labelled_profile_factor','comb(10,5)*profile_count==labelled','profile_count==labelled','fixed root deletion labelled score profile identity'),
 ('duplicate_hidden_transitive_control',"[r for r in inherited['q3_candidates'] if r['error']]","inherited['q3_candidates']",'eighteen inherited nonrigid representatives without duplicate control'),
 ('confuse_group_order_with_orbits',"nonrigid_mass=sum(factorial(11)//r['full_automorphism_order'] for r in classes)","nonrigid_mass=sum(factorial(11)//len(r['vertex_orbit_sizes']) for r in classes)",'remaining regular labelled mass is rigid classes'),
 ('use_one_pointed_type_per_rigid_class','pointed=11*rigid+pointed_nonrigid','pointed=rigid+pointed_nonrigid','regular error from pointed class difference')]
def test(case):
 name,old,new,label=case;out=[];manifest=json.loads((HERE/'MANIFEST.json').read_text())
 for opt in (False,True):
  with tempfile.TemporaryDirectory(prefix='s51-score-mass-hostile-') as tmp:
   root=Path(tmp)
   for item in manifest['files']:
    p=root/item['path'];p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(HERE/item['path'],p)
   shutil.copyfile(HERE/'MANIFEST.json',root/'MANIFEST.json')
   p=root/'check.py';text=p.read_text()
   if text.count(old)!=1:raise ValueError('mutation anchor not unique: '+name)
   p.write_text(text.replace(old,new));flags=['-O'] if opt else []
   run=subprocess.run([sys.executable]+flags+['check.py'],cwd=root,capture_output=True,text=True)
   if run.returncode==0 or 'ValueError: '+label not in run.stderr:raise ValueError('semantic mutation survived: '+name+'\n'+run.stderr)
   verify=subprocess.run([sys.executable]+flags+['verify_manifest.py'],cwd=root,capture_output=True,text=True)
   if verify.returncode==0 or 'manifest byte mismatch: check.py' not in verify.stderr:raise ValueError('manifest accepted mutation: '+name)
   out.append({'optimized':opt,'semantic_rejection':label,'manifest_rejection':True})
 return {'name':name,'runs':out}
if __name__=='__main__':
 with F.ThreadPoolExecutor(max_workers=6) as pool:result=list(pool.map(test,CASES))
 print(json.dumps({'status':'PASS','hostiles':result},sort_keys=True,indent=2))
