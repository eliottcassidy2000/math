#!/usr/bin/env python3
"""Independent semantic hostile controls on isolated physical copies."""
import concurrent.futures
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
HERE=Path(__file__).resolve().parent
CASES=[
 ('omit_distinct_degree_four','if d not in (0,3):','if d not in (0,3,4):','a.fives()','five true strong class histogram'),
 ('include_regular_eleven_presentations','if d not in (0,3):','if d != 0:','a.fives()','five literal defect'),
 ('include_coalesced_nine_cell',"if name=='cyclic' and a==7 and b==0:coalesced+=1;continue",'if False:coalesced+=1;continue','a.triangles()','true triangle skeleton has no residual permutation'),
 ('forget_cyclic_skeleton_quotient','qaut=[p for p in perms3 if relabel(q,p)==q]','qaut=[(0,1,2)]','a.triangles()','seventeen cyclic strong skeletons minus one'),
 ('use_wrong_SCC_index','square[11-k]*e','square[min(4,12-k)]*e','a.run()','SCC padding coefficient10314'),
 ('erase_full_orbit_error','return sum(size[o[0]]-len(o) for o in orbits),tuple(sorted(map(len,orbits)))','return 0,tuple(sorted(map(len,orbits)))','a.fives()','five literal defect'),
]
def check(item):
 name,old,new,command,expected=item;answers=[]
 for optimized in (False,True):
  with tempfile.TemporaryDirectory(prefix='s51-eleven-independent-hostile-') as tmp:
   root=Path(tmp)
   for file in ('audit.py','SOURCE_PINS.json'):
    shutil.copyfile(HERE/file,root/file)
   shutil.copytree(HERE/'inherited',root/'inherited')
   path=root/'audit.py';source=path.read_text()
   if source.count(old)!=1:raise ValueError('mutation anchor not unique: '+name)
   path.write_text(source.replace(old,new))
   if hashlib.sha256(path.read_bytes()).hexdigest()==hashlib.sha256((HERE/'audit.py').read_bytes()).hexdigest():raise ValueError('unchanged mutant')
   run=subprocess.run([sys.executable]+(['-O'] if optimized else [])+['-c','import audit as a; '+command],cwd=root,capture_output=True,text=True)
   if run.returncode==0 or ('ValueError: '+expected) not in run.stderr:raise ValueError('mutant not rejected as expected: '+name+'\n'+run.stderr)
   answers.append({'optimized':optimized,'semantic_rejection':expected,'pinned_script_hash_changed':True})
 return {'name':name,'runs':answers}
if __name__=='__main__':
 with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:result=list(pool.map(check,CASES))
 print(json.dumps({'status':'PASS','hostiles':result},sort_keys=True,indent=2))
