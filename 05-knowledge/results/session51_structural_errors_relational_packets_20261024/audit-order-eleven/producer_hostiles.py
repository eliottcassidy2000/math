#!/usr/bin/env python3
"""Independent mutations of the frozen producer, never the live producer."""
import concurrent.futures
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
HERE=Path(__file__).resolve().parent
SOURCE=HERE/'producer'
CASES=[
 ('omit_degree_two','for d in (1, 2):','for d in (1,):','c.check_five_sector()','all positive five decorated classes'),
 ('include_regular_collapse',"actual_orientations = [(0, 1), (1, 0)] if d == 1 else [(0, 1)]","actual_orientations = [(0, 1), (1, 0)]",'c.check_five_sector()','five decorated orbit error'),
 ('retain_coalesced_triangle_skeletons','if actual == (1, 1, 9):','if False:','c.check_triangles()','triangle actual and collapsed strong skeleton counts'),
 ('omit_triangle_block_permutations','for p in permutations(range(3)) for s in permutations(range(2)))','for p in [(0, 1, 2)] for s in permutations(range(2)))','c.check_triangles()','triangle skeleton free block orbit partition'),
 ('wrong_lower_SCC_coefficient','padded[e] += a2[2]*count','padded[e] += a2[3]*count','c.run()','known total eleven positive distribution'),
 ('erase_actual_orbit_error','return len(actions), sum(sizes[O[0]] - len(O) for O in orbits), tuple(sorted(map(len, orbits)))','return len(actions), 0, tuple(sorted(map(len, orbits)))','c.check_five_sector()','five decorated orbit error'),
]
def check(item):
 name,old,new,command,expected=item;answers=[]
 manifest=json.loads((SOURCE/'MANIFEST.json').read_text())
 for optimized in (False,True):
  with tempfile.TemporaryDirectory(prefix='s51-eleven-producer-hostile-') as tmp:
   root=Path(tmp)
   for e in manifest['files']:
    target=root/e['path'];target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(SOURCE/e['path'],target)
   shutil.copyfile(SOURCE/'MANIFEST.json',root/'MANIFEST.json')
   path=root/'check.py';source=path.read_text()
   if source.count(old)!=1:raise ValueError('mutation anchor not unique: '+name)
   path.write_text(source.replace(old,new))
   run=subprocess.run([sys.executable]+(['-O'] if optimized else [])+['-c','import check as c; '+command],cwd=root,capture_output=True,text=True)
   if run.returncode==0 or ('RuntimeError: '+expected) not in run.stderr:raise ValueError('mutant not rejected as expected: '+name+'\n'+run.stderr)
   verify=subprocess.run([sys.executable]+(['-O'] if optimized else [])+['verify_manifest.py'],cwd=root,capture_output=True,text=True)
   if verify.returncode==0 or not any(label+'check.py' in verify.stderr for label in ('size mismatch: ','digest mismatch: ')):raise ValueError('manifest failed to reject mutation: '+name+'\n'+verify.stderr)
   answers.append({'optimized':optimized,'semantic_rejection':expected,'manifest_rejection':True})
 return {'name':name,'runs':answers}
if __name__=='__main__':
 with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:result=list(pool.map(check,CASES))
 print(json.dumps({'status':'PASS','hostiles':result},sort_keys=True,indent=2))
