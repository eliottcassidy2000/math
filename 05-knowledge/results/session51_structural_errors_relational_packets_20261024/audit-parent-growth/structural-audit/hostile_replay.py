#!/usr/bin/env python3
"""Physically copied mutation controls; untouched source remains pinned."""
import concurrent.futures as cf,json,shutil,subprocess,sys,tempfile
from pathlib import Path
base=Path(__file__).resolve().parent;source=base/'producer';manifest=json.loads((source/'MANIFEST.json').read_text())
cases=[
 ('omit_degree_two','check.py','for d in (1, 2):','for d in (1,):','import check; check.check_five_cells()','five sector error'),
 ('omit_cell_permutations','check.py','for pi in permutations(range(3)):','for pi in (tuple(range(3)),):','import check; check.check_three_triangles()','triangle rigid class split'),
 ('omit_cyclic_majority','check.py','if not is_strong:','if not is_strong or cyclic:','import check; check.check_three_triangles()','constant attachment hostile actual partition'),
 ('zero_actual_error','check.py','return sum(size[o[0]] - len(o) for o in orbits)','return 0','import check; check.check_five_cells()','five actual orbit error'),
 ('omit_scc_side','check.py','total = es10 + 2 * es9_inherited','total = es10 + es9_inherited','import check; check.run()','SCC transfer order ten with inherited Es9'),
 ('wrong_distribution_padding','distribution_check.py','total10[e] += 2 * count','total10[e] += count','import distribution_check; distribution_check.run()','total ten exact defect distribution')]
def run(case):
 name,file,old,new,command,expected=case;answers=[]
 for optimized in (False,True):
  with tempfile.TemporaryDirectory(prefix='order10-hostile-') as tmp:
   dst=Path(tmp)
   for rec in manifest['files']:
    target=dst/rec['path'];target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source/rec['path'],target)
   shutil.copyfile(source/'MANIFEST.json',dst/'MANIFEST.json');p=dst/file;text=p.read_text()
   if text.count(old)!=1:raise RuntimeError('mutation anchor not unique: '+name)
   p.write_text(text.replace(old,new));args=[sys.executable]+(['-O'] if optimized else [])+['-c',command]
   out=subprocess.run(args,cwd=dst,capture_output=True,text=True)
   if out.returncode==0 or 'RuntimeError: '+expected not in out.stderr:raise RuntimeError('mutation not detected as intended: '+name+' '+out.stderr)
   verify=subprocess.run([sys.executable,'verify_manifest.py'],cwd=dst,capture_output=True,text=True)
   if verify.returncode==0 or 'manifest byte mismatch: '+file not in verify.stderr:raise RuntimeError('changed bytes passed manifest: '+name)
   answers.append({'optimized':optimized,'semantic_rejection':expected,'manifest_rejection':True})
 return {'name':name,'runs':answers}
with cf.ThreadPoolExecutor(max_workers=6) as ex: result=list(ex.map(run,cases))
print(json.dumps({'status':'PASS','hostiles':result},indent=2,sort_keys=True))
