#!/usr/bin/env python3
"""Independent decoder mutation tests; all writes occur in fresh copies."""
import concurrent.futures as cf,json,shutil,subprocess,sys,tempfile
from pathlib import Path
base=Path(__file__).resolve().parent;original=(base/'source/decoder.py').read_text()
cases=[
 ('forget_rigid_pins','pinned.update((i,j)); edges.append((i,j,\'rigid\'))','pinned.update(()); edges.append((i,j,\'rigid\'))'),
 ('pin_one_rigid_endpoint','pinned.update((i,j)); edges.append((i,j,\'rigid\'))','pinned.update((i,)); edges.append((i,j,\'rigid\'))'),
 ('ignore_cycle_frustration','balanced &= values[j]==v','balanced = True'),
 ('fail_inverse_gain','graph[j].append((i,pow(slope,-1,p)))','graph[j].append((i,slope))'),
 ('wrong_generator_rotation','g[C[x]]=C[(x+k)%p]','g[C[x]]=C[(x+1)%p]'),
 ('skip_cyclic_recognition','orders = [C if len(C)==1 else cyclic_order(rows,C) for C in cells]','orders = cells'),
 ('wrong_frozen_error','error=sum(len(cells[i])*(len(cells[i])-1) for i in frozen)','error=sum(len(cells[i]) for i in frozen)')]
def run(case):
 name,old,new=case
 if original.count(old)!=1:raise RuntimeError('mutation anchor not unique: '+name)
 results=[]
 for optimized in (False,True):
  with tempfile.TemporaryDirectory(prefix='cyclic-decoder-hostile-') as temp:
   dst=Path(temp)
   for file in ('independent.py','reference.py'):shutil.copyfile(base/file,dst/file)
   (dst/'mutant.py').write_text(original.replace(old,new))
   command=[sys.executable]+(['-O'] if optimized else [])+['independent.py','--quick','--decoder','mutant.py']
   out=subprocess.run(command,cwd=dst,capture_output=True,text=True)
   if out.returncode==0 or 'RuntimeError:' not in out.stderr:raise RuntimeError('mutation did not fail semantically: '+name+' '+out.stderr)
   reason=out.stderr.strip().splitlines()[-1]
   results.append({'optimized':optimized,'reason':reason})
 return {'name':name,'rejections':results}
with cf.ThreadPoolExecutor(max_workers=7) as ex:result=list(ex.map(run,cases))
print(json.dumps({'status':'PASS','mutants':result},indent=2,sort_keys=True))
