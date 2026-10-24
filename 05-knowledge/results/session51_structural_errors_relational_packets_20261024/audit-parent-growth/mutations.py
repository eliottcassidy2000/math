#!/usr/bin/env python3
"""Physical semantic mutation controls, normal and NDEBUG; no removable asserts."""
from pathlib import Path
import subprocess,tempfile,json,hashlib
B=Path(__file__).resolve().parent
cases=[('label_dependent_singleton','colors[v]<colors[selected]','v<unsigned(selected)','complete order-nine labelled orbit mass'),('omit_residual_classes','}else{++without_singleton;','}else{++without_singleton;continue;','complete order-nine labelled orbit mass'),('drop_root_moving_cosets','for(unsigned root=0;root<9;++root)if(points[root].key==own)','for(unsigned root=8;root<9;++root)if(points[root].key==own)','complete order-nine labelled orbit mass'),('omit_attachment','mask<256','mask<255','all parent attachments covered')]
source=(B/'producer/build_nine.cpp').read_text();results=[]
with tempfile.TemporaryDirectory(prefix='parent-growth-mutations-') as tmp:
 D=Path(tmp)
 for name,old,new,failure in cases:
  if source.count(old)!=1:raise RuntimeError('mutation anchor '+name)
  altered=source.replace(old,new);p=D/(name+'.cpp');p.write_text(altered)
  for label,flags in [('O2',['-O2']),('O3-ndebug',['-O3','-DNDEBUG'])]:
   exe=D/(name+label);subprocess.run(['g++','-std=c++17',*flags,str(p),'-o',str(exe)],check=True,capture_output=True)
   r=subprocess.run([str(exe),str(B/'producer/parents-eight.txt'),str(D/'atlas.txt'),str(D/'receipt.json')],capture_output=True,text=True)
   if r.returncode==0 or failure not in r.stderr:raise RuntimeError(name+' unexpected '+str(r.returncode)+' '+r.stderr)
   results.append({'mutation':name,'mode':label,'exit_code':r.returncode,'failure':r.stderr.strip(),'changed_source_sha256':hashlib.sha256(altered.encode()).hexdigest()})
 # Corrupt the independently evaluated root size, which must fail at its first
 # per-parent histogram comparison, not merely global trace divisibility.
 source=(B/'incoming_trace.cpp').read_text();old='unsigned b=root(child);';new='unsigned b=3;';
 if source.count(old)!=1:raise RuntimeError('incoming mutation anchor')
 altered=source.replace(old,new);p=D/'bad-refinement.cpp';p.write_text(altered)
 for label,flags in [('O2',['-O2']),('O3-ndebug',['-O3','-DNDEBUG'])]:
  exe=D/('bad-refinement'+label);subprocess.run(['g++','-std=c++17',*flags,str(p),'-o',str(exe)],check=True,capture_output=True)
  r=subprocess.run([str(exe),str(B/'producer/parents-nine.txt'),str(B/'producer/trace-normal-per-parent.tsv'),str(D/'trace.tsv'),str(D/'trace.json')],capture_output=True,text=True)
  if r.returncode==0 or 'per-parent incoming histogram matches' not in r.stderr:raise RuntimeError('root corruption not caught')
  results.append({'mutation':'corrupt_root_size','mode':label,'exit_code':r.returncode,'failure':r.stderr.strip(),'changed_source_sha256':hashlib.sha256(altered.encode()).hexdigest()})
(B/'mutation-results.json').write_text(json.dumps({'status':'PASS','semantic_failures':len(results),'results':results},indent=2)+'\n')
print(json.dumps({'status':'PASS','semantic_failures':len(results)}))
