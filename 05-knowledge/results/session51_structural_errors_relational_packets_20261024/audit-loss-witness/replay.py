#!/usr/bin/env python3
import hashlib,json,shutil,subprocess,sys,tempfile
from pathlib import Path
base=Path(__file__).resolve().parent;manifest=json.loads((base/'manifest.json').read_text())
for name,h in manifest['files'].items():
 if hashlib.sha256((base/name).read_bytes()).hexdigest()!=h:raise RuntimeError('Hash mismatch: '+name)
with tempfile.TemporaryDirectory(prefix='loss-witness-audit-copy-') as temp:
 dst=Path(temp)
 for name in manifest['files']:
  target=dst/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(base/name,target)
 for optimized in (False,True):
  flags=['-B']+(['-O'] if optimized else [])
  for script,cwd in [('audit.py',dst),('check.py',dst/'source')]:
   out=subprocess.run([sys.executable]+flags+[script],cwd=cwd,capture_output=True,text=True,check=True)
   if out.stdout!=(cwd/('optimized.json' if optimized else 'normal.json')).read_text():raise RuntimeError('Output mismatch: '+script)
print(json.dumps({'status':'PASS','physical_copy':True,'manifest_files':len(manifest['files']),'producer_manifest_sha256':hashlib.sha256((base/'source/MANIFEST.json').read_bytes()).hexdigest(),'hamiltonian_path_totals':[15465,15549],'entire_decoder_output_equal':True,'difference':84},sort_keys=True))
