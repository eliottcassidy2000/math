#!/usr/bin/env python3
import hashlib,json,shutil,subprocess,sys,tempfile
from pathlib import Path
HERE=Path(__file__).resolve().parent
m=json.loads((HERE/'MANIFEST.json').read_text())
for r in m['files']:
 if hashlib.sha256((HERE/r['path']).read_bytes()).hexdigest()!=r['sha256']:raise ValueError('source pin mismatch: '+r['path'])
with tempfile.TemporaryDirectory(prefix='s51-nonrigid-relocation-') as tmp:
 packet=Path(tmp)/'packet';work=Path(tmp)/'unrelated';packet.mkdir();work.mkdir()
 for r in m['files']:
  p=packet/r['path'];p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(HERE/r['path'],p)
 shutil.copyfile(HERE/'MANIFEST.json',packet/'MANIFEST.json')
 for opt in (False,True):
  flags=['-O'] if opt else []
  run=subprocess.run([sys.executable]+flags+[str(packet/'check.py')],cwd=work,capture_output=True,check=True)
  if run.stdout!=(packet/'normal.json').read_bytes():raise ValueError('relocated receipt differs')
  subprocess.run([sys.executable]+flags+[str(packet/'verify_manifest.py')],cwd=work,capture_output=True,check=True)
 run=subprocess.run([sys.executable,str(packet/'hostiles.py')],cwd=work,capture_output=True,check=True)
 if run.stdout!=(packet/'hostile-results.json').read_bytes():raise ValueError('relocated hostile receipt differs')
print(json.dumps({'status':'PASS','manifest_only_physical_copy':True,'unrelated_working_directory':True,'normal_optimized_identical':True,'semantic_rejections':12,'manifest_rejections':12,'nonrigid_regular11_classes':18,'pointed_nonrigid_states':78,'nonrigid_error':660,'nonrigid_unpointing_deficit':120,'manifest_sha256':hashlib.sha256((HERE/'MANIFEST.json').read_bytes()).hexdigest()},sort_keys=True,indent=2))
