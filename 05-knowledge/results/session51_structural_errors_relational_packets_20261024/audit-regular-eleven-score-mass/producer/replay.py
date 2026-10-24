#!/usr/bin/env python3
import hashlib,json,shutil,subprocess,sys,tempfile
from pathlib import Path
HERE=Path(__file__).resolve().parent
manifest=json.loads((HERE/'MANIFEST.json').read_text())
for item in manifest['files']:
 if hashlib.sha256((HERE/item['path']).read_bytes()).hexdigest()!=item['sha256']:raise ValueError('source manifest mismatch')
with tempfile.TemporaryDirectory(prefix='s51-score-mass-replay-') as tmp:
 packet=Path(tmp)/'packet';work=Path(tmp)/'unrelated';packet.mkdir();work.mkdir()
 for item in manifest['files']:
  p=packet/item['path'];p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(HERE/item['path'],p)
 shutil.copyfile(HERE/'MANIFEST.json',packet/'MANIFEST.json')
 for opt in (False,True):
  flags=['-O'] if opt else []
  run=subprocess.run([sys.executable]+flags+[str(packet/'check.py')],cwd=work,capture_output=True,check=True)
  if run.stdout!=(packet/'normal.json').read_bytes():raise ValueError('relocated receipt differs')
  subprocess.run([sys.executable]+flags+[str(packet/'verify_manifest.py')],cwd=work,capture_output=True,check=True)
 run=subprocess.run([sys.executable,str(packet/'hostiles.py')],cwd=work,capture_output=True,check=True)
 if run.stdout!=(packet/'hostile-results.json').read_bytes():raise ValueError('relocated hostile receipt differs')
print(json.dumps({'status':'PASS','manifest_only_physical_copy':True,'normal_optimized_identical':True,'unrelated_working_directory':True,'semantic_rejections':12,'manifest_rejections':12,'manifest_sha256':hashlib.sha256((HERE/'MANIFEST.json').read_bytes()).hexdigest(),'regular11_error':133210,'total11_error':152690,'unrestricted_order11_count':'not computed'},sort_keys=True,indent=2))
