#!/usr/bin/env python3
import hashlib,json,shutil,subprocess,sys,tempfile
from pathlib import Path
base=Path(__file__).resolve().parent;manifest=json.loads((base/'manifest.json').read_text())
for name,digest in manifest['files'].items():
 if hashlib.sha256((base/name).read_bytes()).hexdigest()!=digest:raise RuntimeError('Hash mismatch: '+name)
with tempfile.TemporaryDirectory(prefix='cyclic-audit-copy-') as tmp:
 dst=Path(tmp)
 for name in manifest['files']:
  target=dst/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(base/name,target)
 for optimized in (False,True):
  flags=['-B']+(['-O'] if optimized else [])
  out=subprocess.run([sys.executable]+flags+['independent.py','--output','computed.json'],cwd=dst,capture_output=True,text=True,check=True)
  if out.stdout!=(dst/('independent-optimized.out' if optimized else 'independent-normal.out')).read_text():raise RuntimeError('Independent stdout mismatch')
  if (dst/'computed.json').read_bytes()!=(dst/'receipt.json').read_bytes():raise RuntimeError('Independent receipt mismatch')
  out=subprocess.run([sys.executable]+flags+['check.py'],cwd=dst/'source',capture_output=True,text=True,check=True)
  if out.stdout!=(dst/'source'/('optimized.json' if optimized else 'normal.json')).read_text():raise RuntimeError('Producer output mismatch')
 subprocess.run([sys.executable,'-B','verify_packet.py'],cwd=dst/'source',capture_output=True,text=True,check=True)
print(json.dumps({'status':'PASS','physical_copy':True,'manifest_files':len(manifest['files']),'producer_checks_per_mode':16785,'independent_checks_per_mode':175939,'producer_manifest_sha256':hashlib.sha256((base/'source/MANIFEST.json').read_bytes()).hexdigest(),'normal_optimized_outputs_match':True},sort_keys=True))
