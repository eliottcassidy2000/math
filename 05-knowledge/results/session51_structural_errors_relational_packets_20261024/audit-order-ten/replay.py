#!/usr/bin/env python3
"""Manifest-only physical-copy replay, not execution in producer directory."""
import hashlib,json,shutil,subprocess,sys,tempfile
from pathlib import Path
base=Path(__file__).resolve().parent;manifest=json.loads((base/'manifest.json').read_text())
for name,digest in manifest['files'].items():
 if hashlib.sha256((base/name).read_bytes()).hexdigest()!=digest:raise RuntimeError('Hash mismatch: '+name)
with tempfile.TemporaryDirectory(prefix='order10-independent-copy-') as tmp:
 dst=Path(tmp)
 for name in manifest['files']:
  target=dst/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(base/name,target)
 for args,outname in [([], 'normal.out'),(['-O'],'optimized.out')]:
  result=subprocess.run([sys.executable]+args+['audit.py','--output','computed.json'],cwd=dst,capture_output=True,text=True,check=True)
  if result.stdout!=(dst/outname).read_text():raise RuntimeError('stdout mismatch: '+outname)
  if (dst/'computed.json').read_bytes()!=(dst/'receipt.json').read_bytes():raise RuntimeError('receipt mismatch: '+outname)
 for optimized in (False,True):
  flags=['-O'] if optimized else []
  tests=[('distribution_audit.py',dst,'distribution-optimized.json' if optimized else 'distribution-normal.json'),('check.py',dst/'producer','optimized.json' if optimized else 'normal.json'),('distribution_check.py',dst/'producer','distribution-optimized.json' if optimized else 'distribution-normal.json')]
  for script,cwd,expected in tests:
   result=subprocess.run([sys.executable]+flags+[script],cwd=cwd,capture_output=True,text=True,check=True)
   if result.stdout!=(cwd/expected).read_text():raise RuntimeError('replay mismatch: '+script)
  subprocess.run([sys.executable]+flags+['verify_manifest.py'],cwd=dst/'producer',capture_output=True,text=True,check=True)
print(json.dumps({'physical_copy':True,'manifest_files':len(manifest['files']),'normal_optimized_receipts_match':True,'producer_manifest_sha256':hashlib.sha256((base/'producer/MANIFEST.json').read_bytes()).hexdigest(),'strong_error_10':3992,'full_error_10':5508,'errorful_classes_10':160,'raw_second_moment_10':229416},sort_keys=True))
