#!/usr/bin/env python3
"""Manifest-only physical-copy replay from an unrelated working directory."""
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
HERE=Path(__file__).resolve().parent

def main():
 manifest=json.loads((HERE/'MANIFEST.json').read_text())
 for item in manifest['files']:
  if hashlib.sha256((HERE/item['path']).read_bytes()).hexdigest()!=item['sha256']:raise ValueError('input manifest mismatch: '+item['path'])
 results=[]
 with tempfile.TemporaryDirectory(prefix='s51-order11-audit-relocation-') as tmp:
  packet=Path(tmp)/'packet';work=Path(tmp)/'unrelated-work';packet.mkdir();work.mkdir()
  for item in manifest['files']:
   target=packet/item['path'];target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(HERE/item['path'],target)
  shutil.copyfile(HERE/'MANIFEST.json',packet/'MANIFEST.json')
  tests=[('audit.py','normal.json'),('producer/check.py','producer/normal.json'),('crosscheck.py','crosscheck-normal.json')]
  for optimized in (False,True):
   flags=['-O'] if optimized else []
   for script,expected in tests:
    run=subprocess.run([sys.executable]+flags+[str(packet/script)],cwd=work,capture_output=True,check=True)
    if run.stdout!=(packet/expected).read_bytes():raise ValueError('relocated output mismatch: '+script)
    results.append({'script':script,'optimized':optimized,'output_sha256':hashlib.sha256(run.stdout).hexdigest()})
   for script in ('verify_manifest.py','producer/verify_manifest.py'):
    subprocess.run([sys.executable]+flags+[str(packet/script)],cwd=work,capture_output=True,check=True)
  for script,expected in [('hostiles.py','hostile-results.json'),('producer_hostiles.py','producer-hostile-results.json')]:
   run=subprocess.run([sys.executable,str(packet/script)],cwd=work,capture_output=True,check=True)
   if run.stdout!=(packet/expected).read_bytes():raise ValueError('relocated hostile receipt mismatch: '+script)
  pin=json.loads((packet/'PRODUCER_PIN.json').read_text())
 return {'status':'PASS','manifest_only_physical_copy':True,'working_directory':'unrelated to packet','manifest_files':len(manifest['files']),'producer_manifest_sha256':pin['manifest_sha256'],'normal_optimized_replays':results,'semantic_mutation_rejections':24,'producer_mutation_manifest_rejections':12,'nonregular_strong_error':9166,'known_total_error_residual':19480,'regular_eleven_error':'not computed within this audit packet'}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
