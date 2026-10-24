#!/usr/bin/env python3
"""Rebuild and rerun the packet using ONLY manifest-listed artifacts.
The nested original input manifest is data and must be delivered too.
"""
import hashlib,json,shutil,subprocess,sys,tempfile
from pathlib import Path
H=Path(__file__).resolve().parent
M=json.loads((H/'MANIFEST.json').read_text());listed={r['path'] for r in M['files']}
pins=json.loads((H/'SOURCE_PINS.json').read_text())['files']
for r in pins:
 if r['snapshot'] not in listed:raise RuntimeError('Undelivered source snapshot '+r['snapshot'])
pairs=[('','check.py'),('producer-native.','native570'),('producer-audit.','audit570'),('guard-completeness.','replay_guard_completeness.py'),('odd-part.','check_odd_part.py'),('helper-raw.','check_helper_raw.py')]
with tempfile.TemporaryDirectory(prefix='fresh570-manifest-only-') as temp:
 R=Path(temp)/'renamed-audit-package';R.mkdir()
 for r in M['files']:
  rel=Path(r['path'])
  if rel.is_absolute() or '..' in rel.parts:raise RuntimeError('Unsafe manifest path')
  data=(H/rel).read_bytes()
  if len(data)!=r['bytes'] or hashlib.sha256(data).hexdigest()!=r['sha256']:raise RuntimeError('Source manifest hash mismatch '+str(rel))
  dst=R/rel;dst.parent.mkdir(parents=True,exist_ok=True);dst.write_bytes(data)
 shutil.copyfile(H/'MANIFEST.json',R/'MANIFEST.json')
 def run(args):
  proc=subprocess.run([sys.executable,*args],cwd=R,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
  if proc.returncode:raise RuntimeError('Isolated copied-package replay failed: '+repr(args)+'\n'+proc.stdout+'\n'+proc.stderr)
  return proc.stdout
 run(['-I','-B','verify_manifest.py'])
 # Replay mathematical modules with isolated Python: no cwd/PYTHONPATH fallback.
 for stem,program in pairs:
  for opt,mode in [(False,'normal'),(True,'optimized')]:
   flags=['-I']+(['-O'] if opt else [])+['-B']
   dest=stem+mode+'.json'
   if stem.startswith('producer-'):
    text=run(flags+['replay_native.py',program]);(R/dest).write_text(text)
   else:run(flags+[program,'--out',dest])
   if (R/dest).read_bytes()!=(H/dest).read_bytes():raise RuntimeError('Changed receipt in isolated manifest-only package '+dest)
 run(['-I','-B','verify_manifest.py']);run(['-I','-O','-B','verify_manifest.py'])
result={'status':'PASS: standalone manifest-only renamed-package replay','source_pin_snapshots_delivered':len(pins),'normal_optimized_pairs':6,'python_isolation':'Every child uses -I and cwd=the temporary renamed copy; no original workspace import fallback. Only manifest-listed files plus the root manifest are copied.','receipt_sha256':{stem+mode+'.json':hashlib.sha256((H/(stem+mode+'.json')).read_bytes()).hexdigest() for stem,program in pairs for mode in ['normal','optimized']},'nested_original_manifest_delivered':'inputs/native/MANIFEST.json','scope':'Packaging and path closure only; mathematical proofs and receipt bytes unchanged.'}
print(json.dumps(result,indent=2,sort_keys=True))
