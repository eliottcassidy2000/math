#!/usr/bin/env python3
from pathlib import Path
import json,hashlib,subprocess,sys,tempfile,shutil
D=Path(__file__).resolve().parent
def require(b,m):
 if not b:raise RuntimeError(m)
def validate(root):
 manifest=json.loads((root/'MANIFEST.json').read_text())
 files={str(p.relative_to(root)) for p in root.rglob('*') if p.is_file() and p!=root/'MANIFEST.json'}
 require(files==set(manifest),'exact file set')
 for name,sha in manifest.items():require(hashlib.sha256((root/name).read_bytes()).hexdigest()==sha,'hash '+name)
validate(D)
for flags,result in (([],'result.json'),(['-O'],'result-O.json')):
 r=subprocess.run([sys.executable,*flags,str(D/'audit.py')],cwd='/tmp',capture_output=True,check=True)
 require(r.stdout==(D/result).read_bytes(),'fresh replay '+result)
require((D/'result.json').read_bytes()==(D/'result-O.json').read_bytes(),'optimization agreement')
r=subprocess.run([sys.executable,str(D/'inputs/producer/replay.py'),'--self-test'],cwd='/tmp',capture_output=True,check=True)
require(r.stdout==(D/'producer-replay.txt').read_bytes(),'producer replay')
if '--self-test' in sys.argv:
 for kind in ('change','add','remove'):
  with tempfile.TemporaryDirectory() as t:
   root=Path(t)/'packet';shutil.copytree(D,root)
   if kind=='change':(root/'AUDIT.txt').write_text('corruption')
   elif kind=='add':(root/'extra.txt').write_text('unexpected')
   else:(root/'inputs/producer/PROOF.txt').unlink()
   try:validate(root)
   except RuntimeError:pass
   else:raise RuntimeError('accepted '+kind)
print('PASS: closed independent endpoint audit; normal/-O; producer replay; content/extra/missing hostiles')
