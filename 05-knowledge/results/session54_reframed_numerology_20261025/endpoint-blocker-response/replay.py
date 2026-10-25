#!/usr/bin/env python3
import hashlib,json,pathlib,subprocess,sys,tempfile,shutil
D=pathlib.Path(__file__).resolve().parent

def ck(x,msg):
 if not x:raise RuntimeError(msg)
def check(root):
 manifest=json.loads((root/'MANIFEST.json').read_text())
 actual={str(p.relative_to(root)) for p in root.rglob('*') if p.is_file() and p.name!='MANIFEST.json'}
 ck(actual==set(manifest),'exact file set mismatch')
 for name,sha in manifest.items():ck(hashlib.sha256((root/name).read_bytes()).hexdigest()==sha,'digest mismatch: '+name)
 for source in json.loads((root/'SOURCES.json').read_text()):
  ck(hashlib.sha256((root/source['bundled']).read_bytes()).hexdigest()==source['sha256'],'source pin mismatch')
 return manifest
check(D)
for flags,filename in (([],'result.json'),(['-O'],'result-O.json')):
 p=subprocess.run([sys.executable]+flags+[str(D/'verify.py')],cwd='/tmp',capture_output=True,check=True)
 ck(p.stdout==(D/filename).read_bytes(),'replay mismatch: '+filename)
ck((D/'result.json').read_bytes()==(D/'result-O.json').read_bytes(),'normal optimized mismatch')
if '--self-test' in sys.argv:
 for kind in ('content','extra','missing'):
  with tempfile.TemporaryDirectory() as tmp:
   root=pathlib.Path(tmp)/'packet';shutil.copytree(D,root)
   if kind=='content':(root/'PROOF.txt').write_text((root/'PROOF.txt').read_text()+'\nmutation\n')
   elif kind=='extra':(root/'extra.txt').write_text('mutation')
   else:(root/'inputs/receipt.json').unlink()
   try:check(root)
   except RuntimeError:pass
   else:raise RuntimeError('mutation accepted: '+kind)
print('PASS: source pins, exact manifest, normal/-O replay'+', content/extra/missing hostiles' * ('--self-test' in sys.argv))
