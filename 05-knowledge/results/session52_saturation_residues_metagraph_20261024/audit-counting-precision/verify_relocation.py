#!/usr/bin/env python3
"""Replay only declared audit dependencies from an actual isolated copy."""
import hashlib,json,os,shutil,subprocess,sys,tempfile
from pathlib import Path
HERE=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def require(x,m):
 if not x:raise RuntimeError(m)
manifest=json.loads((HERE/'relocation-manifest.json').read_text())
for e in manifest['files']:
 p=Path(e['path']);require(not p.is_absolute() and '..' not in p.parts,'unsafe path')
 require(not (HERE/p).is_symlink() and sha(HERE/p)==e['sha256'],'input digest mismatch')
with tempfile.TemporaryDirectory(prefix='precision-manifest-only-') as td:
 root=Path(td)/'packet';root.mkdir()
 for e in manifest['files']:
  p=Path(e['path']);(root/p).parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(HERE/p,root/p)
 require(sorted(p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file())==sorted(e['path'] for e in manifest['files']),'file closure')
 results=[]
 env=dict(os.environ);env.pop('PYTHONPATH',None)
 for opt in [False,True]:
  for kind in ['producer','independent']:
   receipt=Path(td)/(kind+str(opt)+'.json')
   script=root/('inputs/producer/verify_precision.py' if kind=='producer' else 'audit_precision.py')
   args=['--repo',str(root/'inputs/sources')] if kind=='producer' else []
   r=subprocess.run([sys.executable,'-B','-I']+(['-O'] if opt else [])+[str(script)]+args+['--receipt',str(receipt)],cwd='/tmp',env=env,capture_output=True,text=True)
   require(r.returncode==0,r.stderr)
   expected=HERE/(kind+('-optimized.json' if opt else '-normal.json'))
   require(receipt.read_bytes()==expected.read_bytes(),'receipt mismatch')
   results.append({'checker':kind,'optimized':opt,'checks':json.loads(receipt.read_text())['checks'],'receipt_sha256':sha(receipt)})
 # Bad pinned source must fail before any receipt, in both execution modes.
 target=root/'inputs/sources'/json.loads((root/'inputs/producer/source-pins.json').read_text())['sources'][0]['path']
 target.write_bytes(target.read_bytes()+b'\n')
 hostiles=[]
 for opt in [False,True]:
  receipt=Path(td)/('hostile'+str(opt)+'.json')
  r=subprocess.run([sys.executable,'-B','-I']+(['-O'] if opt else [])+[str(root/'inputs/producer/verify_precision.py'),'--repo',str(root/'inputs/sources'),'--receipt',str(receipt)],cwd='/tmp',env=env,capture_output=True,text=True)
  require(r.returncode!=0 and 'changed source' in r.stderr and not receipt.exists(),'source hostile accepted')
  hostiles.append({'optimized':opt,'rejected':True,'no_success_receipt':True})
for e in manifest['files']:require(sha(HERE/e['path'])==e['sha256'],'audit input mutated')
print(json.dumps({'status':'PASS','actual_manifest_only_copy':True,'files_copied':len(manifest['files']),'cwd':'/tmp','isolated_python':True,'runtime_dependencies':'Python standard library and system patch executable; only two pinned theorem inputs; no repo/network imports','runs':results,'changed_source_hostiles':hostiles},indent=2))
