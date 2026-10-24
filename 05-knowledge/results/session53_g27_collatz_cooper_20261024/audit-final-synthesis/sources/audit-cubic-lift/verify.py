import hashlib,json,subprocess,sys
from pathlib import Path
root=Path(__file__).resolve().parent

def unique(pairs):
 d={}
 for k,v in pairs:
  if k in d:raise RuntimeError('duplicate key')
  d[k]=v
 return d
m=json.loads((root/'MANIFEST.json').read_text(),object_pairs_hook=unique)
expected=set(m['files'])|{'MANIFEST.json'}
actual={p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file() or p.is_symlink()}
if expected!=actual:raise RuntimeError('membership mismatch')
for rel,h in m['files'].items():
 p=root/rel
 if Path(rel).is_absolute() or '..' in Path(rel).parts or p.is_symlink() or not p.resolve().is_relative_to(root):raise RuntimeError('unsafe member')
 if hashlib.sha256(p.read_bytes()).hexdigest()!=h:raise RuntimeError('hash mismatch')
for flag,fixture in [([], 'normal.json'),(['-O'],'optimized.json')]:
 run=subprocess.run([sys.executable,*flag,str(root/'check.py')],capture_output=True,check=True)
 if run.stdout!=(root/fixture).read_bytes():raise RuntimeError('fixture mismatch')
print('PASS: independent exact checker in both modes, closed hashes and membership')
