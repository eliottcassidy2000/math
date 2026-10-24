import hashlib,json,subprocess,sys
from pathlib import Path
root=Path(__file__).resolve().parent
manifest=json.loads((root/'manifest.json').read_text())
expected=set(manifest['files'])|{'manifest.json'}
actual={p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file()}
if actual!=expected:
    raise RuntimeError(f'membership mismatch: missing={sorted(expected-actual)}, extra={sorted(actual-expected)}')
for rel,want in manifest['files'].items():
    p=root/rel
    if p.is_symlink() or not p.resolve().is_relative_to(root):
        raise RuntimeError('unsafe member '+rel)
    if hashlib.sha256(p.read_bytes()).hexdigest()!=want:
        raise RuntimeError('hash mismatch '+rel)
for flag,fixture in [([], 'normal.json'),(['-O'],'optimized.json')]:
    run=subprocess.run([sys.executable,*flag,str(root/'verify.py')],capture_output=True,check=True)
    if run.stdout!=(root/fixture).read_bytes(): raise RuntimeError('replay fixture mismatch '+fixture)
print('PASS: closed membership, hashes, normal and optimized replay')
