import hashlib,json,subprocess,sys
from pathlib import Path
p=Path(__file__).resolve().parent
manifest=json.loads((p/'MANIFEST.json').read_text())
expected=set(manifest['files'])|{'MANIFEST.json'}
actual={f.relative_to(p).as_posix() for f in p.rglob('*') if f.is_file()}
if actual!=expected:raise RuntimeError('closed file set mismatch')
for name,h in manifest['files'].items():
    if hashlib.sha256((p/name).read_bytes()).hexdigest()!=h:raise RuntimeError('hash mismatch '+name)
for source,normal,optimized in [('check.py', 'check.normal.json', 'check.optimized.json'), ('check_modular.py', 'modular.normal.json', 'modular.optimized.json')]:
    a=subprocess.check_output([sys.executable,str(p/source)],cwd=p)
    b=subprocess.check_output([sys.executable,'-O',str(p/source)],cwd=p)
    if a!=b or a!=(p/normal).read_bytes() or b!=(p/optimized).read_bytes():raise RuntimeError('replay mismatch '+source)
print('PASS: manifest and normal/-O replay')
