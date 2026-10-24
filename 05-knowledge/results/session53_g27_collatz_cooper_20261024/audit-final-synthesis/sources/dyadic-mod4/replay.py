import hashlib,json,subprocess,sys
from pathlib import Path
p=Path(__file__).resolve().parent
m=json.loads((p/'MANIFEST.json').read_text())
actual={str(f.relative_to(p)) for f in p.rglob('*') if f.is_file()}
if actual!=set(m['files'])|{'MANIFEST.json'}:raise RuntimeError('closed recursive membership mismatch')
for name,h in m['files'].items():
    if hashlib.sha256((p/name).read_bytes()).hexdigest()!=h:raise RuntimeError('hash mismatch: '+name)
a=subprocess.check_output([sys.executable,str(p/'check.py')],cwd=p)
b=subprocess.check_output([sys.executable,'-O',str(p/'check.py')],cwd=p)
if a!=b or a!=(p/'normal.json').read_bytes() or b!=(p/'optimized.json').read_bytes():raise RuntimeError('independent replay mismatch')
print('PASS: independent normal/-O and recursive closed manifest')
