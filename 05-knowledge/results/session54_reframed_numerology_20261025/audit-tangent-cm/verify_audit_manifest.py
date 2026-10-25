import json,hashlib
from pathlib import Path
r=Path(__file__).resolve().parent
m=json.loads((r/'MANIFEST.json').read_text())
a={str(p.relative_to(r)) for p in r.rglob('*') if p.is_file() and p.name!='MANIFEST.json'}
if a!=set(m['files']):raise RuntimeError('audit file set mismatch')
for p,h in m['files'].items():
 if hashlib.sha256((r/p).read_bytes()).hexdigest()!=h:raise RuntimeError('audit hash mismatch '+p)
print('Audit file set and hashes verified.')
