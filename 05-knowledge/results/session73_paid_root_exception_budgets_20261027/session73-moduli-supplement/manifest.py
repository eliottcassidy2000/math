#!/usr/bin/env python3
from pathlib import Path
import hashlib,json,sys
root=Path(__file__).resolve().parent
files=sorted(p for p in root.rglob('*') if p.is_file() and p.name not in {'MANIFEST.json','SHA256SUMS'} and '__pycache__' not in p.parts)
rows=[{'path':str(p.relative_to(root)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in files]
obj={'hash_semantics':'native exact bytes','status':'SEALED','files':rows}
s=json.dumps(obj,indent=2,sort_keys=True)+'\n'
h=''.join(r['sha256']+'  '+r['path']+'\n' for r in rows)
if '--verify' in sys.argv:
    if (root/'MANIFEST.json').read_text()!=s or (root/'SHA256SUMS').read_text()!=h:
        raise RuntimeError('supplement manifest mismatch')
    for a,b in [('normal.json','optimized.json'),('patterns-normal.json','patterns-optimized.json')]:
        if (root/a).read_bytes()!=(root/b).read_bytes():
            raise RuntimeError('receipt mismatch: '+a)
    for p in json.loads((root/'SOURCE_PINS.json').read_text()):
        if hashlib.sha256((root/p['snapshot']).read_bytes()).hexdigest()!=p['sha256']:
            raise RuntimeError('source mismatch: '+p['snapshot'])
    print(json.dumps({'status':'PASS','sealed_files':len(rows),'receipt_pairs':'byte-identical'},sort_keys=True))
else:
    (root/'MANIFEST.json').write_text(s)
    (root/'SHA256SUMS').write_text(h)
    print(json.dumps({'status':'SEALED','files':len(rows)},sort_keys=True))
