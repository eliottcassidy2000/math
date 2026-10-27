#!/usr/bin/env python3
"""Create or verify an exact-byte manifest. Use --verify for read-only replay."""
from pathlib import Path
import hashlib,json,sys
ROOT=Path(__file__).resolve().parent
EXCLUDE={'MANIFEST.json','SHA256SUMS'}
files=sorted(p for p in ROOT.rglob('*') if p.is_file() and p.name not in EXCLUDE and '__pycache__' not in p.parts)
rows=[{'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in files]
obj={'status':'SEALED','hash_semantics':'native exact bytes, no newline normalization','files':rows}
encoded=json.dumps(obj,indent=2,sort_keys=True)+'\n'
sums=''.join(r['sha256']+'  '+r['path']+'\n' for r in rows)
if '--verify' in sys.argv:
    if (ROOT/'MANIFEST.json').read_text()!=encoded or (ROOT/'SHA256SUMS').read_text()!=sums:
        raise RuntimeError('manifest mismatch')
    for row in json.loads((ROOT/'SOURCE_PINS.json').read_text()):
        p=ROOT/row['snapshot']
        if hashlib.sha256(p.read_bytes()).hexdigest()!=row['sha256']:
            raise RuntimeError('source pin mismatch: '+row['snapshot'])
    if (ROOT/'receipt.json').read_bytes()!=(ROOT/'receipt.optimized.json').read_bytes():
        raise RuntimeError('normal/-O receipt mismatch')
    print(json.dumps({'status':'PASS','sealed_files':len(rows),'normal_optimized':'byte-identical'},sort_keys=True))
else:
    (ROOT/'MANIFEST.json').write_text(encoded)
    (ROOT/'SHA256SUMS').write_text(sums)
    print(json.dumps({'status':'SEALED','files':len(rows)},sort_keys=True))
