#!/usr/bin/env python3
from pathlib import Path
import hashlib,json
root=Path(__file__).resolve().parent
items=json.loads((root/'MANIFEST.json').read_text())
for p,want in items.items():
 if hashlib.sha256((root/p).read_bytes()).hexdigest()!=want:raise RuntimeError('hash mismatch: '+p)
print(json.dumps({'status':'PASS','files':len(items)},sort_keys=True))
