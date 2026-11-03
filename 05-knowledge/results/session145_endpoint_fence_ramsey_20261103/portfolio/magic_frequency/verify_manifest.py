#!/usr/bin/env python3
from pathlib import Path
import hashlib,json
root=Path(__file__).resolve().parent
manifest=json.loads((root/'MANIFEST.json').read_text())
for rel,record in manifest.items():
    data=(root/rel).read_bytes()
    if len(data)!=record['bytes'] or hashlib.sha256(data).hexdigest()!=record['sha256']:
        raise ValueError('manifest mismatch: '+rel)
print(json.dumps({'status':'PASS','verified_files':len(manifest)},sort_keys=True))
