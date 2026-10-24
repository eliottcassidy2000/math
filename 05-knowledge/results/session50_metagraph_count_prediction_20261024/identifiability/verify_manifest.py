#!/usr/bin/env python3
from pathlib import Path
import hashlib,json
root=Path(__file__).resolve().parent
manifest=json.loads((root/'MANIFEST.json').read_text())
for row in manifest['files']:
    path=root/row['path']
    if not path.is_file():raise RuntimeError('missing '+row['path'])
    data=path.read_bytes()
    if len(data)!=row['bytes'] or hashlib.sha256(data).hexdigest()!=row['sha256']:
        raise RuntimeError('hash mismatch '+row['path'])
for row in json.loads((root/'SOURCE_PINS.json').read_text()):
    data=(root/row['archive_path']).read_bytes()
    if len(data)!=row['bytes'] or hashlib.sha256(data).hexdigest()!=row['sha256']:
        raise RuntimeError('source archive mismatch')
print(json.dumps({'manifest':'PASS','files':len(manifest['files'])},sort_keys=True))
