#!/usr/bin/env python3
"""Check the packet's declared relative files without using assertions."""
from pathlib import Path
import hashlib
import json
HERE=Path(__file__).resolve().parent
manifest=json.loads((HERE/'MANIFEST.json').read_text())
seen=set()
for record in manifest['files']:
    relative=Path(record['path'])
    if relative.is_absolute() or '..' in relative.parts or record['path'] in seen:
        raise RuntimeError('unsafe or duplicate manifest path: '+record['path'])
    seen.add(record['path'])
    path=HERE/relative
    if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=record['sha256']:
        raise RuntimeError('manifest mismatch: '+record['path'])
print(json.dumps({'status':'PASS','files_verified':len(manifest['files'])},sort_keys=True))
