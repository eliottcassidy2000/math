#!/usr/bin/env python3
from pathlib import Path
import hashlib,json
HERE=Path(__file__).resolve().parent
manifest=json.loads((HERE/'MANIFEST.json').read_text())
for record in manifest['files']:
    path=HERE/record['path']
    if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=record['sha256']:
        raise RuntimeError('manifest mismatch: '+record['path'])
print(json.dumps({'status':'PASS','files_verified':len(manifest['files'])},sort_keys=True))
