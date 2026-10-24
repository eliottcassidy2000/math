#!/usr/bin/env python3
from pathlib import Path
from hashlib import sha256
import json
H=Path(__file__).resolve().parent
manifest=json.loads((H/'MANIFEST.json').read_text())
for r in manifest['files']:
    if sha256((H/r['path']).read_bytes()).hexdigest()!=r['sha256']:raise RuntimeError('manifest mismatch: '+r['path'])
print(json.dumps({'status':'PASS','files_verified':len(manifest['files'])},sort_keys=True))
