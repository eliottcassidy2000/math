#!/usr/bin/env python3
from hashlib import sha256
import json
from pathlib import Path
root=Path(__file__).resolve().parent
manifest=json.loads((root/'MANIFEST.json').read_text())
for row in manifest['files']:
    p=root/row['path']
    if not p.is_file() or sha256(p.read_bytes()).hexdigest()!=row['sha256']:
        raise RuntimeError('manifest mismatch: '+row['path'])
print('PASS:',len(manifest['files']),'manifest hashes')
