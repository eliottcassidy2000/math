#!/usr/bin/env python3
"""Check portable bundle hashes relative to this script, without a checkout."""
from pathlib import Path
from hashlib import sha256
import json

root = Path(__file__).resolve().parent
manifest = json.loads((root/'PORTABLE_MANIFEST.json').read_text())
for item in manifest['files']:
    path = root/item['path']
    data = path.read_bytes()
    if len(data) != item['bytes'] or sha256(data).hexdigest() != item['sha256']:
        raise RuntimeError('manifest mismatch: '+item['path'])
print('PASS: portable manifest, '+str(len(manifest['files']))+' files')
