#!/usr/bin/env python3
"""Verify published packet bytes, not mathematical truth."""
from pathlib import Path
import hashlib,json
root=Path(__file__).resolve().parent
manifest=json.loads((root/'MANIFEST.json').read_text())
for row in manifest['files']:
 p=root/row['path'];data=p.read_bytes()
 if len(data)!=row['bytes'] or hashlib.sha256(data).hexdigest()!=row['sha256']:
  raise SystemExit('Changed or missing: '+row['path'])
print('PASS:',len(manifest['files']),'packet files match their SHA-256 manifest.')
