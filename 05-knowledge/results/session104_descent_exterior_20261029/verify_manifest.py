#!/usr/bin/env python3
"""Verify every delivered file against the packet's SHA-256 manifest."""
import hashlib,json
from pathlib import Path
root=Path(__file__).resolve().parent
manifest=json.loads((root/'MANIFEST.json').read_text())
for path,want in manifest.items():
 if hashlib.sha256((root/path).read_bytes()).hexdigest()!=want:
  raise RuntimeError('Manifest mismatch: '+path)
print(json.dumps({'status':'PASS','manifest_files':len(manifest)},sort_keys=True))
