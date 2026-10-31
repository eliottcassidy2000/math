#!/usr/bin/env python3
"""Verify the explicit release whitelist and SHA-256 manifest."""
from pathlib import Path
import json,hashlib
ROOT=Path(__file__).resolve().parent
release=json.loads((ROOT/'RELEASE_FILES.json').read_text())
manifest=json.loads((ROOT/'MANIFEST.json').read_text())['sha256']
if set(manifest)!=set(release)-{'MANIFEST.json'}:raise RuntimeError('manifest/whitelist discrepancy')
for name,want in manifest.items():
 p=ROOT/name
 if not p.is_file():raise RuntimeError('missing release file '+name)
 got=hashlib.sha256(p.read_bytes()).hexdigest()
 if got!=want:raise RuntimeError('hash mismatch '+name)
print(json.dumps({'status':'PASS','checked_files':len(manifest)},sort_keys=True))
