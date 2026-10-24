#!/usr/bin/env python3
"""Verify an exact closed path inventory and all bytes; no assert statements."""
import hashlib,json
from pathlib import Path
root=Path(__file__).resolve().parent
manifest=json.loads((root/'MANIFEST.json').read_text())
actual={str(p.relative_to(root)) for p in root.rglob('*') if p.is_file() and p!=root/'MANIFEST.json'}
expected=set(manifest['files'])
if actual!=expected:raise RuntimeError({'unexpected':sorted(actual-expected),'missing':sorted(expected-actual)})
for name,record in manifest['files'].items():
    p=root/name
    if p.is_symlink():raise RuntimeError('Symlink forbidden: '+name)
    raw=p.read_bytes()
    if len(raw)!=record['bytes'] or hashlib.sha256(raw).hexdigest()!=record['sha256']:raise RuntimeError('Mismatched bytes: '+name)
print(json.dumps({'status':'PASS','closed_manifest_files':len(actual)},sort_keys=True))
