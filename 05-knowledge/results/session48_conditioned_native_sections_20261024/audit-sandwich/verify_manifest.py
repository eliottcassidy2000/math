#!/usr/bin/env python3
"""Verify this audit's closed, portable, LF-normalized manifest."""
from hashlib import sha256
from pathlib import Path
import json

root=Path(__file__).resolve().parent
data=json.loads((root/'MANIFEST.json').read_text())
expected=data['files']
actual={str(p.relative_to(root)) for p in root.rglob('*')
        if p.is_file() and '__pycache__' not in p.parts and p!=root/'MANIFEST.json'}
if actual!=set(expected):
    raise RuntimeError('manifest closure: missing='+repr(sorted(set(expected)-actual))
                       +'; extra='+repr(sorted(actual-set(expected))))
for rel,want in expected.items():
    p=root/rel
    if root not in p.resolve().parents:
        raise RuntimeError('path escapes package: '+rel)
    got=sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
    if got!=want:
        raise RuntimeError('hash mismatch: '+rel)
print('PASS closed audit manifest:',len(expected),'files')
