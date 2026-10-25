#!/usr/bin/env python3
"""Reject missing, altered, or extra payload files, including nested extras."""
import hashlib,json
from pathlib import Path
root=Path(__file__).resolve().parent
manifest=json.loads((root/'MANIFEST.json').read_text())
expected=set(manifest['files'])
# No cache whitelist is needed: executable scripts do not import local modules.
actual={str(p.relative_to(root)) for p in root.rglob('*') if p.is_file() and p!=root/'MANIFEST.json'}
if actual!=expected:
 raise RuntimeError('file set mismatch: missing='+repr(sorted(expected-actual))+' extra='+repr(sorted(actual-expected)))
for name,digest in manifest['files'].items():
 path=root/name
 if path.is_symlink() or not path.resolve().is_relative_to(root):raise RuntimeError('invalid payload path '+name)
 if hashlib.sha256(path.read_bytes()).hexdigest()!=digest:raise RuntimeError('hash mismatch '+name)
for snapshot,record in manifest.get('inherited_sources',{}).items():
 if snapshot not in expected:raise RuntimeError('unlisted inherited source '+snapshot)
 if manifest['files'][snapshot]!=record['sha256']:raise RuntimeError('inherited hash mismatch '+snapshot)
print('Verified exact file set and '+str(len(expected))+' payload hashes.')
