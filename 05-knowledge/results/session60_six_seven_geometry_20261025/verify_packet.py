#!/usr/bin/env python3
"""Verify the exact published packet; standard library only, no removed asserts."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parent
manifest=json.loads((ROOT/'MANIFEST.json').read_text())
actual={p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name!='MANIFEST.json'}
expected={r['path'] for r in manifest['files']}
if actual!=expected:raise RuntimeError(f'File set mismatch: missing={expected-actual}, extra={actual-expected}')
for row in manifest['files']:
 p=ROOT/row['path']
 if p.is_symlink():raise RuntimeError(f'Symlink: {p}')
 b=p.read_bytes()
 if len(b)!=row['bytes'] or hashlib.sha256(b).hexdigest()!=row['sha256']:raise RuntimeError(f'Byte mismatch: {p}')
print(f"PASS: {len(expected)} exact files, {sum(r['bytes'] for r in manifest['files'])} bytes")
