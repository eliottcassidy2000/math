#!/usr/bin/env python3
"""Exact byte/file-set check; all gates remain active under python -O."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parent
symlinks=[p for p in ROOT.rglob('*') if p.is_symlink()]
if symlinks:raise RuntimeError('Symlinks: '+str(symlinks))
manifest=json.loads((ROOT/'MANIFEST.json').read_text())
actual={p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p!=ROOT/'MANIFEST.json'}
expected={row['path'] for row in manifest['files']}
if len(expected)!=len(manifest['files']):raise RuntimeError('Duplicate manifest path')
if actual!=expected:raise RuntimeError(f'File set mismatch: missing={expected-actual}, extra={actual-expected}')
for row in manifest['files']:
 p=ROOT/row['path']
 if p.is_symlink():raise RuntimeError(f'Symlink: {p}')
 b=p.read_bytes()
 if len(b)!=row['bytes'] or hashlib.sha256(b).hexdigest()!=row['sha256']:raise RuntimeError(f'Byte mismatch: {p}')
print(f"PASS: {len(expected)} exact files, {sum(row['bytes'] for row in manifest['files'])} bytes")
