#!/usr/bin/env python3
"""Verify the complete audit payload inventory and SHA-256 pins."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
expected = json.loads((ROOT / 'MANIFEST.json').read_text())['files']
actual = {p.relative_to(ROOT).as_posix(): p for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p != ROOT / 'MANIFEST.json'}
if set(expected) != set(actual):
    raise RuntimeError('audit inventory mismatch')
for name, path in actual.items():
    if hashlib.sha256(path.read_bytes()).hexdigest() != expected[name]:
        raise RuntimeError('audit hash mismatch: ' + name)
print('PASS: ' + str(len(actual)) + ' audit payload files')
