#!/usr/bin/env python3
"""Verify packet exact file membership, bytes, source pins, and replay equality."""
from pathlib import Path
import hashlib
import json

root = Path(__file__).resolve().parent
manifest = json.loads((root / 'MANIFEST.json').read_text())
expected = {r['path']: r for r in manifest['files']}
actual = {str(p.relative_to(root)) for p in root.rglob('*') if p.is_file()
          and p.name != 'MANIFEST.json' and '__pycache__' not in p.parts
          and not p.name.startswith('verification-receipt')}
if actual != set(expected):
    raise RuntimeError('manifest file membership mismatch: ' + str(actual ^ set(expected)))
for name, rec in expected.items():
    data = (root / name).read_bytes()
    if len(data) != rec['bytes'] or hashlib.sha256(data).hexdigest() != rec['sha256']:
        raise RuntimeError('manifest byte mismatch: ' + name)
pins = json.loads((root / 'SOURCE_PINS.json').read_text())
for rec in pins['sources']:
    data = (root / rec['archived_path']).read_bytes()
    if len(data) != rec['bytes'] or hashlib.sha256(data).hexdigest() != rec['sha256']:
        raise RuntimeError('source pin mismatch: ' + rec['archived_path'])
normal = json.loads((root / 'normal.json').read_text())
optimized = json.loads((root / 'optimized.json').read_text())
if normal != optimized or normal['status'] != 'PASS':
    raise RuntimeError('normal/optimized outputs disagree')
if normal['strong_error_order_ten'] != 3992 or normal['total_error_order_ten'] != 5508:
    raise RuntimeError('replay arithmetic mismatch')
distribution = json.loads((root / 'distribution-normal.json').read_text())
if distribution != json.loads((root / 'distribution-optimized.json').read_text()) or distribution['status'] != 'PASS':
    raise RuntimeError('distribution normal/optimized outputs disagree')
if distribution['total_ten_errorful_classes'] != 160 or distribution['total_ten_error'] != 5508:
    raise RuntimeError('distribution arithmetic mismatch')
print(json.dumps({'status':'PASS','files':len(expected),'sources':len(pins['sources']),
                  'normal_optimized_identical':True,'strong_error':normal['strong_error_order_ten'],
                  'total_error':normal['total_error_order_ten']},sort_keys=True))
