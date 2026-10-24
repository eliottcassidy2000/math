#!/usr/bin/env python3
"""Verify every frozen artifact hash and matching normal/optimized receipts."""
from pathlib import Path
import hashlib,json
H=Path(__file__).resolve().parent;M=json.loads((H/'MANIFEST.json').read_text())
for name,expected in M['artifacts'].items():
 actual=hashlib.sha256((H/name).read_bytes()).hexdigest()
 if actual!=expected:raise RuntimeError('Artifact hash mismatch: '+name)
for row in M['receipt_pairs']:
 a=(H/row['normal']).read_bytes();b=(H/row['optimized']).read_bytes()
 if a!=b:raise RuntimeError('Receipt mismatch: '+row['normal'])
 r=json.loads(a)
 if r['status']!='PASS' or r['checks']!=row['checks']:raise RuntimeError('Receipt status/count mismatch')
print(json.dumps(dict(status='PASS',artifacts=len(M['artifacts']),receipt_pairs=len(M['receipt_pairs'])),sort_keys=True))
