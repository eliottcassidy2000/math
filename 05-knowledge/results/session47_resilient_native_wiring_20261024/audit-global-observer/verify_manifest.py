#!/usr/bin/env python3
import hashlib
import json
from pathlib import Path
P = Path(__file__).resolve().parent

def need(value, reason):
    if not value:
        raise RuntimeError(reason)

manifest = json.loads((P / 'MANIFEST.json').read_text())
for row in manifest['artifacts']:
    q = P / row['path']
    need(q.is_file() and q.stat().st_size == row['bytes'] and hashlib.sha256(q.read_bytes()).hexdigest() == row['sha256'], row['path'])
for row in json.loads((P / 'INPUT_PINS.json').read_text())['inputs']:
    q = P / row['snapshot']
    need(q.stat().st_size == row['bytes'] and hashlib.sha256(q.read_bytes()).hexdigest() == row['sha256'], row['source'])
for a, b, count, key in [('normal.json', 'optimized.json', 494223, 'explicit_checks'),
                         ('producer-replay-normal.json', 'producer-replay-optimized.json', 145809, 'checks')]:
    left, right = (json.loads((P / name).read_text()) for name in (a, b))
    for d in (left, right):
        d.pop('elapsed_seconds', None)
        d.pop('python_optimization', None)
    need(left == right and left['status'] == 'PASS' and left[key] == count, a)
print(json.dumps({'status': 'PASS', 'files_verified': len(manifest['artifacts']),
                  'independent_checks_per_mode': 494223, 'producer_checks_per_mode': 145809}, sort_keys=True))
