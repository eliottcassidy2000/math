#!/usr/bin/env python3
"""Verify the frozen audit and matching normal/optimized results."""
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def need(condition, reason):
    if not condition:
        raise RuntimeError(reason)


def semantic(path):
    data = json.loads(path.read_text())
    data.pop('python_optimization', None)
    data.pop('elapsed_seconds', None)
    return data


manifest = json.loads((HERE / 'MANIFEST.json').read_text())
for row in manifest['artifacts']:
    path = HERE / row['path']
    need(path.is_file(), 'missing file: ' + row['path'])
    need(path.stat().st_size == row['bytes'], 'size mismatch: ' + row['path'])
    need(hashlib.sha256(path.read_bytes()).hexdigest() == row['sha256'], 'hash mismatch: ' + row['path'])
for pin_name in ('INPUT_PINS.json', 'GENERATING_INPUT_PINS.json'):
    for row in json.loads((HERE / pin_name).read_text())['inputs']:
        path = HERE / row['snapshot']
        need(path.stat().st_size == row['bytes'] and hashlib.sha256(path.read_bytes()).hexdigest() == row['sha256'],
             'input mismatch: ' + row['source'])
pairs = [('receipt-normal.json', 'receipt-optimized.json', 2340879, 'explicit_checks'),
         ('generating-normal.json', 'generating-optimized.json', 6740, 'explicit_checks'),
         ('replay-bottleneck-normal.json', 'replay-bottleneck-optimized.json', 2386018, 'explicit_checks'),
         ('replay-suspension-normal.json', 'replay-suspension-optimized.json', 12118, 'checks'),
         ('replay-generating-normal.json', 'replay-generating-optimized.json', 3019, 'checks')]
for normal, optimized, expected, key in pairs:
    left, right = semantic(HERE / normal), semantic(HERE / optimized)
    need(left == right, 'normal/optimized results differ: ' + normal)
    need(left['status'] == 'PASS' and left[key] == expected, 'wrong checker outcome: ' + normal)
print(json.dumps({'status': 'PASS', 'files_verified': len(manifest['artifacts']),
                  'matched_normal_optimized_pairs': len(pairs),
                  'independent_checks_per_mode': 2347619}, sort_keys=True))
