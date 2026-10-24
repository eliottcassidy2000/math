#!/usr/bin/env python3
"""Verify the exact delivered packet without rewriting it."""
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def need(value, message):
    if not value:
        raise RuntimeError(message)


def main():
    manifest = json.loads((HERE / 'MANIFEST.json').read_text())
    for row in manifest['files']:
        path = HERE / row['path']
        need(path.is_file(), f'Missing file: {row["path"]}')
        data = path.read_bytes()
        need(len(data) == row['bytes'], f'Size mismatch: {row["path"]}')
        need(hashlib.sha256(data).hexdigest() == row['sha256'],
             f'Hash mismatch: {row["path"]}')
    left = json.loads((HERE / 'receipt-normal.json').read_text())
    right = json.loads((HERE / 'receipt-optimized.json').read_text())
    need(left['python_optimization'] == 0 and right['python_optimization'] == 1,
         'Receipts do not cover normal and optimized execution')
    for receipt in (left, right):
        need(receipt['status'] == 'PASS', 'Receipt did not pass')
        need(receipt['checker_sha256'] == hashlib.sha256((HERE / 'check.py').read_bytes()).hexdigest(),
             'Receipt does not match current checker')
        need(receipt['source_inputs']['sourcepins_sha256'] == hashlib.sha256((HERE / 'SOURCEPINS.json').read_bytes()).hexdigest(),
             'Receipt does not match current source pins')
    for key in ('elapsed_seconds', 'python_optimization'):
        left.pop(key)
        right.pop(key)
    need(left == right, 'Normal and optimized mathematical results differ')
    print(json.dumps({'status': 'PASS', 'files_verified': len(manifest['files']),
                      'normal_optimized_results_identical': True,
                      'checks_per_full_run': left['explicit_checks']}, sort_keys=True))


if __name__ == '__main__':
    main()
