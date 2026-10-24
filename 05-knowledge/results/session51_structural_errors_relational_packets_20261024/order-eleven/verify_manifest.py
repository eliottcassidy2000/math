#!/usr/bin/env python3
"""Strict relocatable membership, source pin, receipt and manifest verification."""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parent
EXEMPT = {'MANIFEST.json', 'verification-receipt.json', 'relocation-receipt.json'}

def demand(test, message):
    if not test:
        raise RuntimeError(message)

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    manifest = json.loads((ROOT / 'MANIFEST.json').read_text())
    actual = {p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.relative_to(ROOT).as_posix() not in EXEMPT}
    expected = {r['path'] for r in manifest['files']}
    demand(actual == expected, 'exact packet membership mismatch')
    for r in manifest['files']:
        path = ROOT / r['path']
        demand(path.stat().st_size == r['bytes'], 'size mismatch: ' + r['path'])
        demand(sha(path) == r['sha256'], 'digest mismatch: ' + r['path'])
    sources = json.loads((ROOT / 'SOURCE_PINS.json').read_text())
    for r in sources['sources']:
        path = ROOT / r['archived_path']
        demand(path.stat().st_size == r['bytes'] and sha(path) == r['sha256'], 'archived source mismatch')
    normal = ROOT / 'normal.json'
    optimized = ROOT / 'optimized.json'
    demand(normal.read_bytes() == optimized.read_bytes(), 'normal and optimized disagree')
    data = json.loads(normal.read_text())
    demand(data['status'] == 'PASS' and data['nonregular_strong_error'] == 9166 and data['known_total_error'] == 19480, 'result arithmetic mismatch')
    demand(data['unknown']['positive_regular_eleven_error_support'] == [22,44,66,88,110], 'unknown regular support mismatch')
    demand(data['upper_trace_modulo_eleven'] == 10, 'prospective residue mismatch')
    return {'status':'PASS', 'files':len(actual), 'sources':len(sources['sources']), 'normal_optimized_identical':True,
      'normal_sha256':sha(normal), 'nonregular_strong_error':9166, 'known_total_error':19480,
      'regular_eleven_error':'not computed by this packet', 'full_order_eleven_count':'not computed'}

if __name__ == '__main__':
    print(json.dumps(main(), indent=2, sort_keys=True))
