#!/usr/bin/env python3
"""Check copied-source hashes and the literal primitive ROOT affine edge."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
pins = json.loads((ROOT / 'SOURCE_PINS.json').read_text())
for source in pins['sources']:
    if hashlib.sha256((ROOT / source['local_copy']).read_bytes()).hexdigest() != source['sha256']:
        raise RuntimeError(source['path'])

count = 0
for s in range(1, 100002, 2):
    if s % 8 == 5:
        continue
    a = 2 if s % 8 == 1 else 1
    y = (3 * s + 1) // (2 ** a)
    j = 0
    # B(z)=4z+1 has positive odd inverse exactly on residue5 mod8.
    while y % 8 == 5:
        y = (y - 1) // 4
        j += 1
    t = a + 2 * j
    if y <= 0 or y % 2 != 1 or y % 8 == 5:
        raise RuntimeError(('wrong primitive target', s, y))
    if 2 ** t * (3 * y + 1) != 3 * (3 * s + 1) + 2 ** a:
        raise RuntimeError(('wrong affine carry', s, y, a, t))
    count += 1

print(json.dumps({
    'base': pins['source_head'],
    'literal_native_edges_checked': count,
    'maximum_seed': 100001,
    'source_hashes_verified': len(pins['sources']),
    'first3_witness': {'s': 3, 'y': 1, 'a': 1, 't': 3, 'left': 32, 'right': 32},
    'scope': 'Source recovery and literal affine edge check only; no all-depth energy or entry theorem.'
}, indent=2, sort_keys=True))
