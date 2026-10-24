#!/usr/bin/env python3
"""Separate parent-only input validator; never imported by blind predictor."""
import hashlib
import importlib.machinery
import importlib.util
import json
import math
from pathlib import Path


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def main():
    here = Path(__file__).resolve().parent
    source = here / 'archive/growth-transfer.py.txt'
    pins = json.loads((here / 'source-pins.json').read_text())
    require(hashlib.sha256(source.read_bytes()).hexdigest() == pins['sources'][0]['sha256'],
            'pinned input-generation source')
    loader = importlib.machinery.SourceFileLoader('frozen_parent_validator', str(source))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    payload = json.loads((here / 'parents-through-seven.json').read_text())
    reports = []
    for level in payload['levels']:
        n = level['order']
        keys = set()
        mass = 0
        for record in level['parents']:
            rows = tuple(record['rows'])
            actions = tuple(tuple(g) for g in record['actions'])
            require(set(actions) == set(module.automorphisms(rows)),
                    'supplied actions exhaust parent automorphism group')
            key = module.canonical(rows)[0]
            require(key not in keys, 'parent representatives are pairwise nonisomorphic')
            keys.add(key)
            mass += math.factorial(n) // len(actions)
        require(mass == 1 << (n * (n-1) // 2), 'full labelled tournament mass at parent order')
        reports.append({'order': n, 'classes': len(keys), 'labelled_mass': mass})
    require(max(module.CANON_ORDERS, default=0) <= 7, 'validator never sees target-order objects')
    print(json.dumps({'status': 'PASS', 'scope': 'input validation only, separate from blind prediction',
                      'all_parent_actions_complete': True,
                      'all_parent_atlases_pairwise_nonisomorphic_and_complete': True,
                      'maximum_canonicalized_order': max(module.CANON_ORDERS, default=0),
                      'levels': reports}, sort_keys=True))


if __name__ == '__main__':
    main()
