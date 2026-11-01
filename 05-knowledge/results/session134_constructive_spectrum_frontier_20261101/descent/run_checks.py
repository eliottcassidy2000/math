#!/usr/bin/env python3
"""Reproduce every final receipt and compare optimization modes."""
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import sys


def require(ok, label):
    if not ok:
        raise RuntimeError(label)


def main():
    root = Path(__file__).resolve().parent
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
    records = []
    parsed = {}
    for script, prefix in [('verify.py', ''), ('independent_check.py', 'independent.'),
                           ('family_check.py', 'family.')]:
        outputs = []
        for mode, flags in [('normal', []), ('optimized', ['-O'])]:
            result = subprocess.run([sys.executable, *flags, str(root / script)],
                                    cwd=root, env=env, capture_output=True, check=False)
            require(result.returncode == 0, (script, mode, result.stderr.decode()))
            require(not result.stderr, (script, mode, 'unexpected stderr'))
            receipt = root / f'{prefix}{mode}.json'
            require(result.stdout == receipt.read_bytes(), (script, mode, 'retained receipt mismatch'))
            value = json.loads(result.stdout)
            require(value['status'] == 'PASS', (script, mode, 'not PASS'))
            require(value['source_sha256'] == sha256((root / script).read_bytes()).hexdigest(),
                    (script, mode, 'source hash mismatch'))
            outputs.append(result.stdout)
        require(outputs[0] == outputs[1], (script, 'normal versus optimized mismatch'))
        parsed[script] = json.loads(outputs[0])
        records.append({'script': script, 'modes': ['normal', 'optimized'],
                        'receipt_sha256': sha256(outputs[0]).hexdigest(),
                        'source_sha256': sha256((root / script).read_bytes()).hexdigest()})
    require(parsed['verify.py']['one_vertex_pairs_by_half_gain'] ==
            parsed['independent_check.py']['one_vertex_pairs_by_half_gain'],
            'independent pairsets disagree')
    require(parsed['family_check.py']['literal_source_state'] == [33, 32, 2, 9],
            'actual source state')
    verification = {'status': 'PASS', 'records': records,
                    'independent_pairsets_agree': True,
                    'optimization_mode_outputs_identical': True,
                    'source_realizability_verified': True,
                    'all_order_scope': 'Theorems in THEORY.md; no global spectrum coverage claim.',
                    'harness_sha256': sha256(Path(__file__).read_bytes()).hexdigest()}
    encoded = json.dumps(verification, indent=2, sort_keys=True) + '\n'
    (root / 'verification.json').write_text(encoded)
    print(encoded, end='')

if __name__ == '__main__':
    main()
