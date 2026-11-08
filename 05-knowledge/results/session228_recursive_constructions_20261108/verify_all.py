#!/usr/bin/env python3
"""Portable exact producer/auditor replays; every gate survives optimization."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
CASES = [
    ('mobile_triangle/verify.py', 'mobile_triangle/CHECK.json'),
    ('mobile_triangle/verify_membership.py', 'mobile_triangle/MEMBERSHIP_CHECK.json'),
    ('mobile_triangle/verify_multi.py', 'mobile_triangle/MULTI_CHECK.json'),
    ('mobile_audit/verify.py', 'mobile_audit/normal.json'),
    ('mobile_audit/classify.py', 'mobile_audit/classify.normal.json'),
    ('mobile_audit/extensions/check.py', 'mobile_audit/extensions/normal.json'),
    ('cyclic_inverse/verify.py', 'cyclic_inverse/CHECK.json'),
    ('cyclic_audit/check.py', 'cyclic_audit/normal.out'),
    ('recursive/verify.py', 'recursive/CHECK.json'),
    ('recursive_audit/audit.py', 'recursive_audit/AUDIT.json'),
    ('physical_state/verify.py', 'physical_state/receipt.json'),
    ('physical_state/verify_oracle.py', 'physical_state/oracle-receipt.json'),
    ('physical_state_audit/audit.py', 'physical_state_audit/normal.json'),
]

def main():
    results = []
    for script, receipt in CASES:
        expected = (ROOT / receipt).read_bytes()
        for optimized in (False, True):
            cmd = [sys.executable, '-B'] + (['-O'] if optimized else []) + [str(ROOT / script)]
            proc = subprocess.run(cmd, cwd=(ROOT / script).parent, stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE, check=False)
            if proc.returncode != 0:
                raise RuntimeError(f'{script}, optimized={optimized}: {proc.stderr.decode()}')
            if proc.stdout != expected:
                raise RuntimeError(f'{script}, optimized={optimized}: deterministic receipt mismatch')
            results.append({'script': script, 'optimized': optimized,
                            'receipt_sha256': hashlib.sha256(expected).hexdigest()})
    print(json.dumps({'status': 'PASS', 'runs': len(results), 'results': results},
                     sort_keys=True, indent=2))

if __name__ == '__main__':
    main()
