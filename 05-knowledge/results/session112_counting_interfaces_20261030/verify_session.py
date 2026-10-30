#!/usr/bin/env python3
"""Replay every Session 112 finite checker, normally and with optimization."""
from pathlib import Path
import json
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
CASES = (
    ('arithmetic/verify_modular.py', 'arithmetic/modular_checks.json'),
    ('arithmetic/verify_aliquot_modular.py', 'arithmetic/aliquot_modular_checks.json'),
    ('arithmetic/verify_power_observers.py', 'arithmetic/power_checks.json'),
    ('norms/verify_norms.py', 'norms/checks.json'),
    ('norms/verify_defect.py', 'norms/defect_checks.json'),
    ('symmetry/check.py', 'symmetry/results.json'),
    ('symmetry/prime_check.py', 'symmetry/prime_results.json'),
    ('counting/check.py', 'counting/results.json'),
    ('counting/bank_rank.py', 'counting/bank_results.json'),
)


def main():
    for script, receipt in CASES:
        expected = json.loads((ROOT / receipt).read_text(encoding='utf-8'))
        if expected.get('status') != 'PASS':
            raise RuntimeError(f'{receipt}: published result is not PASS')
        for optimized in (False, True):
            command = [sys.executable, '-B']
            if optimized:
                command.append('-O')
            command.append(str(ROOT / script))
            result = subprocess.run(command, cwd=(ROOT / script).parent,
                                    text=True, capture_output=True, check=False)
            mode = 'optimized' if optimized else 'normal'
            if result.returncode:
                raise RuntimeError(f'{script} ({mode}) failed: {result.stderr}')
            try:
                actual = json.loads(result.stdout)
            except json.JSONDecodeError as exc:
                raise RuntimeError(f'{script} ({mode}) did not emit JSON') from exc
            if actual != expected:
                raise RuntimeError(f'{script} ({mode}) differs from {receipt}')
        print(f'PASS {script}: normal and optimized match', flush=True)
    print(f'PASS: {len(CASES)} mathematical checkers, {2 * len(CASES)} executions')


if __name__ == '__main__':
    main()
