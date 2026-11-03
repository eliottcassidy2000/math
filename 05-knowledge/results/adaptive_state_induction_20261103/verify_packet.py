#!/usr/bin/env python3
"""Replay the declared experiments once normally and once with optimization.

No new scan universe is introduced. Internally written receipts are compared
to their frozen bytes, as are stdout receipts. Failures raise under -O too.
"""
from pathlib import Path
import hashlib
import json
import subprocess
import sys

HERE = Path(__file__).resolve().parent
JOBS = [
    ('block_switch.py', 'block_switch.json', False),
    ('prime_boundary.py', 'prime_boundary.json', False),
    ('active_rank_probe.py', 'active_rank_modular.json', False),
    ('active_rank_exact.py', 'active_rank_exact.json', False),
    ('metagraph_reset_probe.py', 'metagraph_reset_probe.json', False),
    ('metagraph_reset_contexts.py', 'metagraph_reset_contexts.json', False),
    ('metagraph_reset_rank_audit.py', 'metagraph_reset_rank_audit.json', False),
    ('return_bank_probe.py', 'return_bank.json', True),
    ('live_port.py', 'live_port.json', False),
    ('live_port_independent_audit.py', 'live_port_independent_audit.json', True),
    ('block_switch_independent_audit.py', 'block_switch_independent_audit.json', True),
    ('two_port_envelope.py', 'two_port_envelope.json', False),
    ('two_port_envelope_independent_audit.py', 'two_port_envelope_independent_audit.json', True),
]


def main():
    receipts = []
    for script, filename, internal in JOBS:
        expected = (HERE / filename).read_bytes()
        for optimized in (False, True):
            command = [sys.executable, '-I', '-B']
            if optimized:
                command += ['-O']
            command += [str(HERE / script)]
            run = subprocess.run(command, cwd=HERE, capture_output=True, check=True)
            actual = (HERE / filename).read_bytes() if internal else run.stdout
            if actual != expected:
                raise RuntimeError(f'{script}: frozen receipt mismatch, optimized={optimized}')
        receipts.append({'script': script, 'receipt': filename,
                         'sha256': hashlib.sha256(expected).hexdigest(),
                         'normal_exact': True, 'optimized_exact': True})
        print(f'PASS {script}', flush=True)
    (HERE / 'verification.json').write_text(json.dumps({
        'status': 'PASS', 'scope': 'Exact replay of thirteen frozen jobs in normal and optimized Python; no expanded universe',
        'jobs': receipts}, indent=2, sort_keys=True) + '\n')


if __name__ == '__main__':
    main()
