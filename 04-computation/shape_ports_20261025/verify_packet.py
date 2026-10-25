#!/usr/bin/env python3
"""Replay the scoped research controls against the frozen JSON receipts.

Default: normal Python. --optimized: optimized Python. --both: both modes.
Only named source files are executed; no inherited producer code is imported.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
RESULT = ROOT / '05-knowledge/results/shape_ports_20261025'
CASES = [
 ('shape_ports_20261025/low_order_detours.py', 'tournaments_low_order_detours.json', True),
 ('shape_ports_20261025/weighted_endpoint_boundary.py', 'tournaments_weighted_endpoint_boundary.json', True),
 ('shape_ports_20261025/endpoint_four_seven_control.py', 'tournaments_endpoint_four_seven_control.json', False),
 ('shape_anchors_20261025/verify_quartic_sextic.py', 'anchors_quartic_sextic_receipt.json', False),
 ('shape_anchors_20261025/anchors_independent_audit.py', 'anchors_audit.json', False),
 ('shape_anchors_20261025/verify_collatz_phase.py', 'anchors_collatz_phase_receipt.json', False),
 ('shape_anchors_20261025/collatz27_native_refill.py', 'collatz_native_refill.json', True),
 ('shape_anchors_20261025/collatz_independent_audit.py', 'collatz_audit.json', False),
 ('shape_spectral_20261025/verify_shapes.py', 'shapes-verification.json', False),
]


def main():
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--optimized', action='store_true')
    group.add_argument('--both', action='store_true')
    args = parser.parse_args()
    modes = [False, True] if args.both else [args.optimized]
    manifest = json.loads((RESULT / 'MANIFEST.json').read_text())
    for relative, expected in manifest['sha256'].items():
        actual = hashlib.sha256((ROOT / relative).read_bytes()).hexdigest()
        if actual != expected:
            raise RuntimeError('Artifact hash mismatch: ' + relative)
    with tempfile.TemporaryDirectory(prefix='math-shape-replay-') as scratch:
        for optimized in modes:
            for source, receipt, writes_file in CASES:
                command = [sys.executable, '-B'] + (['-O'] if optimized else [])
                command += [str(ROOT / '04-computation' / source)]
                output = Path(scratch) / receipt
                if writes_file:
                    command += ['--output', str(output)]
                run = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
                if run.returncode:
                    raise RuntimeError(f'{source} failed:\n{run.stderr}')
                produced = json.loads(output.read_text() if writes_file else run.stdout)
                expected = json.loads((RESULT / receipt).read_text())
                if produced != expected:
                    raise RuntimeError('Receipt mismatch: ' + source)
                print(('optimized' if optimized else 'normal') + ' PASS ' + source, flush=True)
    print(f'PASS: {len(CASES) * len(modes)} receipt replays; artifact hashes verified.')


if __name__ == '__main__':
    main()
