#!/usr/bin/env python3
"""Replay the bounded six/seven continuation, checking its frozen bytes."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
CODE = ROOT / '04-computation/six_seven_gale_cube_20261025'
RESULT = ROOT / '05-knowledge/results/six_seven_gale_cube_20261025'
CASES = [
    ('verify_budget_gale.py', 'budget_gale.json', False),
    ('check_character_decomposition.py', 'characters.json', False),
    ('check_geometry.py', 'geometry.json', False),
    ('check_free_deck_codes.py', 'free_deck_codes.json', True),
    ('audit_character.py', 'character_audit.json', False),
    ('audit_free_quotients.py', 'free_quotients_audit.json', False),
    ('check_reflection_markings.py', 'reflection_markings.json', True),
]


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def main():
    manifest = json.loads((RESULT / 'MANIFEST.json').read_text())
    for rel, digest in manifest['sha256'].items():
        actual = hashlib.sha256((ROOT / rel).read_bytes()).hexdigest()
        require(actual == digest, f'Hash mismatch: {rel}')
    rows = []
    with tempfile.TemporaryDirectory(prefix='six-seven-replay-') as temp:
        for script, receipt, writes_file in CASES:
            expected = (RESULT / receipt).read_bytes()
            json.loads(expected)
            observed = []
            for optimized in (False, True):
                command = [sys.executable, '-B']
                if optimized:
                    command.append('-O')
                command.append(str(CODE / script))
                output = Path(temp) / receipt
                if writes_file:
                    command.extend(['--output', str(output)])
                process = subprocess.run(command, cwd=ROOT, capture_output=True,
                                         check=False, timeout=180)
                require(process.returncode == 0,
                        f'Checker failed: {script}: {process.stderr.decode()}')
                data = output.read_bytes() if writes_file else process.stdout
                require(json.loads(data) == json.loads(expected),
                        f'Receipt mismatch: {script}, optimized={optimized}')
                require(data == expected,
                        f'Byte mismatch: {script}, optimized={optimized}')
                observed.append(data)
            require(observed[0] == observed[1], f'Optimization mismatch: {script}')
            rows.append({'script': script, 'normal': 'PASS', 'optimized': 'PASS'})
    print(json.dumps({'status': 'PASS', 'manifest_files': len(manifest['sha256']),
                      'executions': 2 * len(rows), 'replays': rows},
                     indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
