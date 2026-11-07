#!/usr/bin/env python3
"""Replay the small exact certificates; this does not assert target coverage."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
for line in (ROOT / 'SHA256SUMS').read_text().splitlines():
    digest, relative = line.split('  ', 1)
    actual = hashlib.sha256((ROOT / relative).read_bytes()).hexdigest()
    if actual != digest:
        raise RuntimeError(f'Hash mismatch: {relative}')

cases = [
    ('extrema/verify.py', 'extrema/VERIFY.json'),
    ('selectors/verify.py', 'selectors/RESULT.json'),
    ('audit/audit.py', 'audit/RESULT.json'),
    ('audit/audit_cylinders.py', 'audit/CYLINDERS.json'),
    ('audit/ribbon/verify_audit.py', 'audit/ribbon/RESULT.json'),
    ('exploration/verify_switch.py', 'exploration/RESULT.json'),
]
for script, receipt in cases:
    proc = subprocess.run([sys.executable, str(ROOT / script)], cwd=ROOT,
                          capture_output=True, text=True, check=True)
    actual = json.loads(proc.stdout)
    expected = json.loads((ROOT / receipt).read_text())
    if actual != expected:
        raise RuntimeError(f'Receipt mismatch: {script}')
    print(f'PASS {script}', flush=True)
proc = subprocess.run([sys.executable, str(ROOT / 'ribbon/verify.py')], cwd=ROOT,
                      capture_output=True, text=True, check=True)
if proc.stdout != (ROOT / 'ribbon/verify.log').read_text():
    raise RuntimeError('Receipt mismatch: ribbon/verify.py')
print('PASS ribbon/verify.py: full finite census replay', flush=True)
print('PASS: seven normal-mode exact replays; all-target coverage remains OPEN.')
