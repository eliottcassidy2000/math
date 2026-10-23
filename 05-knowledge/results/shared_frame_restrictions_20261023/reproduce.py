#!/usr/bin/env python3
"""Replay the frozen packet; do not regenerate expected results on failure."""
from hashlib import sha256
from pathlib import Path
import json
import subprocess
import sys

root = Path(__file__).resolve().parent
manifest = json.loads((root/'manifest.json').read_text())
for name, digest in manifest['files'].items():
    if sha256((root/name).read_bytes()).hexdigest() != digest:
        raise RuntimeError('hash mismatch: '+name)
for script, receipt in [('check.py','receipt.json'), ('native.py','native_receipt.json'),
                        ('algebra.py','algebra_receipt.json')]:
    expected = (root/receipt).read_bytes()
    for mode in ([], ['-O']):
        actual = subprocess.check_output([sys.executable, *mode, str(root/script)])
        if actual != expected:
            raise RuntimeError('receipt mismatch: '+script+' '+str(mode))
    print(script+': normal/-O exact replay PASS', flush=True)
print('All packet hashes and receipts PASS')
