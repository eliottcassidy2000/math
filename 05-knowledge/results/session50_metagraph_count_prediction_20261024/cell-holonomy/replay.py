#!/usr/bin/env python3
"""Manifest-only isolated normal and optimized replay."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

root = Path(__file__).resolve().parent
manifest = json.loads((root / 'manifest.json').read_text())
for relative, expected in manifest['files'].items():
    actual = hashlib.sha256((root / relative).read_bytes()).hexdigest()
    if actual != expected:
        raise RuntimeError('SHA256 mismatch: ' + relative)
with tempfile.TemporaryDirectory(prefix='session50-cell-holonomy-') as directory:
    target = Path(directory)
    for relative in manifest['files']:
        destination = target / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(root / relative, destination)
    expected_receipt = (root / 'receipt.json').read_bytes()
    for optimized in (False, True):
        command = [sys.executable] + (['-O'] if optimized else []) + ['check.py']
        result = subprocess.run(command, cwd=target, check=True, capture_output=True)
        if (target / 'receipt.json').read_bytes() != expected_receipt:
            raise RuntimeError('isolated receipt changed')
        expected_output = (root / ('optimized.out' if optimized else 'normal.out')).read_bytes()
        if result.stdout != expected_output:
            raise RuntimeError('isolated stdout changed')
        if result.stderr:
            raise RuntimeError('unexpected isolated stderr')
print(json.dumps({'result': 'PASS', 'manifest_files': len(manifest['files']),
                  'isolated_runs': ['normal', 'optimized'],
                  'receipt_sha256': hashlib.sha256(expected_receipt).hexdigest()},
                 sort_keys=True))
