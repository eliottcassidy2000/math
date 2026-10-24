#!/usr/bin/env python3
"""Replay from only the manifest closure, without ambient repository files."""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent


def main():
    manifest = json.loads((ROOT / 'MANIFEST.json').read_text())
    with tempfile.TemporaryDirectory(prefix='session49-nonreturn-manifest-') as scratch:
        destination = Path(scratch)
        for name in list(manifest['files']) + ['MANIFEST.json']:
            target = destination / name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(ROOT / name, target)
        verify = subprocess.run([sys.executable, '-B', 'verify_manifest.py'], cwd=destination,
                                check=True, capture_output=True)
        normal = subprocess.run([sys.executable, '-B', 'check.py'], cwd=destination,
                                check=True, capture_output=True)
        optimized = subprocess.run([sys.executable, '-O', '-B', 'check.py'], cwd=destination,
                                   check=True, capture_output=True)
        if normal.stdout != optimized.stdout or normal.stdout != (destination / 'receipt.json').read_bytes():
            raise RuntimeError('portable normal, optimized and frozen receipts differ')
        if verify.stderr or normal.stderr or optimized.stderr:
            raise RuntimeError('unexpected replay stderr')
        print(json.dumps(dict(status='passed', manifest_only_files=len(manifest['files']) + 1,
                              python_modes=['normal', 'optimized'],
                              receipt_stdout_sha256=hashlib.sha256(normal.stdout).hexdigest()), sort_keys=True))


if __name__ == '__main__':
    main()
