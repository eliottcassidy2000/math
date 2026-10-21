#!/usr/bin/env python3
"""Authenticate the packet, then compare both normal/optimized consequences."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from palette_check import repo_root, need


def main():
    m = json.loads((HERE/'manifest.json').read_text()); root = repo_root()
    for name, digest in m['dependencies'].items():
        need(hashlib.sha256((root/name).read_bytes()).hexdigest() == digest, 'dependency '+name)
    for name, digest in m['files'].items():
        need(hashlib.sha256((HERE/name).read_bytes()).hexdigest() == digest, 'packet '+name)
    result = []
    for checker, receipt in [('zero_heavy_check.py','zero-heavy-check.json'),
                             ('palette_check.py','palette-check.json')]:
        expected = (HERE/receipt).read_bytes()
        for mode in ([], ['-O']):
            got = subprocess.check_output([sys.executable]+mode+[str(HERE/checker)],cwd=root)
            need(got == expected, 'replay '+checker+' '+str(mode))
            result.append({'checker':checker, 'optimized':bool(mode), 'bytes':len(got),
                           'sha256':hashlib.sha256(got).hexdigest()})
    print(json.dumps({'status':'PASS', 'runs':result},indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
