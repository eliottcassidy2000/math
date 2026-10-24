#!/usr/bin/env python3
"""Verify exact raw-byte closure and the vendored input source pins."""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parent


def verify():
    manifest = json.loads((ROOT / 'MANIFEST.json').read_text())
    expected = manifest['files']
    observed = {str(p.relative_to(ROOT)) for p in ROOT.rglob('*')
                if p.is_file() and '__pycache__' not in p.parts and p != ROOT / 'MANIFEST.json'}
    if observed != set(expected):
        raise RuntimeError('manifest file closure mismatch: ' + repr(sorted(observed ^ set(expected))))
    for name, digest in expected.items():
        path = ROOT / name
        if path.is_symlink() or ROOT not in path.resolve().parents:
            raise RuntimeError('unsafe manifest path: ' + name)
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise RuntimeError('exact-byte manifest mismatch: ' + name)
    pins = json.loads((ROOT / 'SOURCEPINS.json').read_text())
    for name, info in pins['files'].items():
        data = (ROOT / name).read_bytes()
        if len(data) != info['bytes'] or hashlib.sha256(data).hexdigest() != info['sha256']:
            raise RuntimeError('source pin mismatch: ' + name)
    return len(expected), len(pins['files'])


if __name__ == '__main__':
    files, archives = verify()
    print(json.dumps(dict(status='verified', manifest_files=files, exact_byte_source_archives=archives), sort_keys=True))
