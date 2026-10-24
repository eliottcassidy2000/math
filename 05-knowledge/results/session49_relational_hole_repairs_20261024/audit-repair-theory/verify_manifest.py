#!/usr/bin/env python3
"""Verify byte hashes and exact file closure, independent of working directory."""
from pathlib import Path
from hashlib import sha256
import json

HERE = Path(__file__).resolve().parent


def main():
    manifest = json.loads((HERE / 'MANIFEST.json').read_text())
    actual = {str(p.relative_to(HERE)) for p in HERE.rglob('*') if p.is_file() and p != HERE / 'MANIFEST.json'}
    if actual != set(manifest['files']):
        raise RuntimeError('Manifest closure differs: ' + repr(sorted(actual ^ set(manifest['files']))))
    for name, row in manifest['files'].items():
        data = (HERE / name).read_bytes()
        if len(data) != row['bytes'] or sha256(data).hexdigest() != row['sha256']:
            raise RuntimeError('Manifest bytes differ: ' + name)
    print('PASS exact manifest closure: ' + str(len(actual)) + ' files')


if __name__ == '__main__':
    main()
