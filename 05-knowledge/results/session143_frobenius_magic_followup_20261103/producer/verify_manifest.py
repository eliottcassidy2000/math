#!/usr/bin/env python3
"""Portable SHA256 integrity gate; manifest itself is not self-hashed."""
from hashlib import sha256
from pathlib import Path
import json


def main():
    root = Path(__file__).resolve().parent
    manifest = json.loads((root/'manifest.json').read_text())
    if manifest.get('schema') != 1:
        raise RuntimeError('unknown manifest schema')
    files = manifest['files']
    actual = sorted(str(p.relative_to(root)) for p in root.rglob('*')
                    if p.is_file() and p.name != 'manifest.json' and '__pycache__' not in p.parts)
    if sorted(files) != actual:
        raise RuntimeError('manifest inventory mismatch')
    for name, expected in sorted(files.items()):
        path = (root/name).resolve()
        if not path.is_relative_to(root):
            raise RuntimeError('manifest path escapes root')
        if sha256(path.read_bytes()).hexdigest() != expected:
            raise RuntimeError('SHA256 mismatch: '+name)
    print(json.dumps({'status':'PASS','files_verified':len(files)},sort_keys=True))


if __name__ == '__main__':
    main()
