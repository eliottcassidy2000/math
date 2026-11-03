#!/usr/bin/env python3
"""Portable exact inventory and SHA-256 gate; no Python assert statements."""
from hashlib import sha256
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parent

def require(condition, message):
    if not condition:
        raise RuntimeError(message)

def verify():
    manifest = json.loads((ROOT / 'MANIFEST.json').read_text())
    require(manifest.get('schema') == 1, 'unknown manifest schema')
    expected = manifest['files']
    actual = sorted(p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*')
                    if p.is_file() and p != ROOT/'MANIFEST.json'
                    and '__pycache__' not in p.parts)
    require(sorted(expected) == actual, 'manifest inventory mismatch')
    for name, item in sorted(expected.items()):
        path = (ROOT / name).resolve()
        require(ROOT in path.parents, 'manifest path escapes root')
        data = path.read_bytes()
        require(len(data) == item['bytes'], 'byte count mismatch: ' + name)
        require(sha256(data).hexdigest() == item['sha256'], 'SHA256 mismatch: ' + name)
    return {'status': 'PASS', 'files_verified': len(actual)}

if __name__ == '__main__':
    print(json.dumps(verify(), sort_keys=True))
