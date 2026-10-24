#!/usr/bin/env python3
"""Portable closed LF-normalized manifest verifier; no runtime repository input."""
from pathlib import Path
from hashlib import sha256
import json
HERE=Path(__file__).resolve().parent

def digest(path):
    return sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()

def main():
    data=json.loads((HERE/'MANIFEST.json').read_text())
    expected=data['files']
    actual={str(p.relative_to(HERE)) for p in HERE.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name!='MANIFEST.json'}
    if actual!=set(expected):
        raise RuntimeError('Manifest closure differs: missing='+repr(sorted(set(expected)-actual))+'; extra='+repr(sorted(actual-set(expected))))
    for rel,hexhash in expected.items():
        path=HERE/rel
        if path.resolve().parent!=HERE and HERE not in path.resolve().parents:
            raise RuntimeError('Manifest path leaves package: '+rel)
        if digest(path)!=hexhash:raise RuntimeError('Manifest mismatch: '+rel)
    print('PASS portable closed manifest:',len(expected),'files; UTF-8/LF-normalized SHA256')

if __name__=='__main__':main()
