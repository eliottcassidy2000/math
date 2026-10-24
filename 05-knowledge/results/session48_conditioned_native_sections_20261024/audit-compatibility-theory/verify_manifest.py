#!/usr/bin/env python3
"""Strict portable byte-exact closure; excludes only root manifest and caches."""
from pathlib import Path
from hashlib import sha256
import json
HERE=Path(__file__).resolve().parent

def main():
    data=json.loads((HERE/'MANIFEST.json').read_text())
    found={str(p.relative_to(HERE)) for p in HERE.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p!=HERE/'MANIFEST.json'}
    expected=data['files']
    if found!=set(expected):raise RuntimeError('Manifest closure mismatch: '+repr(found^set(expected)))
    for relative,entry in expected.items():
        p=HERE/relative
        if HERE not in p.resolve().parents:raise RuntimeError('Out-of-package manifest path')
        if p.stat().st_size!=entry['bytes'] or sha256(p.read_bytes()).hexdigest()!=entry['sha256']:raise RuntimeError('Digest mismatch: '+relative)
    print('PASS closed byte-exact audit manifest:',len(expected),'files')

if __name__=='__main__':main()
