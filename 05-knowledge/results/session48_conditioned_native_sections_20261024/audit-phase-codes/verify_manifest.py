#!/usr/bin/env python3
"""Exact package closure; no symlinks or unlisted files except bytecode caches."""
from pathlib import Path
from hashlib import sha256
import json
HERE=Path(__file__).resolve().parent

def main():
    expected=json.loads((HERE/'MANIFEST.json').read_text())['files']
    paths=list(HERE.rglob('*'))
    if any(p.is_symlink() for p in paths):raise RuntimeError('Symlink in closed audit packet')
    found={p.relative_to(HERE).as_posix() for p in paths if p.is_file() and '__pycache__' not in p.parts and p!=HERE/'MANIFEST.json'}
    if found!=set(expected):raise RuntimeError('Manifest closure mismatch: '+repr(sorted(found^set(expected))))
    for name,entry in expected.items():
        p=HERE/name
        if HERE not in p.resolve().parents:raise RuntimeError('Manifest path escape')
        if p.stat().st_size!=entry['bytes'] or sha256(p.read_bytes()).hexdigest()!=entry['sha256']:raise RuntimeError('Byte mismatch: '+name)
    print('PASS closed byte-exact audit manifest:',len(expected),'files')

if __name__=='__main__':main()
