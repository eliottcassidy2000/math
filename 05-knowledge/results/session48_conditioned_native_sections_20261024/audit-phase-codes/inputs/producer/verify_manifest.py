#!/usr/bin/env python3
"""Verify exact standalone package closure and LF-normalized file hashes."""
from pathlib import Path
from hashlib import sha256
import json
HERE=Path(__file__).resolve().parent

def main():
    expected=json.loads((HERE/'MANIFEST.json').read_text())['files']
    entries=list(HERE.rglob('*'))
    if any(p.is_symlink() for p in entries):raise RuntimeError('Package contains a symlink')
    actual={p.relative_to(HERE).as_posix() for p in entries if p.is_file() and '__pycache__' not in p.parts and p.relative_to(HERE).as_posix()!='MANIFEST.json'}
    if actual!=set(expected):raise RuntimeError('Manifest closure differs: missing='+repr(sorted(set(expected)-actual))+'; extra='+repr(sorted(actual-set(expected))))
    for name,want in expected.items():
        path=HERE/name
        if HERE not in path.resolve().parents:raise RuntimeError('Manifest path escapes package: '+name)
        got=sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
        if got!=want:raise RuntimeError('Manifest hash mismatch: '+name)
    print('PASS closed manifest:',len(expected),'files; LF-normalized SHA256; no external paths')

if __name__=='__main__':main()
