#!/usr/bin/env python3
"""Strict exact-file inventory and SHA-256 verification; valid under -O."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def payload_files(root):
    return {p.relative_to(root).as_posix(): p for p in root.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p != root / 'MANIFEST.json'}


def verify(root=ROOT):
    m = json.loads((root / 'MANIFEST.json').read_text())
    actual = payload_files(root)
    expected = m['files']
    if set(actual) != set(expected):
        raise RuntimeError('manifest inventory mismatch: missing=' + repr(sorted(set(expected)-set(actual))) + '; extra=' + repr(sorted(set(actual)-set(expected))))
    for name, p in actual.items():
        if hashlib.sha256(p.read_bytes()).hexdigest() != expected[name]:
            raise RuntimeError('manifest SHA-256 mismatch: ' + name)
    return len(actual)


if __name__ == '__main__':
    print('VERIFIED ' + str(verify()) + ' SHA-256-pinned payload files')
