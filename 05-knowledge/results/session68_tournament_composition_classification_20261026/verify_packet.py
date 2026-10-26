#!/usr/bin/env python3
"""Verify frozen packet bytes, then optionally replay in an isolated copy.

The manifest is an integrity ledger, not an external trust anchor. Replays
never update the frozen expectations. All gates survive Python -O.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def safe_path(name):
    p = PurePosixPath(name)
    require(name == p.as_posix() and not p.is_absolute()
            and '..' not in p.parts and name not in ('', '.'),
            'unsafe manifest path: ' + repr(name))
    return ROOT.joinpath(*p.parts)


def verify_files():
    manifest = json.loads((ROOT / 'FILE_MANIFEST.json').read_text())
    require(manifest['schema'] == 'session68-frozen-packet-v1', 'manifest schema')
    rows = manifest['files']
    names = [row['path'] for row in rows]
    require(names == sorted(set(names)), 'manifest entries must be unique and sorted')
    for row in rows:
        p = safe_path(row['path'])
        require(p.is_file() and not p.is_symlink(), 'missing/nonregular file: ' + row['path'])
        data = p.read_bytes()
        require(len(data) == row['bytes'], 'size mismatch: ' + row['path'])
        require(hashlib.sha256(data).hexdigest() == row['sha256'],
                'SHA256 mismatch: ' + row['path'])
    actual = sorted(p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*')
                    if p.is_file() and '__pycache__' not in p.parts
                    and p.name != 'FILE_MANIFEST.json')
    # Nested native manifests also have evidence value and must be indexed.
    actual.extend(p.relative_to(ROOT).as_posix() for p in ROOT.rglob('FILE_MANIFEST.json')
                  if p != ROOT / 'FILE_MANIFEST.json')
    require(sorted(actual) == names, 'unlisted or absent packet file')
    return manifest


def replay(manifest):
    results = []
    with tempfile.TemporaryDirectory(prefix='session68-portable-') as temp:
        packet = Path(temp) / 'packet'
        shutil.copytree(ROOT, packet, ignore=shutil.ignore_patterns('__pycache__'))
        for job in manifest['replays']:
            script = safe_path(job['script']).relative_to(ROOT)
            expected = (ROOT / job['normal']).read_bytes()
            require(expected == (ROOT / job['optimized']).read_bytes(),
                    'native normal/-O receipts differ: ' + job['script'])
            for optimized in (False, True):
                flags = ['-I', '-B'] + (['-O'] if optimized else [])
                proc = subprocess.run([sys.executable, *flags, str(packet / script)],
                                      cwd=packet / script.parent,
                                      capture_output=True, check=False)
                require(proc.returncode == 0,
                        'replay failed: ' + job['script'] + '\n' + proc.stderr.decode(errors='replace'))
                require(proc.stdout == expected,
                        'replay output mismatch: ' + job['script'] + (' (-O)' if optimized else ''))
            results.append({'script': job['script'], 'normal': 'PASS', 'optimized': 'PASS',
                            'receipt_sha256': hashlib.sha256(expected).hexdigest()})
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--replay', action='store_true')
    args = parser.parse_args()
    manifest = verify_files()
    result = {'status': 'PASS', 'files_verified': len(manifest['files']),
              'source_pin': manifest['source_pin'],
              'replays': replay(manifest) if args.replay else 'NOT_RUN'}
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
