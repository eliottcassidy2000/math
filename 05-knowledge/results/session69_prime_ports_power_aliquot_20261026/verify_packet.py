#!/usr/bin/env python3
"""Validate exact frozen bytes, then replay isolated copies without updating expectations."""
from pathlib import Path, PurePosixPath
import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def safe_path(name):
    p = PurePosixPath(name)
    need(name == p.as_posix() and not p.is_absolute() and '..' not in p.parts
         and name not in ('', '.'), 'unsafe manifest path: ' + repr(name))
    result = ROOT.joinpath(*p.parts)
    need(not any(q.is_symlink() for q in (result, *result.parents)
                 if q != ROOT and ROOT in q.parents), 'symlink path: ' + name)
    return result


def verify_files():
    manifest = json.loads((ROOT / 'FILE_MANIFEST.json').read_text())
    need(manifest['schema'] == 'session69-frozen-packet-v1', 'manifest schema')
    names = [row['path'] for row in manifest['files']]
    need(names == sorted(set(names)), 'manifest entries must be unique and sorted')
    for row in manifest['files']:
        p = safe_path(row['path'])
        need(p.is_file(), 'missing/nonregular file: ' + row['path'])
        data = p.read_bytes()
        need(len(data) == row['bytes'], 'size mismatch: ' + row['path'])
        need(hashlib.sha256(data).hexdigest() == row['sha256'],
             'SHA256 mismatch: ' + row['path'])
    actual = sorted(p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*')
                    if (p.is_file() or p.is_symlink()) and '__pycache__' not in p.parts
                    and p != ROOT / 'FILE_MANIFEST.json')
    need(actual == names, 'unlisted or absent packet file')
    return manifest


def replay(manifest):
    results = []
    with tempfile.TemporaryDirectory(prefix='session69-replay-') as temp:
        packet = Path(temp) / 'packet'
        shutil.copytree(ROOT, packet, ignore=shutil.ignore_patterns('__pycache__'))
        for job in manifest['replays']:
            script = safe_path(job['script']).relative_to(ROOT)
            expected = safe_path(job['normal']).read_bytes()
            need(expected == safe_path(job['optimized']).read_bytes(),
                 'native normal/optimized mismatch: ' + job['script'])
            args = [str(packet / safe_path(a).relative_to(ROOT))
                    if a.startswith('session69-') else a for a in job.get('args', [])]
            for optimized in (False, True):
                flags = ['-I', '-B'] + (['-O'] if optimized else [])
                proc = subprocess.run([sys.executable, *flags, str(packet / script), *args],
                                      cwd=packet / script.parent, capture_output=True)
                need(proc.returncode == 0, 'replay failed: ' + job['script'] + '\n'
                     + proc.stderr.decode(errors='replace'))
                need(proc.stdout == expected, 'replay receipt mismatch: ' + job['script'])
            results.append({'script': job['script'], 'normal': 'PASS', 'optimized': 'PASS',
                            'receipt_sha256': hashlib.sha256(expected).hexdigest()})
        for job in manifest.get('negative_controls', []):
            script = safe_path(job['script']).relative_to(ROOT)
            for optimized in (False, True):
                flags = ['-I', '-B'] + (['-O'] if optimized else [])
                proc = subprocess.run([sys.executable, *flags, str(packet / script), *job['args']],
                                      cwd=packet / script.parent, capture_output=True)
                need(proc.returncode == job['exit_code'] and not proc.stdout
                     and job['stderr_contains'].encode() in proc.stderr,
                     'negative control missed intended gate: ' + job['script'])
            results.append({'script': job['script'], 'negative_control': 'PASS',
                            'normal_and_optimized': True})
        # A native cross-check may write its own input ledger. It must reproduce
        # the original ledger, not silently mutate a frozen expectation.
        for row in manifest['files']:
            p = packet / row['path']
            need(p.is_file() and hashlib.sha256(p.read_bytes()).hexdigest() == row['sha256'],
                 'replay changed a frozen input/expectation: ' + row['path'])
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--replay', action='store_true')
    args = parser.parse_args()
    manifest = verify_files()
    print(json.dumps({'status': 'PASS', 'files_verified': len(manifest['files']),
                      'source_pin': manifest['source_pin'],
                      'replays': replay(manifest) if args.replay else 'NOT_RUN'},
                     indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
