#!/usr/bin/env python3
"""Replay the original 24-file random draw and check frozen source snapshots.

This script only reads text/bytes; none of the sampled historical code runs.
"""
from pathlib import Path
import hashlib
import json
import random

ROOT = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def main():
    arithmetic = ROOT / 'session68-arithmetic'
    source = ROOT / 'random-source-snapshots'
    draw = json.loads((arithmetic / 'random_manifest.json').read_text())
    manifest = json.loads((source / 'SOURCE_MANIFEST.json').read_text())
    require(draw['seed'] == 2026102668, 'seed changed')
    require(draw['pin'] == manifest['source_pin']
            == 'fc3768ac5fcd87022fde9e6d7e3f06c0e5bc22e8', 'source pin changed')
    rng = random.Random(draw['seed'])
    selected = []
    native_blob_hashes = {}
    rows = []
    for stratum in draw['strata']:
        data = (arithmetic / (stratum['name'] + '_universe.txt')).read_bytes()
        require(hashlib.sha256(data).hexdigest() == stratum['universe_sha256'],
                'universe digest: ' + stratum['name'])
        paths = data.decode().splitlines()
        require(paths == sorted(set(paths)), 'universe not sorted and unique')
        require(len(paths) == stratum['count'], 'universe size mismatch')
        require(all(p.startswith(stratum['prefix']) and p.endswith(stratum['suffix'])
                    for p in paths), 'stratum path filter')
        actual = rng.sample(paths, 6)
        require(actual == stratum['draws'], 'random draw mismatch')
        selected.extend(actual)
        native_blob_hashes.update(stratum['draw_blob_sha256'])
        rows.append({'stratum': stratum['name'], 'universe_size': len(paths), 'draws': actual})
    require(selected == [row['source_path'] for row in manifest['files']],
            'snapshot order differs from original draw')
    excerpts = ''
    for row in manifest['files']:
        data = (source / row['path']).read_bytes()
        require(len(data) == row['bytes'], 'snapshot size mismatch')
        require(hashlib.sha256(data).hexdigest() == row['sha256'], 'snapshot hash mismatch')
        require(row['sha256'] == native_blob_hashes[row['source_path']],
                'snapshot differs from native draw blob digest')
        blob = hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest()
        require(blob == row['git_blob'], 'snapshot Git blob mismatch')
        excerpts += '\n===== ' + row['source_path'] + ' =====\n'
        excerpts += ''.join(data.decode().splitlines(keepends=True)[:140])
    require(excerpts.encode() == (arithmetic / 'draw_excerpts.txt').read_bytes(),
            'native 140-line excerpts differ from pinned source snapshots')
    print(json.dumps({'status': 'PASS', 'seed': draw['seed'], 'source_pin': draw['pin'],
                      'sampled_files': len(selected), 'historical_code_executed': False,
                      'excerpt_rule': 'first 140 lines of each frozen source, preserving line endings',
                      'strata': rows}, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
