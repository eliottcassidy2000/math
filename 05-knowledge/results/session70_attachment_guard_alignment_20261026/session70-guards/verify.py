#!/usr/bin/env python3
"""Exact file universe, immutable receipt replays, and active negative controls."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent

def need(ok, why):
    if not ok:
        raise RuntimeError(why)

def main():
    manifest = json.loads((ROOT / 'MANIFEST.json').read_text())
    expected = {x['name'] for x in manifest['files']} | {'MANIFEST.json'}
    actual = {p.name for p in ROOT.iterdir() if p.is_file()}
    need(actual == expected, 'missing or extra packet files')
    for record in manifest['files']:
        data = (ROOT / record['name']).read_bytes()
        need(len(data) == record['bytes'], 'byte count: ' + record['name'])
        need(hashlib.sha256(data).hexdigest() == record['sha256'], 'hash: ' + record['name'])
    receipt = (ROOT / 'normal.json').read_bytes()
    need(receipt == (ROOT / 'optimized.json').read_bytes(), 'frozen normal/optimized receipts differ')
    for optimized in (False, True):
        command = [sys.executable, '-I', '-B'] + (['-O'] if optimized else []) + [str(ROOT / 'check.py')]
        result = subprocess.run(command, capture_output=True, check=False)
        need(result.returncode == 0, 'replay failed: ' + result.stderr.decode())
        need(result.stdout == receipt, 'replay receipt changed')
        hostile = subprocess.run(command + ['--mutate'], capture_output=True, check=False)
        need(hostile.returncode != 0, 'corrupted decoded arc was accepted')
        need(b'decoded relation differs from original' in hostile.stderr, 'wrong negative-control failure')
    print(json.dumps({'status': 'PASS', 'exact_files': len(expected), 'normal_replay': 'PASS',
                      'optimized_replay': 'PASS', 'decoded_arc_corruption_normal': 'REJECTED',
                      'decoded_arc_corruption_optimized': 'REJECTED'}, sort_keys=True, indent=2))

if __name__ == '__main__':
    main()
