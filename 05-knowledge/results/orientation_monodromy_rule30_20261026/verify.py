#!/usr/bin/env python3
"""Closed inventory, immutable source pins, and optional exact replay."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def digest(data):
    return hashlib.sha256(data.replace(b'\r\n', b'\n')).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--replay', action='store_true')
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    repo = here.parents[2]
    manifest = json.loads((here/'MANIFEST.json').read_text())
    actual = {p.name for p in here.iterdir() if p.is_file() and p.name != 'MANIFEST.json'}
    need(actual == set(manifest['packet']), 'closed inventory')
    for name, pin in manifest['packet'].items():
        need(digest((here/name).read_bytes()) == pin, ('packet hash', name))
    for name, pin in manifest['sources'].items():
        data = subprocess.check_output(['git', 'show', pin['revision']+':'+name], cwd=repo)
        need(digest(data) == pin['sha256'], ('immutable source hash', name))
    expected = (here/'normal.json').read_bytes()
    need(expected == (here/'optimized.json').read_bytes(), 'normal/optimized equality')
    runs = 0
    if args.replay:
        for flags in ([], ['-O']):
            data = subprocess.check_output([sys.executable, *flags, '-B', str(here/'check.py')], cwd=repo)
            need(data == expected, ('fresh replay', flags))
            runs += 1
    print(json.dumps({'status':'PASS', 'packet_files':len(actual),
                      'source_pins':len(manifest['sources']), 'fresh_runs':runs,
                      'proof_review':'internal; no independent agent or formal audit'}, sort_keys=True))


if __name__ == '__main__':
    main()
