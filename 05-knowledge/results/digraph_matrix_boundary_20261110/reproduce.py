#!/usr/bin/env python3
"""Verify frozen evidence; --rerun freshly runs each normal/optimized checker.

--write-manifest records the completed packet after all evidence checks.
Default verification does not rerun the small censuses; it freshly checks
the displayed matrix witnesses, script/source hashes and local links.
"""
from hashlib import sha256
from itertools import permutations
from pathlib import Path
import argparse
import json
import re
import subprocess
import sys
from urllib.parse import unquote

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
JOBS = [
    ('construction/check.py', 'construction/check.json', 'construction/check.optimized.json'),
    ('history/compact_audit.py', 'history/compact_audit.json', 'history/compact_audit.optimized.json'),
    ('completion/check.py', 'completion/check.json', 'completion/check.optimized.json'),
    ('enumeration/check.py', 'enumeration/normal.json', 'enumeration/optimized.json'),
    ('history/final_audit_check.py', 'history/final_audit_check.json', 'history/final_audit_check.optimized.json'),
]


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def paths(matrix):
    n = len(matrix)
    return [p for p in permutations(range(n))
            if all(matrix[u][v] == '1' for u, v in zip(p, p[1:]))]


def literal_controls():
    samples = [
        (7, ['0111', '0011', '1001', '1000']),
        (21, ['01111', '00110', '00011', '10001', '11000']),
    ]
    for h, matrix in samples:
        n = len(matrix)
        need(len(paths(matrix)) == h, 'ordinary displayed witness')
        need(all(matrix[i][i] == '0' for i in range(n)), 'loopless witness')
        need(all(matrix[i][j] == '1' or matrix[j][i] == '1'
                 for i in range(n) for j in range(i + 1, n)), 'no00 witness')
    raw = ['01110', '00111', '10011', '00001', '10000']
    decoded = [''.join('1' if i != j and raw[i][j] == '1' and raw[j][i] == '0' else '0'
                       for j in range(5)) for i in range(5)]
    need(len(paths(raw)) == 11 and len(paths(decoded)) == 7, 'raw/decoded control')
    return dict(ordinary_witnesses=[7, 21], raw_cancelled_pair=[11, 7])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rerun', action='store_true')
    parser.add_argument('--write-manifest', action='store_true')
    args = parser.parse_args()
    checked_inputs = set()
    for script, normal, optimized in JOBS:
        data = (HERE / normal).read_bytes()
        need(data == (HERE / optimized).read_bytes(), 'normal/-O mismatch: ' + script)
        obj = json.loads(data)
        need(obj.get('status') == 'PASS', 'nonpassing evidence: ' + normal)
        self_pin = obj.get('script_sha256', obj.get('checker_sha256'))
        if self_pin:
            need(digest(HERE / script) == self_pin, 'stale checker: ' + script)
        for rel, expected in obj.get('source_sha256', {}).items():
            matches = {p.resolve() for p in (ROOT / rel, HERE / rel) if p.is_file()}
            need(len(matches) == 1, 'unique pinned source path: ' + rel)
            source = matches.pop()
            need(digest(source) == expected, 'stale pinned source: ' + rel)
            checked_inputs.add(str(source.relative_to(ROOT)))
        if args.rerun:
            for opt in (False, True):
                command = [sys.executable, '-B'] + (['-O'] if opt else []) + [str(HERE / script)]
                fresh = subprocess.run(command, cwd=ROOT, capture_output=True, check=True)
                need(fresh.stdout == data, 'fresh result differs: ' + script)
            print('PASS fresh normal/-O ' + script, flush=True)
    link_count = 0
    for doc in HERE.rglob('*.md'):
        for raw in re.findall(r'\]\(([^)]+)\)', doc.read_text()):
            target = unquote(raw.strip('<>').split('#', 1)[0])
            if not target or re.match(r'\w+://', target):
                continue
            need((doc.parent / target).exists(), 'broken link in ' + str(doc) + ': ' + target)
            link_count += 1
    controls = literal_controls()
    files = {str(p.relative_to(HERE)): digest(p) for p in sorted(HERE.rglob('*'))
             if p.is_file() and '__pycache__' not in p.parts and p.name != 'manifest.json'}
    record = dict(scope='Packet artifact hashes; excludes this manifest and Python caches', files=files)
    manifest = HERE / 'manifest.json'
    if args.write_manifest:
        manifest.write_text(json.dumps(record, indent=2, sort_keys=True) + '\n')
    else:
        need(json.loads(manifest.read_text()) == record, 'packet manifest differs')
    print(json.dumps(dict(status='PASS', mode='fresh checkers' if args.rerun else 'saved evidence and literal witnesses',
        normal_optimized_pairs=len(JOBS), distinct_pinned_inputs=len(checked_inputs),
        local_links=link_count, artifacts=len(files), fresh_literal_controls=controls), sort_keys=True))


if __name__ == '__main__':
    main()
