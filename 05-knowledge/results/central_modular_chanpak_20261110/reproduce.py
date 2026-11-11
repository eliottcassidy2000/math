#!/usr/bin/env python3
"""Replay explicit finite universes, source pins, links and frozen artifacts."""
from hashlib import sha256
from pathlib import Path
import argparse
import json
import re
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BASE = '07dc31fae3'
JOBS = (
    ('central/check.py', 'central/evidence.json', 'central/evidence.optimized.json'),
    ('information/check.py', 'information/check.json', 'information/check.optimized.json'),
    ('compiler/check.py', 'compiler/check.json', 'compiler/check.optimized.json'),
    ('gauge/check.py', 'gauge/check.json', 'gauge/check.optimized.json'),
    ('central/audit_information.py', 'central/audit_information_evidence.json',
     'central/audit_information_evidence.optimized.json'),
    ('information/audit_compiler.py', 'information/audit_compiler.json',
     'information/audit_compiler.optimized.json'),
    ('compiler/incoming_height_check.py', 'compiler/incoming_height_check.json',
     'compiler/incoming_height_check.optimized.json'),
)

def require(ok, detail):
    if not ok:
        raise RuntimeError(detail)

def digest(path):
    return sha256(path.read_bytes()).hexdigest()

def dump(path, data):
    path.write_text(json.dumps(data, indent=2, sort_keys=True)+'\n')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--rerun', action='store_true')
    parser.add_argument('--refresh', action='store_true',
                        help='Explicitly regenerate receipts, source pins and manifest')
    args = parser.parse_args()
    source_paths = set()
    for script, normal, optimized in JOBS:
        for opt, output in ((False, normal), (True, optimized)):
            if args.rerun or args.refresh:
                command = [sys.executable, '-B'] + (['-O'] if opt else []) + [str(HERE/script)]
                run = subprocess.run(command, cwd=ROOT, check=True, capture_output=True)
                data = json.loads(run.stdout)
                if args.refresh:
                    dump(HERE/output, data)
                else:
                    require(data == json.loads((HERE/output).read_text()), ('fresh receipt', script, opt))
        require((HERE/normal).read_bytes() == (HERE/optimized).read_bytes(), ('optimized mismatch', script))
        data = json.loads((HERE/normal).read_text())
        require(str(data['status']).startswith(('PASS', 'FINITE-EXACT')), ('receipt status', script))
        if 'script_sha256' in data:
            require(data['script_sha256'] == digest(HERE/script), ('script changed', script))
        for source, expected in data.get('source_sha256', {}).items():
            path = ROOT/source
            require(digest(path) == expected, ('receipt input changed', source))
            if not path.is_relative_to(HERE):
                source_paths.add(source)
    links = 0
    for md in HERE.rglob('*.md'):
        for link in re.findall(r'\]\(([^\s)]+)\)', md.read_text()):
            if '://' in link or link.startswith('#'):
                continue
            path = (md.parent/link.split('#', 1)[0]).resolve()
            require(path.is_file(), ('broken local link', str(md.relative_to(HERE)), link))
            links += 1
            if path.is_relative_to(ROOT) and not path.is_relative_to(HERE):
                source_paths.add(str(path.relative_to(ROOT)))
    if args.refresh:
        dump(HERE/'SOURCE_PINS.json', {
            'source_HEAD': BASE,
            'source_sha256': {p:digest(ROOT/p) for p in sorted(source_paths)},
            'literature': [{
                'url':'https://www.math.ucla.edu/~pak/papers/DominoPermanent16.pdf',
                'location':'section 6.5, Conjecture 6.2, PDF page 17',
                'checked':'2026-11-10',
                'claim':'odd tournament path-count range except 7 and 21; credits bof/Royle2016'
            }],
        })
    pins = json.loads((HERE/'SOURCE_PINS.json').read_text())
    for name, expected in pins['source_sha256'].items():
        path = ROOT/name
        if path.is_file() and digest(path) == expected:
            continue
        # Mutable source notes remain identifiable at the pinned inherited commit.
        raw = subprocess.run(['git','show',pins['source_HEAD']+':'+name],
                             cwd=ROOT, check=True, capture_output=True).stdout
        require(sha256(raw).hexdigest() == expected, ('source pin', name))
    artifacts = {str(p.relative_to(HERE)):digest(p) for p in sorted(HERE.rglob('*'))
                 if p.is_file() and p.name != 'MANIFEST.json' and '__pycache__' not in p.parts}
    if args.refresh:
        dump(HERE/'MANIFEST.json', artifacts)
    require(json.loads((HERE/'MANIFEST.json').read_text()) == artifacts, 'artifact manifest mismatch')
    print(json.dumps({'status':'PASS', 'checker_pairs':len(JOBS),
                      'source_pins':len(pins['source_sha256']), 'local_links':links,
                      'artifacts':len(artifacts), 'fresh_replay':args.rerun or args.refresh},
                     indent=2, sort_keys=True))

if __name__ == '__main__':
    main()
