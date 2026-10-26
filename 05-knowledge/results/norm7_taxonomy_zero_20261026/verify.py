#!/usr/bin/env python3
"""Check immutable packet/source bytes, optionally replay independent paths."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def demand(condition, message):
    if not condition:
        raise RuntimeError(message)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--replay', action='store_true')
    parser.add_argument('--pari', action='store_true')
    args = parser.parse_args()
    manifest = json.loads((HERE / 'MANIFEST.json').read_text())
    actual = {p.relative_to(HERE).as_posix() for p in HERE.rglob('*')
              if p.is_file() and '__pycache__' not in p.parts
              and p.relative_to(HERE).as_posix() != 'MANIFEST.json'}
    demand(actual == set(manifest['files']), 'closed file manifest mismatch')
    for name, expected in manifest['files'].items():
        demand(digest(HERE / name) == expected, 'file hash: ' + name)
    for name, expected in manifest['source_pins'].items():
        demand(digest(ROOT / name) == expected, 'source hash: ' + name)
    lanes = [
        ('arithmetic/check_quadratic_atlas.py', 'arithmetic/normal.json', 'arithmetic/optimized.json'),
        ('certificates/audit_quadratic_atlas.py', 'certificates/quadratic_audit.normal.json', 'certificates/quadratic_audit.optimized.json'),
        ('motives/check_extensions.py', 'motives/normal.json', 'motives/optimized.json'),
        ('arithmetic/audit_extension_independent.py', 'arithmetic/extension_audit.normal.json', 'arithmetic/extension_audit.optimized.json'),
        ('taxonomy/check_taxonomy.py', 'taxonomy/receipt.json', None),
        ('certificates/check_certificates.py', 'certificates/results.json', 'certificates/results.optimized.json'),
    ]
    for script, normal, optimized in lanes:
        expected = (HERE / normal).read_bytes()
        if optimized:
            demand(expected == (HERE / optimized).read_bytes(), 'mode pair: ' + script)
        if args.replay:
            for flags in ([], ['-O']):
                result = subprocess.run([sys.executable, '-B', *flags, str(HERE / script)],
                                        cwd=(HERE / script).parent, capture_output=True, check=True)
                demand(result.stdout == expected, 'replay: ' + script + ' ' + str(flags))
            print('PASS normal/optimized ' + script, flush=True)
    if args.pari:
        result = subprocess.run(['gp', '-q', 'replay_modular_rows.gp'],
                                cwd=HERE / 'certificates', capture_output=True, check=True)
        demand(result.stdout == (HERE / 'certificates/modular_rows.out').read_bytes(), 'PARI replay')
        print('PASS PARI modular coefficient replay', flush=True)
    print(json.dumps({'status': 'PASS', 'packet_files': len(actual),
                      'source_pins': len(manifest['source_pins']),
                      'replay': args.replay, 'pari': args.pari}, sort_keys=True))


if __name__ == '__main__':
    main()
