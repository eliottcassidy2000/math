"""Verify frozen inputs, then replay every checker normally and under -O.

No assertion in this wrapper is removed by optimization. --write-manifest
is an explicit maintenance operation and does not run or certify checks.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
MANIFEST = HERE / 'manifest.json'
EXTERNAL = (
    '05-knowledge/results/endpoint_overlap_ports_20261023/receipt.json.gz',
    '05-knowledge/results/odd_cell_gluing_20261023/check_fibonacci.py',
)
CHECKERS = (
    'check_networks.py',
    'check_residual.py',
    'five/check.py',
    'five/check_prime_refinement.py',
    'five/audit_networks.py',
    'five/audit_canonical.py',
    'seven/check.py',
    'seven/canonical.py',
    'seven/audit_network.py',
    'fibonacci/check.py',
)


def inventory():
    paths = [p for p in HERE.rglob('*') if p.is_file()
             and p.suffix in {'.py', '.md', '.json'} and p != MANIFEST]
    paths.extend(ROOT / name for name in EXTERNAL)
    records = []
    for path in sorted(paths, key=lambda p: p.relative_to(ROOT).as_posix()):
        binary = path.suffix == '.gz'
        data = path.read_bytes()
        if not binary:
            data = data.replace(b'\r\n', b'\n')
            data.decode('utf-8')
        records.append(dict(path=path.relative_to(ROOT).as_posix(),
                            encoding='raw' if binary else 'utf8-lf',
                            sha256=hashlib.sha256(data).hexdigest()))
    return dict(format=1, files=records, checkers=list(CHECKERS))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-manifest', action='store_true')
    parser.add_argument('--manifest-only', action='store_true')
    args = parser.parse_args()
    actual = inventory()
    if args.write_manifest:
        MANIFEST.write_text(json.dumps(actual, indent=2) + '\n', encoding='utf-8')
        print('Manifest written; checks have NOT been run by this operation.')
        return
    if json.loads(MANIFEST.read_text(encoding='utf-8')) != actual:
        raise ValueError('Manifest differs: inspect source/input changes before freezing again.')
    print(json.dumps(dict(manifest='PASS', files=len(actual['files']))), flush=True)
    if args.manifest_only:
        return
    for optimized in (False, True):
        for relative in CHECKERS:
            command = [sys.executable, '-B'] + (['-O'] if optimized else [])
            command.append(str(HERE / relative))
            print(json.dumps(dict(checker=relative, optimized=optimized)), flush=True)
            subprocess.run(command, cwd=ROOT, check=True)
    print(json.dumps(dict(status='PASS', processes=2 * len(CHECKERS))), flush=True)


if __name__ == '__main__':
    main()
