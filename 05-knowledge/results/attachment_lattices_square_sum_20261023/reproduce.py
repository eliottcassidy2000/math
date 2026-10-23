"""Verify the frozen packet and replay every checker normally and under -O.

--write-manifest is an explicit maintenance operation, not a test result.
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
CHECKERS = (
    'check_kernels_moves.py',
    'check_joint_alignment.py',
    'lattice/check.py',
    'lattice/audit_non_normal_kernel.py',
    'complement/check.py',
    'complement/cut_response.py',
    'complement/check_flow.py',
    'complement/audit_lattice.py',
    'square_sum/check.py',
    'square_sum/audit_complement.py',
    'square_sum/check_shared_clocks.py',
    'square_sum/audit_shared_clocks_independent.py',
)


def inventory():
    paths = [p for p in HERE.rglob('*') if p.is_file()
             and p.suffix in {'.py', '.md', '.json'} and p != MANIFEST]
    records = []
    for path in sorted(paths, key=lambda p: p.relative_to(ROOT).as_posix()):
        data = path.read_bytes().replace(b'\r\n', b'\n')
        data.decode('utf-8')
        records.append(dict(path=path.relative_to(ROOT).as_posix(),
                            encoding='utf8-lf',
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
        print('Manifest written; this operation has NOT run the checks.')
        return
    if json.loads(MANIFEST.read_text(encoding='utf-8')) != actual:
        raise ValueError('Packet differs from manifest; inspect changes before freezing.')
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
