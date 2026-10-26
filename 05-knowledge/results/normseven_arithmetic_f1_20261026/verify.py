"""Verify pinned bytes and optionally replay all independent exact paths."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

sys.stdout.reconfigure(newline="\n")
HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--replay', action='store_true')
    args = parser.parse_args()
    manifest = json.loads((HERE/'MANIFEST.json').read_text())
    actual = {p.relative_to(HERE).as_posix() for p in HERE.rglob('*')
              if p.is_file() and p.name != 'MANIFEST.json' and '__pycache__' not in p.parts}
    need(actual == set(manifest), 'manifest exact file set')
    for relative, digest in manifest.items():
        need(hashlib.sha256((HERE/relative).read_bytes()).hexdigest() == digest, 'hash: '+relative)
    sources = json.loads((HERE/'SOURCES.json').read_text())
    for relative, digest in sources['files'].items():
        need(hashlib.sha256((REPO/relative).read_bytes()).hexdigest() == digest, 'source changed: '+relative)
    jobs = [
        ('boundary/check.py','boundary/normal.json','boundary/optimized.json'),
        ('boundary/audit_gale.py','boundary/audit_gale.json','boundary/audit_gale_optimized.json'),
        ('gale/check.py','gale/normal.json','gale/optimized.json'),
        ('gale/audit_boundary.py','gale/boundary-audit.normal.json','gale/boundary-audit.optimized.json'),
        ('f1/check.py','f1/normal.json','f1/optimized.json'),
        ('reduction/check_reduction.py','reduction/normal.json','reduction/optimized.json'),
        ('reduction/audit_f1.py','reduction/audit_normal.json','reduction/audit_optimized.json'),
    ]
    for script,normal,optimized in jobs:
        expected=(HERE/normal).read_bytes()
        need(expected == (HERE/optimized).read_bytes(), 'paired bytes: '+script)
        if args.replay:
            for mode in ([],['-O']):
                result=subprocess.check_output([sys.executable,'-B']+mode+[str(HERE/script)],cwd=REPO)
                need(result == expected,'replay bytes: '+script+str(mode))
            print('PASS replay '+script,flush=True)
    if args.replay:
        expected=(HERE/'f1/parent_audit.json').read_bytes()
        for mode in ([],['-O']):
            subprocess.check_call([sys.executable,'-B']+mode+[str(HERE/'f1/audit_parent.py')],cwd=REPO)
            need((HERE/'f1/parent_audit.json').read_bytes()==expected,'parent audit replay')
    print(f'PASS: {len(manifest)} packet hashes, {len(sources["files"])} source pins, 7 paired receipts; replay={args.replay}')


if __name__ == '__main__':
    main()
