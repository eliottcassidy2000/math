#!/usr/bin/env python3
"""Replay the three immutable producer archives, separately from our checks."""
from pathlib import Path
from hashlib import sha256
import argparse
import json
import subprocess
import sys

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--write', action='store_true'); args = parser.parse_args()
    answer = {}
    for name in ('hole-theory', 'disjoint-sections', 'nonreturn-repairs'):
        folder = HERE / 'inputs' / name
        results = []
        for flags in (['-I', '-B'], ['-I', '-O', '-B']):
            p = subprocess.run([sys.executable, *flags, str(folder / 'check.py')], cwd=HERE, capture_output=True, text=True, check=True)
            results.append({'stdout_sha256': sha256(p.stdout.encode()).hexdigest(), 'stderr': p.stderr, 'exit_code': p.returncode})
        if results[0] != results[1]:
            raise RuntimeError('Producer normal and optimized outputs differ: ' + name)
        answer[name] = {'manifest_sha256': sha256((folder / 'MANIFEST.json').read_bytes()).hexdigest(), 'normal_and_optimized': results[0]}
    payload = json.dumps(answer, indent=2, sort_keys=True) + '\n'
    if args.write:
        (HERE / 'PRODUCER_REPLAY.json').write_text(payload)
    elif (HERE / 'PRODUCER_REPLAY.json').read_text() != payload:
        raise RuntimeError('Producer replay differs from frozen receipt')
    print('PASS all three producer archives, normal and optimized')


if __name__ == '__main__':
    main()
