#!/usr/bin/env python3
"""Rebuild from pinned order-eight input. No network or external Python packages."""
import argparse
from pathlib import Path
import hashlib
import json
import os
import subprocess
import tempfile
import time

HERE = Path(__file__).resolve().parent

def need(value, label):
    if not value:
        raise RuntimeError(label)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-directory', required=True)
    args = parser.parse_args()
    output = Path(args.output_directory).resolve()
    output.mkdir(parents=True, exist_ok=True)
    need(output != HERE, 'replays use a separate output directory')
    pins = json.loads((HERE / 'SOURCE_PINS.json').read_text())
    parent_pin = next(x for x in pins['sources'] if x['archived_path'] == 'parents-eight.txt')
    need(hashlib.sha256((HERE / 'parents-eight.txt').read_bytes()).hexdigest() == parent_pin['sha256'], 'complete parent-eight source pin')
    environment = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
    receipts = []
    for phase, source, flags in [
        ('build-normal', 'build_nine.cpp', []), ('build-ndebug', 'build_nine.cpp', ['-DNDEBUG']),
        ('trace-normal', 'trace_ten.cpp', []), ('trace-full', 'trace_ten.cpp', ['-DNDEBUG', '-DFULL_REFINEMENT'])]:
        binary = output / phase
        subprocess.run(['g++', '-std=c++17', '-O3', '-Wall', '-Wextra', '-pedantic', *flags,
                        str(HERE / source), '-o', str(binary)], check=True, env=environment)
        if phase.startswith('build'):
            atlas = output / (phase + '-parents.txt')
            command = [str(binary), str(HERE / 'parents-eight.txt'), str(atlas), str(output / (phase + '.json'))]
        else:
            command = [str(binary), str(output / 'build-normal-parents.txt'), str(output / (phase + '.json')),
                       str(output / (phase + '-per-parent.tsv'))]
        start = time.monotonic()
        result = subprocess.run(command, capture_output=True, text=True, env=environment)
        (output / (phase + '.out')).write_text(result.stdout + result.stderr)
        need(result.returncode == 0, phase + ' process succeeded')
        receipts.append({'phase': phase, 'seconds': time.monotonic() - start})
        need(json.loads((output / (phase + '.json')).read_text()) == json.loads((HERE / (phase + '.json')).read_text()), phase + ' exact receipt replay')
    expected = (HERE / 'parents-nine.txt').read_bytes()
    for phase in ('build-normal', 'build-ndebug'):
        need((output / (phase + '-parents.txt')).read_bytes() == expected, phase + ' exact atlas replay')
    detail = (HERE / 'trace-normal-per-parent.tsv').read_bytes()
    for phase in ('trace-normal', 'trace-full'):
        need((output / (phase + '-per-parent.tsv')).read_bytes() == detail, phase + ' exact per-parent replay')
    receipt = {'status': 'PASS', 'runs': receipts,
               'parent_atlas_sha256': hashlib.sha256(expected).hexdigest(),
               'per_parent_trace_sha256': hashlib.sha256(detail).hexdigest()}
    (output / 'replay-receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, sort_keys=True))

if __name__ == '__main__':
    main()
