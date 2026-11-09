#!/usr/bin/env python3
"""Replay Session235 without changing packet artifacts.

Requirements: Python 3 (standard library) and a C++17 compiler named by CXX
(default g++). Run anywhere: python /path/to/audit/replay_all.py.
Optional --receipt PATH writes only that explicitly requested receipt. All
compiler products, scratch tables, and regenerated outputs live in an automatic
temporary directory. The inherited Session234 generator and results6.csv are
copied into the sibling layout expected by the existing scripts.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile

PACKET = Path(__file__).resolve().parent.parent
INHERITED = PACKET.parent / 'session234_digraphs_20261109'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--receipt', type=Path, help='Optional destination for deterministic JSON receipt')
    args = parser.parse_args()
    compiler = shlex.split(os.environ.get('CXX', 'g++'))
    if not compiler or shutil.which(compiler[0]) is None:
        raise RuntimeError('C++ compiler unavailable; set CXX to a C++17 compiler')
    checked = []
    steps = []
    # Snapshot all reference files before running. No source output is overwritten.
    references = {p.relative_to(PACKET).as_posix(): p.read_bytes()
                  for p in PACKET.rglob('*') if p.is_file()
                  and '__pycache__' not in p.parts and p.suffix != '.bin'
                  and p.name != 'replay_all.json'}
    def compare_bytes(relative, actual):
        expected = references[relative]
        if actual != expected:
            raise AssertionError('Byte mismatch: ' + relative)
        checked.append({'file': relative, 'comparison': 'byte-for-byte',
                        'sha256': hashlib.sha256(expected).hexdigest()})
    def compare_json(relative, actual):
        expected = json.loads(references[relative])
        if actual != expected:
            raise AssertionError('JSON mismatch: ' + relative)
        canonical = json.dumps(expected, sort_keys=True, separators=(',', ':')).encode()
        checked.append({'file': relative, 'comparison': 'JSON value equality',
                        'canonical_sha256': hashlib.sha256(canonical).hexdigest()})
    with tempfile.TemporaryDirectory(prefix='session235-replay-') as temporary:
        root = Path(temporary)
        work = root / PACKET.name
        old = root / INHERITED.name / 'census'
        old.mkdir(parents=True)
        for name in ('generator.py', 'results6.csv'):
            shutil.copy2(INHERITED / 'census' / name, old / name)
        for directory in ('census', 'models', 'theory', 'audit'):
            (work / directory).mkdir(parents=True)
            for source in (PACKET / directory).iterdir():
                if source.suffix in ('.py', '.cpp') and source.name != 'replay_all.py':
                    shutil.copy2(source, work / directory / source.name)
        environment = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
        def run(command, cwd=work):
            result = subprocess.run(command, cwd=cwd, env=environment,
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                    check=False)
            if result.returncode:
                # Failure diagnostics may contain temporary paths, but are never
                # put into successful receipts or hashed as result data.
                raise RuntimeError('Command failed: ' + ' '.join(map(str, command))
                                   + '\n' + result.stderr.decode(errors='replace')
                                   + '\n' + result.stdout.decode(errors='replace'))
            return result.stdout
        run(compiler + ['-O3', '-std=c++17', 'census/census.cpp', '-o', str(root / 'census')])
        producer_output = b''
        for n in range(1, 7):
            producer_output += run([str(root / 'census'), str(n)], work / 'census')
            for prefix in ('hist', 'joint', 'classes', 'missing'):
                relative = f'census/{prefix}{n}.csv'
                compare_bytes(relative, (work / relative).read_bytes())
        compare_bytes('census/run.txt', producer_output)
        compare_bytes('audit/producer_replay.txt', producer_output)
        steps.append('C++ census: regenerate and byte-compare all 24 CSVs, orders 1..6')
        run(compiler + ['-O3', '-std=c++17', 'audit/scatter_compare.cpp', '-o', str(root / 'scatter')])
        scatter_output = run([str(root / 'scatter'), str(work / 'census')])
        compare_bytes('audit/scatter_compare.txt', scatter_output)
        steps.append('Independent direct scatter: compare every oriented table entry through order 6')
        run([sys.executable, 'verify.py'], work / 'census')
        for name in ('validation.json', 'witnesses.json'):
            compare_json('census/' + name, json.loads((work / 'census' / name).read_text()))
        steps.append('Census verifier: exhaustive endpoint DP through order 5; 3398 order-6 cases; all target paths; inherited semicomplete marginal')
        conditional = {}
        for n in range(1, 7):
            rows = list(csv.DictReader((work / 'census' / f'joint{n}.csv').open()))
            data = {}
            for target in (7, 21):
                data[f'oriented_H{target}_to_semicomplete'] = {
                    row['H_semicomplete']: int(row['count']) for row in rows
                    if int(row['H_oriented']) == target}
                data[f'semicomplete_H{target}_to_oriented'] = {
                    row['H_oriented']: int(row['count']) for row in rows
                    if int(row['H_semicomplete']) == target}
            conditional[str(n)] = data
        compare_json('census/target_decode_histograms.json', conditional)
        steps.append('Conditional target histograms: reconstruct from regenerated joint tables')
        models = run([sys.executable, 'models/verify_models.py'])
        compare_json('models/validation.json', json.loads(models))
        steps.append('Pair-model verifier: all four-state and ternary arrays through order 4')
        run([sys.executable, 'theory/check_bridges.py'])
        compare_json('theory/checks.json', json.loads((work / 'theory/checks.json').read_text()))
        steps.append('Theory verifier: exact parity/forest and bridge checks through order 4; sparse strong-core controls')
        run([sys.executable, 'census/core_transport.py'])
        compare_json('census/core_transport.json', json.loads((work / 'census/core_transport.json').read_text()))
        steps.append('Core transport: 37 inherited cores, paths and automorphisms; grammar coefficients through order 10')
        for stem in ('independent', 'theory_controls', 'core_transport_check'):
            run([sys.executable, 'audit/' + stem + '.py'])
            compare_json('audit/' + stem + '.json', json.loads((work / 'audit' / (stem + '.json')).read_text()))
        steps.append('Independent audit controls: one-deletion scan through 6; 2460 hub-cycle cases; 36 bouquets; 4165 parity cases; core/GF checks through order 12')
    for relative, original in references.items():
        if (PACKET / relative).read_bytes() != original:
            raise AssertionError('Reference changed during replay: ' + relative)
    receipt = {
        'status': 'PASS', 'packet': PACKET.name,
        'source_outputs_unchanged': True,
        'temporary_layout_cleaned': True,
        'scope': steps,
        'limitations': [
            'All-order claims use the stated proofs and inherited complete Session234 core classification, not extrapolation.',
            'Full second-algorithm table comparison is oriented; semicomplete verification uses the inherited full marginal plus endpoint checks.',
            'The inherited Session234 generator and order-six marginal are inputs, not independently regenerated by this runner.'
        ],
        'verified_artifacts': sorted(checked, key=lambda item: item['file']),
        'inherited_inputs': [
            {'file': 'session234_digraphs_20261109/census/' + name,
             'sha256': hashlib.sha256((INHERITED / 'census' / name).read_bytes()).hexdigest()}
            for name in ('generator.py', 'results6.csv')
        ]
    }
    output = json.dumps(receipt, indent=2) + '\n'
    if args.receipt is not None:
        destination = args.receipt.resolve()
        if destination.is_relative_to(PACKET) and destination != PACKET / 'audit/replay_all.json':
            raise ValueError('Inside this packet, receipt may only be written as audit/replay_all.json')
        destination.write_text(output)
    print(output, end='')


if __name__ == '__main__':
    main()
