"""Reproduce the full independent finite audit; C++17 + Python standard library.

No external representative lists are needed. Optional --producer-input adds a
bijection check against the submitted reps1.txt,...,reps9.txt without publishing
those scratch lists. All products remain under this audit directory.
"""
import argparse
import hashlib
import json
import platform
from pathlib import Path
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
SOURCES = ['independent_census.cpp', 'check_phase.py', 'check_polynomial.py',
           'check_automorphisms.py', 'run_audit.py', 'phase_bank.json',
           'small_representatives.json', 'producer_snapshot/cut_polynomial.py',
           'producer_snapshot/full_binary.py']

def need(condition, message):
    if not condition:
        raise ValueError(message)

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--producer-input', type=Path)
    args = parser.parse_args()
    build = HERE / '_build'
    build.mkdir(exist_ok=True)
    runs = []
    compiler = subprocess.check_output(['g++', '--version'], text=True).splitlines()[0]
    for mode, flags, pyflags in [('normal', ['-O2'], []), ('optimized', ['-O3', '-DNDEBUG'], ['-O'])]:
        binary = build / ('census.' + mode)
        command = ['g++', '-std=c++17', '-Wall', '-Wextra', '-Wpedantic', *flags,
                   str(HERE / 'independent_census.cpp'), '-o', str(binary)]
        subprocess.run(command, check=True, cwd=build)
        census_command = [str(binary)]
        if args.producer_input:
            census_command += ['--producer-input', str(args.producer_input.resolve())]
        commands = [
            ('census', census_command),
            ('phase', [sys.executable, *pyflags, str(HERE / 'check_phase.py')]),
            ('polynomial', [sys.executable, *pyflags, str(HERE / 'check_polynomial.py'), '--binary', str(binary)]),
            ('automorphisms', [sys.executable, *pyflags, str(HERE / 'check_automorphisms.py'), '--census', str(HERE / ('census.' + mode + '.json'))])]
        for name, command in commands:
            start = time.monotonic()
            output = subprocess.check_output(command, cwd=build)
            json.loads(output)
            destination = HERE / (name + '.' + mode + '.json')
            destination.write_bytes(output)
            runs.append({'check': name, 'mode': mode, 'exit': 0,
                         'seconds': round(time.monotonic() - start, 3), 'sha256': sha(destination)})
            print(name + '.' + mode + ': PASS', flush=True)
    for name in ['census', 'phase', 'polynomial', 'automorphisms']:
        need((HERE / (name + '.normal.json')).read_bytes() == (HERE / (name + '.optimized.json')).read_bytes(),
             'normal/optimized mismatch: ' + name)
    # This replay has no producer-input argument and runs outside the package's
    # root. It verifies the standalone generation path through order 7.
    smoke = subprocess.check_output([str(build / 'census.optimized'), '--max-order', '7'], cwd=build)
    smoke_data = json.loads(smoke)
    full = json.loads((HERE / 'census.normal.json').read_text())
    need(smoke_data['orders'] == full['orders'][:7], 'standalone generation replay mismatch')
    (HERE / 'portable_smoke.json').write_bytes(smoke)
    receipt = {
        'status': 'PASS', 'scope': 'all tournaments of orders 1..9; parent H>=9; arbitrary cuts; balanced depth2',
        'compiler': compiler, 'python': platform.python_version(),
        'normal_cpp_flags': ['-O2'], 'optimized_cpp_flags': ['-O3', '-DNDEBUG'],
        'normal_optimized_byte_identical': True,
        'producer_representative_bijection_checked': bool(args.producer_input),
        'standalone_no_representative_input_replay_max_order': 7,
        'all_full_runs_cwd': 'audit-local _build (not the source directory)',
        'direct_child_DP_checks': sum(row['direct_child_dp_checks'] for row in full['orders']),
        'qualifying_child_cuts': sum(row['qualifying_children'] for row in full['orders']),
        'grandchild_cuts_checked_by_direct_DP': sum(row['grandchild_cuts'] for row in full['orders']),
        'sources': {name: {'sha256': sha(HERE / name), 'bytes': (HERE / name).stat().st_size} for name in SOURCES},
        'runs': runs,
        'portable_smoke_sha256': sha(HERE / 'portable_smoke.json')}
    (HERE / 'AUDIT_RECEIPT.json').write_text(json.dumps(receipt, sort_keys=True, indent=2) + '\n')
    print(json.dumps({key: receipt[key] for key in ['status', 'direct_child_DP_checks', 'qualifying_child_cuts', 'grandchild_cuts_checked_by_direct_DP', 'normal_optimized_byte_identical']}, sort_keys=True))

if __name__ == '__main__':
    main()
