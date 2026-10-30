"""LF-normalized integrity and immutable inheritance, not a proof checker."""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
ARTIFACTS = ('README.md', 'PROOFS.md', 'CONNECTIONS.md', 'AUDIT.md',
             'check.py', 'receipt.json', 'seal.py')
RUNTIME = '05-knowledge/results/regular_parent_symmetry_20261029/verify.py'
PINS = {
    '7107d2f144': (
        RUNTIME,
        '05-knowledge/results/session78_fourth_power_aliquot_20261027/ALIQUOT.md',
        '05-knowledge/results/session81_counterexample_families_20261028/ALIQUOT_OBSERVERS.md',
        '05-knowledge/results/session108_sigma_symmetry_20261030/arithmetic/ARITHMETIC.md',
        '05-knowledge/results/session112_counting_interfaces_20261030/arithmetic/ALIQUOT_MODULAR.md',
        '05-knowledge/results/session112_counting_interfaces_20261030/arithmetic/POWER_OBSERVERS.md',
        '05-knowledge/results/overnight_rank_return_20261029/PERIODIC_COMPLEMENT.md',
        '05-knowledge/results/scaled_count_fft_carriers_20261030/THEORY.md',
        '05-knowledge/results/layer_lifts_cyclic_intersections_20261030/COLLATZ.md',
        '05-knowledge/results/pinned_ramified_counts_20261030/PROOFS.md'),
    '46677c32e4': (
        '05-knowledge/results/session113_symmetry_consumers_20261030/composite/COMPOSITE_PHASE.md',
        '05-knowledge/results/session113_symmetry_consumers_20261030/prime_ranking/PRIME_RANKING.md')
}


def digest(blob):
    return hashlib.sha256(blob.replace(b'\r\n', b'\n')).hexdigest()


def git(*args):
    return subprocess.check_output(['git', *args], cwd=ROOT)


def snapshot():
    pins = []
    for revision, paths in PINS.items():
        commit = git('rev-parse', revision).decode().strip()
        for path in paths:
            blob = git('show', commit+':'+path)
            if path == RUNTIME and digest(blob) != digest((ROOT/path).read_bytes()):
                raise RuntimeError('Changed inherited runtime: '+path)
            pins.append({'commit': commit, 'path': path, 'lf_sha256': digest(blob)})
    return {'status': 'Integrity only; same-agent proof review',
            'artifacts': {name: digest((HERE/name).read_bytes()) for name in ARTIFACTS},
            'source_pins': pins, 'live_runtime_inputs': [RUNTIME]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--record', action='store_true')
    args = parser.parse_args()
    current = snapshot()
    path = HERE/'manifest.json'
    if args.record:
        path.write_bytes((json.dumps(current, indent=2, sort_keys=True)+'\n').encode())
    elif current != json.loads(path.read_text(encoding='utf-8')):
        raise RuntimeError('Manifest mismatch')
    print(json.dumps({'match': True, 'artifacts': len(ARTIFACTS),
                      'source_pins': len(current['source_pins']), 'live_runtime_inputs': 1}))


if __name__ == '__main__':
    main()
