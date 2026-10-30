"""Check LF integrity and pinned inheritance; not a mathematical proof checker."""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
ARTIFACTS = ('README.md', 'PROOFS.md', 'ARITHMETIC_AND_SOURCE.md',
             'AUDIT.md', 'check.py', 'receipt.json', 'seal.py')
RUNTIMES = ('05-knowledge/results/universal_ramified_closure_20261030/check.py',
            '05-knowledge/results/regular_parent_symmetry_20261029/verify.py')
PINS = {
    'c830320f16': ('05-knowledge/results/universal_ramified_closure_20261030/PROOFS.md',
                  '05-knowledge/results/universal_ramified_closure_20261030/check.py',
                  '05-knowledge/results/regular_parent_symmetry_20261029/verify.py',
                  '05-knowledge/results/session78_fourth_power_aliquot_20261027/ALIQUOT.md',
                  '05-knowledge/results/session81_counterexample_families_20261028/ALIQUOT_OBSERVERS.md',
                  '05-knowledge/results/session82_decimal_aliquot_square_sum_20261028/ANCHOR_CHARACTERIZATION.md',
                  '05-knowledge/results/session108_sigma_symmetry_20261030/arithmetic/ARITHMETIC.md',
                  '05-knowledge/results/session111_scale_connections_20261030/ALIQUOT_MERGE.md'),
    'b51ef5cd73': ('05-knowledge/results/anchor_content_flat_quotients_20261030/TOURNAMENTS.md',),
    'fea3b07542': ('05-knowledge/results/session112_counting_interfaces_20261030/arithmetic/MODULAR_OBSERVERS.md',),
    'ecd0185d34': ('05-knowledge/results/session112_counting_interfaces_20261030/symmetry/PROOFS.md',
                  '05-knowledge/results/session112_counting_interfaces_20261030/counting/README.md')
}


def digest(blob):
    return hashlib.sha256(blob.replace(b'\r\n', b'\n')).hexdigest()


def git(*args):
    return subprocess.check_output(['git', *args], cwd=ROOT)


def snapshot():
    sources = []
    for ref, paths in PINS.items():
        commit = git('rev-parse', ref).decode().strip()
        for path in paths:
            source = git('show', commit+':'+path)
            if path in RUNTIMES and digest(source) != digest((ROOT/path).read_bytes()):
                raise RuntimeError('Live runtime differs from immutable source: '+path)
            sources.append({'commit': commit, 'path': path, 'lf_sha256': digest(source)})
    return {'status': 'Integrity only; same-agent proof review',
            'artifacts': {p: digest((HERE/p).read_bytes()) for p in ARTIFACTS},
            'source_pins': sources, 'live_runtime_inputs': list(RUNTIMES),
            'external_source': {
                'url': 'https://hexagonmath.org/pdf/2610.00600v1',
                'date_checked': '2026-10-30', 'pages_text_read': [1,2,3,4,5,6],
                'status': 'Preprint claim; no local Lean build or PDF visual inspection'}}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--record', action='store_true')
    args = parser.parse_args()
    current = snapshot()
    file = HERE/'manifest.json'
    if args.record:
        file.write_bytes((json.dumps(current, indent=2, sort_keys=True)+'\n').encode())
    elif current != json.loads(file.read_text(encoding='utf-8')):
        raise RuntimeError('Manifest mismatch')
    print(json.dumps({'match': True, 'artifacts': len(ARTIFACTS),
                      'source_pins': len(current['source_pins']), 'live_runtime_inputs': 2}))


if __name__ == '__main__':
    main()
