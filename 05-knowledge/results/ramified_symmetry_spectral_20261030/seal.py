"""Immutable source pins, live runtime inputs, and LF artifact integrity."""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
ARTIFACTS = ('README.md', 'PROOFS.md', 'AUDIT.md', 'check.py', 'certificate.json',
             'parent_group.py', 'parent_group.json', 'seal.py')
RUNTIME = (
    '05-knowledge/results/regular_parent_symmetry_20261029/verify.py',
    '05-knowledge/results/pinned_discrepancy_reset_tails_20261029/check_pins.py',
    '05-knowledge/results/source_cutoffs_balanced_lifts_20261029/lifts.json',
)
SOURCES = {
    'f31384340a': (
        *RUNTIME,
        '05-knowledge/results/session110_overnight_20261030/LOGFREE_PRIME_TOURNAMENTS.md',
        '05-knowledge/results/session110_overnight_20261030/CHECKPOINT_02.md',
        '05-knowledge/results/session110_overnight_20261030/MODULE_ENUMERATION.md',
        '05-knowledge/results/session110_overnight_20261030/QUINARY_CARRIERS.md',
        '05-knowledge/results/source_cutoffs_balanced_lifts_20261029/LIFTS.md',
        '05-knowledge/results/regular_parent_symmetry_20261029/PROOFS.md',
        '05-knowledge/results/prime_power_plus_one_symmetry_20261030/PROOFS.md',
        '05-knowledge/results/overnight_rank_return_20261029/COMPLEMENT_KERNEL.md',
        '01-canon/theorems/THM-1920-the-spectral-insertion-response.md',
        '01-canon/theorems/THM-4169-prime-parent-one-vertex-augmentation-and-quartic-johnson-transfer.md',
        '01-canon/theorems/THM-4743-converse-extended-gluck-for-tournaments-primitive-classification-paley-tournaments-and-c3-towers.md',
        '01-canon/theorems/THM-4810-h-spectrum-atlas-strong-fibres-exact-low-spectra-complement-parity-invariant-score-monotonicity-paley-budgets.md',
        '01-canon/theorems/THM-4813-collatz-numerator-coordinate-forcing-cocycle-gauss-sum-intertwiner-and-sharp-spectral-exception-budgets.md',
        '01-canon/theorems/THM-4815-dipole-principle-native-move-distance-bounds-endpoint-relation-cut-or-colour-exactness-and-the-p7-parity-correction.md',
    ),
    '5d01b430aa': ('05-knowledge/results/session110_overnight_20261030/README.md',),
    'de77c39fe5': ('05-knowledge/results/endpoint_shadows_path_corrections_20261030/TOURNAMENTS.md',),
}


def git(*args):
    return subprocess.check_output(['git', *args], cwd=ROOT)


def sha(data):
    return hashlib.sha256(data.replace(b'\r\n', b'\n')).hexdigest()


def snapshot():
    sources = []
    for ref, paths in SOURCES.items():
        commit = git('rev-parse', ref).decode().strip()
        for path in paths:
            spec = commit+':'+path
            source = git('show', spec)
            if path in RUNTIME and sha((ROOT/path).read_bytes()) != sha(source):
                raise RuntimeError('Live runtime input differs from immutable pin: '+path)
            sources.append({'commit': commit, 'path': path,
                            'git_blob': git('rev-parse', spec).decode().strip(), 'lf_sha256': sha(source)})
    return {'status': 'Artifact/source integrity, not independent proof review',
            'artifacts_lf_sha256': {p: sha((HERE/p).read_bytes()) for p in ARTIFACTS},
            'runtime_inputs_checked': list(RUNTIME), 'sources': sources}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--record', action='store_true')
    args = parser.parse_args()
    result = snapshot()
    manifest = HERE/'manifest.json'
    if args.record:
        manifest.write_bytes((json.dumps(result, indent=2, sort_keys=True)+'\n').encode())
    elif result != json.loads(manifest.read_text(encoding='utf-8')):
        raise RuntimeError('Manifest mismatch')
    print(json.dumps({'artifacts': len(ARTIFACTS), 'source_pins': len(result['sources']),
                      'live_runtime_inputs': len(RUNTIME), 'match': True}))


if __name__ == '__main__':
    main()
