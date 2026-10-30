"""Check normalized artifact hashes and immutable source Git blobs."""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
MANIFEST = HERE/'manifest.json'
ARTIFACTS = ('README.md', 'PROOFS.md', 'AUDIT.md', 'check.py', 'certificate.json', 'seal.py')
SOURCES = {
    '0ad8de0355': (
        'regular_parent_symmetry_20261029/PROOFS.md',
        'regular_parent_symmetry_20261029/HANDOFF.md',
        'regular_parent_symmetry_20261029/verify.py',
        'stabilizer_necklace_kernels_20261028/STABILIZERS.md',
        'fermat_wieferich_order_20261014e.md',
        'charged_recurrence_memory_20261028/FROBENIUS.md',
        'growing_attachment_families_20261029/PROOFS.md',
        'session110_overnight_20261030/CHECKPOINT_02.md',
    ),
    '1f52ab41e6': (
        'prime_power_plus_one_20261029/GROWTH.md',
        'prime_power_plus_one_20261029/ARITHMETIC.md',
        'prime_power_plus_one_20261029/HANDOFF.md',
        'source_cutoffs_balanced_lifts_20261029/README.md',
    ),
}


def git(*args):
    return subprocess.check_output(['git', *args], cwd=ROOT)


def digest(data):
    return hashlib.sha256(data.replace(b'\r\n', b'\n')).hexdigest()


def snapshot():
    sources = []
    for ref, paths in SOURCES.items():
        commit = git('rev-parse', ref).decode().strip()
        for suffix in paths:
            path = '05-knowledge/results/'+suffix
            spec = commit+':'+path
            sources.append({'commit': commit, 'path': path,
                            'git_blob': git('rev-parse', spec).decode().strip(),
                            'lf_sha256': digest(git('show', spec))})
    return {'status': 'Immutable provenance and artifact integrity, not independent proof review',
            'artifacts_lf_sha256': {p: digest((HERE/p).read_bytes()) for p in ARTIFACTS},
            'sources': sources}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--record', action='store_true')
    args = parser.parse_args()
    current = snapshot()
    if args.record:
        MANIFEST.write_text(json.dumps(current, indent=2, sort_keys=True)+'\n', encoding='utf-8')
    elif current != json.loads(MANIFEST.read_text(encoding='utf-8')):
        raise RuntimeError('Artifact or immutable source differs from the recorded manifest')
    print(json.dumps({'artifacts': len(current['artifacts_lf_sha256']),
                      'source_pins': len(current['sources']), 'match': True}))


if __name__ == '__main__':
    main()
