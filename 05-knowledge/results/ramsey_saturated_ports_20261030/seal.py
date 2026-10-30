"""Seal inherited Git sources and exact owned artifact bytes."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
REV = '98f961a7aa142a739403de1ba343720371d6125c'
SOURCES = [
    '05-knowledge/results/collatz_catalan_ramsey_20260927.md',
    '05-knowledge/results/affine_tree_recurrences_20261028/AUDIT.md',
    '05-knowledge/results/unifying_binary_dependence_free_bits_7_15_56_20261028/README.md',
    '05-knowledge/results/number_bridges_57_27_16_20261028/ATLAS.md',
    '05-knowledge/results/session35_dimensions_dualities_20261023/SINGER_HEISENBERG_E7_OBSTRUCTION.md',
    '01-canon/theorems/THM-4624-finite-beckman-quarles-for-unit-circle-cayley-graphs-affine-rigidity-except-clebsch-hamming-and-ternary-frame.md',
    '01-canon/theorems/THM-4701-kuratowski-graphs-the-petersen-graph-and-the-finite-unit-circle-graphs.md',
    '05-knowledge/results/regular_parent_symmetry_20261029/PROOFS.md',
    '05-knowledge/results/regular_parent_symmetry_20261029/verify.py',
    '05-knowledge/results/odd_eighth_owner_families_20261030/PROOFS.md',
    '07-reflections/criterion-C-formalized-part-A-reduced-to-realization-plus-paley-QR-resonance-klein-S204.md',
]
INCOMING = {
    '05-knowledge/results/cayley_dickson_integral_20261030/README.md': '7f70dcd601',
    '05-knowledge/results/return_kernel_core_transplants_20261030/README.md': '043b760966',
}


def sha(data):
    return hashlib.sha256(data).hexdigest()


def blob(ref, path):
    return subprocess.check_output(['git', 'show', f'{ref}:{path}'], cwd=ROOT)


def longpath(path):
    value = str(path.resolve())
    return Path('\\\\?\\'+value) if os.name == 'nt' and not value.startswith('\\\\?\\') else Path(value)


def build():
    artifacts = {}
    for path in sorted(HERE.rglob('*')):
        if not path.is_file() or path.name == 'manifest.json' or '__pycache__' in path.parts:
            continue
        data = path.read_bytes()
        if b'\r' in data:
            raise RuntimeError(f'Expected LF bytes: {path.name}')
        artifacts[path.relative_to(HERE).as_posix()] = sha(data)
    return {'source_revision': REV, 'sources': {p: sha(blob(REV, p)) for p in SOURCES},
            'incoming_sources': {p: {'revision': rev, 'sha256': sha(blob(rev, p))} for p, rev in INCOMING.items()},
            'artifacts': artifacts,
            'hash_basis': 'Exact Git source blobs and exact LF artifact bytes.',
            'runtime': 'Native verifiers: Python standard library. Independent audit: NetworkX. Discovery only: python-sat.',
            'proof_scope': 'Same-agent proof audit; independent finite algorithms; no formal proof build.',
            'cited_dependency': 'Completeness of McKay graph banks and the named small Ramsey bounds.'}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--write', action='store_true')
    ap.add_argument('--staged', action='store_true')
    args = ap.parse_args()
    result = build()
    target = HERE/'manifest.json'
    if args.write:
        target.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8', newline='\n')
    elif json.loads(target.read_text()) != result:
        raise RuntimeError('Manifest mismatch')
    for name, digest in result['sources'].items():
        if sha(longpath(ROOT/name).read_bytes().replace(b'\r\n', b'\n')) != digest:
            raise RuntimeError(f'Current dependency differs from pin: {name}')
    # Incoming summaries are historical reading pins, not runtime imports.
    # build() verifies their exact Git blobs; later mainline amendments to
    # these README files must not invalidate an otherwise unchanged packet.
    if args.staged:
        for name, digest in result['artifacts'].items():
            relative = (HERE/name).relative_to(ROOT).as_posix()
            if sha(blob('', relative)) != digest:
                raise RuntimeError(f'Staged artifact differs: {name}')
        if blob('', target.relative_to(ROOT).as_posix()) != target.read_bytes():
            raise RuntimeError('Staged manifest differs')
    print(json.dumps({'sources': len(SOURCES), 'artifacts': len(result['artifacts']),
                      'write': args.write, 'staged': args.staged}))


if __name__ == '__main__':
    main()
