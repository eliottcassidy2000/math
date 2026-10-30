"""Pin inherited Git blobs and exact LF packet bytes; checking never mutates."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
REV = '5ab37f308003e8634d13ea4c0dfe8f1558c11054'
SOURCES = [
    '01-canon/theorems/THM-4803-intransitivity-is-the-relation-modulo-its-modular-tree-and-a000568-from-the-prime-cycle-index-census.md',
    '05-knowledge/results/stabilizer_necklace_kernels_20261028/STABILIZERS.md',
    '05-knowledge/results/exception_palette_dynamics_20261021/session24-palette-dynamics/PALETTE_DYNAMICS.md',
    '05-knowledge/results/session68_tournament_composition_classification_20261026/session68-classification/CLASSIFICATION.md',
    '05-knowledge/results/session107_quasirandom_arithmetic_20261029/colored/COLORED_RECURRENCE.md',
    '05-knowledge/results/session108_sigma_symmetry_20261030/tournaments/PALETTE_SYMMETRY.md',
    '05-knowledge/results/session109_pisano_symmetry_20261030/README.md',
    '05-knowledge/results/growing_attachment_families_20261029/README.md',
    '05-knowledge/results/wedge_descent_168_20261029/THEORY.md',
    '05-knowledge/results/tournament_open_problem_routes_20261029/PROOFS.md',
    '05-knowledge/results/overnight_rank_return_20261029/atom_precision.py',
]


def sha(data):
    return hashlib.sha256(data).hexdigest()


def blob(ref, path):
    return subprocess.check_output(['git', 'show', f'{ref}:{path}'], cwd=ROOT)


def local_file(path):
    resolved = str((ROOT/path).resolve())
    return Path('\\\\?\\'+resolved) if os.name == 'nt' and not resolved.startswith('\\\\?\\') else Path(resolved)


def build():
    artifacts = {}
    for file in sorted(HERE.iterdir()):
        if not file.is_file() or file.name == 'manifest.json':
            continue
        data = file.read_bytes()
        if b'\r' in data:
            raise RuntimeError(f'Non-LF artifact: {file.name}')
        artifacts[file.name] = sha(data)
    return {
        'source_revision': REV,
        'hash_basis': 'Inherited sources: exact Git blobs. Packet artifacts: exact LF bytes.',
        'sources': {path: sha(blob(REV, path)) for path in sorted(SOURCES)},
        'artifacts': artifacts,
        'scope': 'Same-agent proof audit; independent native algorithms; no independent referee.',
        'runtime': 'Python standard library only; verify.py imports no local research modules.',
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--staged', action='store_true')
    args = parser.parse_args()
    result = build()
    target = HERE/'manifest.json'
    if args.write:
        target.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8', newline='\n')
    elif json.loads(target.read_text(encoding='utf-8')) != result:
        raise RuntimeError('Saved manifest differs from source or artifact bytes')
    for path, expected in result['sources'].items():
        if sha(local_file(path).read_bytes().replace(b'\r\n', b'\n')) != expected:
            raise RuntimeError(f'Inherited checkout source changed: {path}')
    if args.staged:
        for name, expected in result['artifacts'].items():
            path = (HERE/name).relative_to(ROOT).as_posix()
            if sha(blob('', path)) != expected:
                raise RuntimeError(f'Staged artifact differs: {path}')
    print(json.dumps({'sources': len(result['sources']), 'artifacts': len(result['artifacts']),
                      'write': args.write, 'staged': args.staged}))


if __name__ == '__main__':
    main()
