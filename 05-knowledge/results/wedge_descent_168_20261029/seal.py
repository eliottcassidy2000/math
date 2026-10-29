"""Pin source Git blobs and LF artifact bytes; --check never mutates files."""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
REV = '96e3f3fe347ac3060a09c88988809e8e5aacfaad'
SOURCE_PATHS = [
    '05-knowledge/results/session103_words_necklaces_20261029/README.md',
    '05-knowledge/results/session103_words_necklaces_20261029/triangularity/AVOIDING_WORD_TRIANGULARITY.md',
    '05-knowledge/results/session103_words_necklaces_20261029/bridge/BRIDGE_SINK_STAR_LEMMA.md',
    '05-knowledge/results/operation_closed_synthesis_20261029/THEORY.md',
    '05-knowledge/results/session102_tournament_classification_20261029/switching/SIGNED_CUTS_AND_SIMPLE_SYMMETRY.md',
    '01-canon/theorems/THM-4553-the-point-stripping-ladder-psl27-borel-torus-and-the-collatz-frobenius.md',
    '01-canon/theorems/THM-2753-six-edge-parity-erasure-and-three-matching-resolvent-restoration.md',
    '01-canon/theorems/THM-3045-k4-edge-isotypic-binary-ternary-integral-clutch.md',
    '01-canon/theorems/THM-4507-finite-valuation-collatz-potential-obstruction.md',
    '01-canon/theorems/THM-2422-operation-fibres-summand-closure-and-twin-center-ancestry.md',
    '01-canon/theorems/THM-2433-operation-fibre-deletion-incidence-and-startup-scar.md',
    '05-knowledge/results/arithmetic_braids_20260917_summand.md',
    '05-knowledge/results/collatz_mod6_20260922_w6_square_sum_hamiltonicity.md',
    '05-knowledge/results/ap_square_pythagorean_20261019.md',
    '05-knowledge/results/collatz_inductive_floor_receipts_20261005.md',
    '05-knowledge/results/collatz_atomic_rank_interface_20261011.md',
    '05-knowledge/results/weighted_crt_clocks_20261028/CLOCKS.md',
    '05-knowledge/results/cursor_cm_phase_20261028/TRACE_CARRIER.md',
    '05-knowledge/reference/CORE-PAPERS-ABC-IUT.md',
    '05-knowledge/results/four_vertex_glue_fermat_20261029/verify.py',
    '05-knowledge/results/overnight_rank_return_20261029/two_run.json',
]


def sha(data):
    return hashlib.sha256(data).hexdigest()


def git_blob(ref, path):
    return subprocess.check_output(['git', 'show', f'{ref}:{path}'], cwd=ROOT)


def build():
    # Discover transitive local Python imports actually loaded by the hostile audit.
    spec = importlib.util.spec_from_file_location('wedge_probe_manifest', HERE/'verify.py')
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.inherited_hostile()
    sources = set(SOURCE_PATHS)
    for module in tuple(sys.modules.values()):
        filename = getattr(module, '__file__', None)
        if not filename:
            continue
        file = Path(filename).resolve()
        if file.suffix == '.py' and file.is_relative_to(ROOT) and not file.is_relative_to(HERE):
            sources.add(file.relative_to(ROOT).as_posix())
    artifacts = {}
    for file in sorted(HERE.iterdir()):
        if not file.is_file() or file.name == 'manifest.json':
            continue
        data = file.read_bytes()
        if b'\r' in data:
            raise RuntimeError(f'Artifact is not LF: {file.name}')
        artifacts[file.name] = sha(data)
    return {'source_revision': REV,
            'hash_basis': 'inherited sources: exact Git blob; packet artifacts: exact LF bytes',
            'sources': {p: sha(git_blob(REV, p)) for p in sorted(sources)},
            'artifacts': artifacts,
            'audit_scope': 'same-agent proof and independent algorithms; no independent referee'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--staged', action='store_true')
    args = parser.parse_args()
    result = build()
    if args.check:
        saved = json.loads((HERE/'manifest.json').read_text(encoding='utf-8'))
        if saved != result:
            raise RuntimeError('Source/artifact manifest differs')
        # Also check the local inherited checkout agrees with pinned source blobs
        # after Git's text normalization. Never silently inherit a changed theorem.
        for path, expected in result['sources'].items():
            data = (ROOT/path).read_bytes().replace(b'\r\n', b'\n')
            if sha(data) != expected:
                raise RuntimeError(f'Inherited source changed: {path}')
    else:
        (HERE/'manifest.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8', newline='\n')
    if args.staged:
        for name, expected in result['artifacts'].items():
            path = (HERE/name).relative_to(ROOT).as_posix()
            if sha(git_blob('', path)) != expected:
                raise RuntimeError(f'Staged blob hash mismatch: {path}')
    print(json.dumps({'sources': len(result['sources']), 'artifacts': len(result['artifacts']),
                      'check': args.check, 'staged': args.staged}))


if __name__ == '__main__':
    main()
