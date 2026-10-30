"""Hash current packet bytes and pinned inherited Git blobs; verify without mutation."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BASE = 'cd6d98ba58a2418fb81054659bc37962676c8dd0'
INCOMING = 'db3daffd7c'
LATE = '4b423e24b8'
BASE_PATHS = [
    '01-canon/theorems/LEM-004-tournaments-are-odd-functions.md',
    '01-canon/theorems/THM-466-the-2-adic-digits-of-H-are-odd-cycle-collection-counts.md',
    '01-canon/theorems/THM-1380-free-Z2-but-no-odd-map.md',
    '05-knowledge/results/odd_prime_function_cycle_20261018d.md',
    '05-knowledge/results/attachment_lattices_square_sum_20261023/square_sum/SHARED_CLOCKS.md',
    '05-knowledge/results/aliquot_66_five_fifteen_20261028/PROOFS.md',
    '05-knowledge/results/collatz_mod6_20260917_floor_sums_odd_functions.md',
    '05-knowledge/results/odd_sign_spectra_trijections_20261028/LRC_UNIT_DISTANCE.md',
    '05-knowledge/results/lucas_mirror_gap_spaces_20261028/SYNTHESIS.md',
    '05-knowledge/results/aliquot_power_124_20261026/tournaments/TOURNAMENTS_124.md',
    '05-knowledge/results/liouville_mobius_11_13_20261026/RECOVERY.md',
    '05-knowledge/results/session27_repunits_lattices_20261022/session27-response/ODD_CYCLE_FAITHFULNESS.md',
    '05-knowledge/results/regular_parent_symmetry_20261029/verify.py',
    '05-knowledge/results/terminal_fees_phase_covers_20261030/README.md',
    '05-knowledge/results/terminal_fees_phase_covers_20261030/ADAPTIVE.md',
    '05-knowledge/results/rectangular_counting_compiler_20261030/README.md',
    '05-knowledge/results/scaled_count_fft_carriers_20261030/THEORY.md',
]
INCOMING_PATHS = [
    '05-knowledge/results/session118_shared_gate_deletions_20261030/README.md',
    '05-knowledge/results/session118_shared_gate_deletions_20261030/frontier/FRONTIER.md',
]
SOURCES = {**dict.fromkeys(BASE_PATHS, BASE), **dict.fromkeys(INCOMING_PATHS, INCOMING)}
SOURCES.update(dict.fromkeys([
    '05-knowledge/results/session119_odd_transport_classifications_20261030/README.md',
    '05-knowledge/results/session119_odd_transport_classifications_20261030/cycles/SIGNED_GRAPH_INHOMOGENEOUS.md',
    '05-knowledge/results/session119_odd_transport_classifications_20261030/arithmetic/WEIGHTED_CYCLE_OBSTRUCTIONS.md',
], LATE))


def sha(data):
    return hashlib.sha256(data).hexdigest()


def git(*args):
    return subprocess.check_output(['git', *args], cwd=ROOT)


def long_path(p):
    path = str(p.resolve())
    return Path('\\\\?\\'+path) if os.name == 'nt' and not path.startswith('\\\\?\\') else Path(path)


def build():
    artifacts = {}
    for p in sorted(HERE.rglob('*')):
        if not p.is_file() or p.name == 'manifest.json' or '__pycache__' in p.parts:
            continue
        data = p.read_bytes()
        if b'\r' in data:
            raise RuntimeError('Packet artifact is not LF: '+p.name)
        artifacts[p.relative_to(HERE).as_posix()] = sha(data)
    sources = {path: {'revision': git('rev-parse', rev).decode().strip(),
                      'sha256': sha(git('show', f'{rev}:{path}'))}
               for path, rev in SOURCES.items()}
    return {'sources': sources, 'artifacts': artifacts,
            'basis': 'Exact Git source blobs and exact LF packet bytes.',
            'status': 'PROVED scoped synthesis / FINITE-EXACT / OPEN major applications.',
            'review': 'Same-agent review; separate exact algorithms; no formal build.',
            'runtime': 'Python standard library and pinned inherited native 81/27 helper.'}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--write', action='store_true')
    ap.add_argument('--staged', action='store_true')
    ap.add_argument('--allow-incoming-missing', action='store_true')
    args = ap.parse_args()
    result = build()
    target = HERE/'manifest.json'
    if args.write:
        target.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8', newline='\n')
    elif json.loads(target.read_text(encoding='utf-8')) != result:
        raise RuntimeError('Manifest mismatch')
    for path, entry in result['sources'].items():
        local = long_path(ROOT/path)
        if not local.exists() and args.allow_incoming_missing and path in INCOMING_PATHS:
            continue
        if sha(local.read_bytes().replace(b'\r\n', b'\n')) != entry['sha256']:
            raise RuntimeError('Current dependency differs from pin: '+path)
    if args.staged:
        for name, digest in result['artifacts'].items():
            rel = (HERE/name).relative_to(ROOT).as_posix()
            if sha(git('show', ':'+rel)) != digest:
                raise RuntimeError('Staged artifact differs: '+name)
        if git('show', ':'+target.relative_to(ROOT).as_posix()) != target.read_bytes():
            raise RuntimeError('Staged manifest differs')
    if (HERE/'certificate.json').read_bytes() != (HERE/'certificate.optimized.json').read_bytes():
        raise RuntimeError('Normal and optimized receipts differ')
    print(json.dumps({'sources': len(SOURCES), 'artifacts': len(result['artifacts']),
                      'write': args.write, 'staged': args.staged,
                      'normal_optimized_equal': True}))


if __name__ == '__main__':
    main()
