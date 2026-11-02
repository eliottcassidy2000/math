"""Pin inherited Git objects and seal this packet's LF-normalized bytes."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SOURCE_REV = 'c16226249f1f27a335cf7c84afdeaa90879b9580'
SOURCES = [
    '05-knowledge/results/stability_priority_moduli_20261102/STABILITY_N500.md',
    '05-knowledge/results/stability_priority_moduli_20261102/STABILITY_N276.md',
    '05-knowledge/results/stability_priority_moduli_20261102/FRACTRAN_PRECISION.md',
    '05-knowledge/results/session138_stability_and_actual_contexts_20261102/threshold_refinement/PROOF_N1500.md',
    '05-knowledge/results/session138_stability_and_actual_contexts_20261102/asymptotic_stability/PROOF.md',
    '05-knowledge/results/session138_stability_and_actual_contexts_20261102/stability/MULTISTAR_STABILITY.md',
    '05-knowledge/results/session138_stability_and_actual_contexts_20261102/stability/TRIANGLE_VARIANCE.md',
    '05-knowledge/results/session138_stability_and_actual_contexts_20261102/stability/H1_STABILITY.md',
    '05-knowledge/results/session138_stability_and_actual_contexts_20261102/stability/SHARP_LOW_ORDER_STABILITY.md',
    '05-knowledge/results/session138_stability_and_actual_contexts_20261102/rounding_graph/PROOF.md',
    '05-knowledge/results/session139_incidence_moduli_and_stability_20261102/stability/PROOF_N526.md',
    '05-knowledge/results/session139_incidence_moduli_and_stability_20261102/independent_threshold526/AUDIT.md',
    '05-knowledge/results/session139_incidence_moduli_and_stability_20261102/all_moduli/PROOF.md',
    '05-knowledge/results/session139_incidence_moduli_and_stability_20261102/riemann_five_eleven/WEIGHTED_ZETA_BOUNDARY.md',
    '05-knowledge/results/session136_full_tournament_spectrum_20261102/local_closure/PROOF.md',
    '05-knowledge/results/session135_fractran_tm_structure_20261101/SYNTHESIS.md',
    '05-knowledge/results/primegame_variant_selector_automaton_metagraph_20261031.md',
    '05-knowledge/results/cursor_cm_phase_20261028/TRACE_CARRIER.md',
    '05-knowledge/results/cursor_cm_phase_20261028/RAMANUJAN.md',
    '05-knowledge/results/regular_parent_symmetry_20261029/verify.py',
    '05-knowledge/results/ninja_pairs_square_sum_number_system_20261031.md',
    '01-canon/theorems/THM-002-ocf.md',
    '01-canon/theorems/THM-209-independence-polynomial.md',
    '01-canon/theorems/THM-326-staircase-hp-equals-ip.md',
    '01-canon/theorems/THM-466-the-2-adic-digits-of-H-are-odd-cycle-collection-counts.md',
    '01-canon/theorems/THM-505-the-ocf-non-spectral-defect-H-equals-spectral-skeleton-plus-witt-defect.md',
    '01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md',
    '01-canon/theorems/THM-1780-H-leaves-the-spectral-ladder-at-n6.md',
    '01-canon/theorems/THM-4729-closed-and-open-tilings-hamiltonian-cycle-bipartition-of-the-metagraph-lines-and-the-forced-arrangement.md',
]


def normalized(data):
    return data.replace(b'\r\n', b'\n')


def digest(data):
    return hashlib.sha256(normalized(data)).hexdigest()


def local_bytes(path):
    raw = str(path.resolve())
    if os.name == 'nt' and not raw.startswith('\\\\?\\'):
        raw = '\\\\?\\'+raw
    return Path(raw).read_bytes()


def git(*args):
    return subprocess.check_output(['git', *args], cwd=ROOT)


def source_records():
    records = []
    for path in SOURCES:
        data = git('show', SOURCE_REV+':'+path)
        records.append({'path': path, 'revision': SOURCE_REV,
                        'git_blob': git('rev-parse', SOURCE_REV+':'+path).decode().strip(),
                        'sha256_lf': digest(data),
                        'current_file_must_match': not path.startswith('01-canon/theorems/THM-1780-')})
    return {'normalization': 'CRLF to LF', 'sources': records,
            'historical_exception': 'THM-1780 is explicitly repaired in this session; its original is pinned for lineage.'}


def artifact_records():
    paths = sorted(p for p in HERE.iterdir() if p.is_file() and p.name != 'SHA256SUMS.json')
    return {p.name: digest(local_bytes(p)) for p in paths}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    source_path, seal_path = HERE/'SOURCE_PINS.json', HERE/'SHA256SUMS.json'
    if not args.check:
        source_path.write_text(json.dumps(source_records(), indent=2)+'\n', encoding='utf-8')
        seal_path.write_text(json.dumps(artifact_records(), indent=2)+'\n', encoding='utf-8')
    actual, expected = source_records(), json.loads(source_path.read_text(encoding='utf-8'))
    if actual != expected:
        raise RuntimeError('pinned Git objects do not match source manifest')
    for item in expected['sources']:
        if item['current_file_must_match'] and digest(local_bytes(ROOT/item['path'])) != item['sha256_lf']:
            raise RuntimeError('inherited working input changed: '+item['path'])
    if artifact_records() != json.loads(seal_path.read_text(encoding='utf-8')):
        raise RuntimeError('artifact seal mismatch')
    print(json.dumps({'status': 'PASS', 'pinned_sources': len(actual['sources']),
                      'sealed_artifacts': len(artifact_records()), 'revision': SOURCE_REV}))


if __name__ == '__main__':
    main()
