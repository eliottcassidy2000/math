#!/usr/bin/env python3
"""Packet integrity, exact mode comparisons, source pins and arithmetic guards."""
from pathlib import Path
import hashlib
import json
import re

ROOT = Path(__file__).resolve().parent
IGNORE = {'MANIFEST.json', 'build_nine', 'build_nine_ndebug', 'trace_ten', 'trace_ten_full'}

def need(value, label):
    if not value:
        raise RuntimeError(label)


def main():
    manifest = json.loads((ROOT / 'MANIFEST.json').read_text())
    records = {r['path']: r for r in manifest['files']}
    actual = {str(p.relative_to(ROOT)) for p in ROOT.rglob('*') if p.is_file()
              and p.name not in IGNORE and '__pycache__' not in p.parts
              and not p.name.startswith('verification-')}
    need(set(records) == actual, 'exact payload file membership')
    for path, rec in records.items():
        data = (ROOT / path).read_bytes()
        need(len(data) == rec['bytes'] and hashlib.sha256(data).hexdigest() == rec['sha256'], 'payload hash: ' + path)
    pins = json.loads((ROOT / 'SOURCE_PINS.json').read_text())
    for rec in pins['sources']:
        data = (ROOT / rec['archived_path']).read_bytes()
        need(len(data) == rec['bytes'] and hashlib.sha256(data).hexdigest() == rec['sha256'], 'source pin: ' + rec['archived_path'])
    parent_json = (ROOT / 'sources/certified-parents-eight.json').read_bytes()
    parent_receipt = json.loads((ROOT / 'sources/certified-parents-eight-receipt.json').read_text())
    need(hashlib.sha256(parent_json).hexdigest() == parent_receipt['atlas_sha256'], 'input matches inherited completeness certificate')
    certified = json.loads(parent_json)['parents']
    parent_lines = iter((ROOT / 'parents-eight.txt').read_text().splitlines())
    need(int(next(parent_lines)) == len(certified) == parent_receipt['ordinary_classes'], 'text input matches certified parent cardinality')
    for parent in certified:
        values = list(map(int, next(parent_lines).split()))
        need(values[:8] == parent['rows'] and values[8] == len(parent['actions']), 'text parent agrees with certified JSON')
        actions = [list(map(int, next(parent_lines).split())) for _ in range(values[8])]
        need(actions == parent['actions'], 'text full actions agree with certified JSON')
    need(next(parent_lines, None) is None, 'certified text input exhausted exactly')
    build = json.loads((ROOT / 'build-normal.json').read_text())
    need(build == json.loads((ROOT / 'build-ndebug.json').read_text()), 'normal/ndebug builder equality')
    normal = json.loads((ROOT / 'trace-normal.json').read_text())
    full = json.loads((ROOT / 'trace-full.json').read_text())
    common = set(normal) - {'refinement_rounds', 'singleton_early_exits'}
    need(all(normal[k] == full[k] for k in common), 'early/full trace arithmetic equality')
    need(full['singleton_early_exits'] == 0, 'full mode forbids singleton shortcut')
    need((ROOT / 'trace-normal-per-parent.tsv').read_bytes() == (ROOT / 'trace-full-per-parent.tsv').read_bytes(), 'all per-parent traces equal')
    for name in ('build_nine.cpp', 'trace_ten.cpp'):
        source = (ROOT / name).read_text()
        need(not re.search(r'\b(?:191536|9733056|97336068)\b', source), 'no output-count/trace literals: ' + name)
        need('assert(' not in source and 'assert (' not in source, 'guards are not disabled assertions: ' + name)
    trace_source = (ROOT / 'trace_ten.cpp').read_text()
    need('canonical_eight' not in trace_source and 'pointed_key' not in trace_source, 'trace has no parent reroot engine')
    need(build['order_nine_canonicalizations'] == 0 and build['maximum_canonicalized_order'] == 8, 'builder canonicalization boundary')
    need(all(normal[k] == 0 for k in ('target_canonicalizations', 'target_deletions', 'target_automorphism_computations')), 'trace target boundary')
    correction = json.loads((ROOT / 'sources/structural-correction-receipt.json').read_text())['total_error_order_ten']
    need(normal['structural_error_correction'] == correction, 'correction source agreement')
    histogram = dict(normal['stable_root_cell_histogram'])
    need(sum(histogram.values()) == normal['pointed_extension_states'], 'pointed histogram sum')
    need(sum(k * v for k, v in histogram.items()) == normal['upper_refinement_trace'], 'weighted root-cell histogram')
    need((normal['upper_refinement_trace'] - correction) == 10 * normal['predicted_tournament_count'], 'exact prediction identity')
    need(build['labelled_orbit_mass'] == normal['parent_labelled_mass'] == 1 << 36, 'parent labelled universe mass')
    lines = iter((ROOT / 'parents-nine.txt').read_text().splitlines())
    count = int(next(lines)); need(count == build['constructed_classes'] == normal['parent_classes'], 'constructed parent count')
    mass = 0
    for _ in range(count):
        row = list(map(int, next(lines).split())); need(len(row) == 10, 'parent output row length')
        actions = row[-1]; need(actions > 0 and 362880 % actions == 0, 'parent action order')
        mass += 362880 // actions
        for _ in range(actions):
            action = list(map(int, next(lines).split())); need(sorted(action) == list(range(9)), 'output action permutation')
    need(next(lines, None) is None and mass == 1 << 36, 'complete generated input and mass')
    print(json.dumps({'status': 'PASS', 'payload_files': len(records), 'sources': len(pins['sources']),
                      'constructed_parent_classes': count, 'upper_trace': normal['upper_refinement_trace'],
                      'prediction': normal['predicted_tournament_count'], 'correction': correction,
                      'normal_ndebug_build_equal': True, 'all_per_parent_traces_equal': True}, sort_keys=True))

if __name__ == '__main__':
    main()
