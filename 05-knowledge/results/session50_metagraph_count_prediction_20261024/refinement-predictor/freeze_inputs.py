#!/usr/bin/env python3
"""One-time INPUT preparation, explicitly NOT the blind prediction entrypoint.

Loads the literal pinned older growth source and stops at order seven.
The original engine may perform parent/child isomorphism and root-transport
operations here. None of this code is imported by predictor.py.
"""
import argparse
import hashlib
import importlib.machinery
import importlib.util
import json
from pathlib import Path


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', default='parents-through-seven.json')
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    source = here / 'archive/growth-transfer.py.txt'
    source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
    pins = json.loads((here / 'source-pins.json').read_text())
    require(source_hash == pins['sources'][0]['sha256'], 'frozen growth source hash')
    loader = importlib.machinery.SourceFileLoader('frozen_growth_input_only', str(source))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    reps = [()]
    levels = [{'order': 0, 'parents': [{'rows': [], 'actions': [[]]}]}]
    for order in range(1, 8):
        reps, report, unused = module.predict(reps)
        require(report['order'] == order, 'input generation order cap')
        levels.append({'order': order, 'parents': [
            {'rows': list(rows), 'actions': [list(g) for g in module.automorphisms(rows)]}
            for rows in sorted(reps)]})
    require([len(level['parents']) for level in levels] == [1, 1, 1, 2, 4, 12, 56, 456],
            'validated inherited prefix through seven')
    payload = {'schema': 'tournament-parent-atlas-actions-v1',
               'maximum_parent_order': 7, 'source_sha256': source_hash,
               'levels': levels}
    output = Path(args.output)
    output.write_text(json.dumps(payload, sort_keys=True, separators=(',', ':')) + '\n')
    print(json.dumps({'status': 'PASS', 'maximum_generated_order': 7,
                      'level_sizes': [len(level['parents']) for level in levels],
                      'atlas_sha256': hashlib.sha256(output.read_bytes()).hexdigest(),
                      'frozen_growth_source_sha256': source_hash}, sort_keys=True))


if __name__ == '__main__':
    main()
