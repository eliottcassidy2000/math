#!/usr/bin/env python3
"""Pinned semantic replay: no saved result receipt is used as an input."""
from hashlib import sha256
from pathlib import Path
import argparse
import json
import sys
import check_pairs as pairs
import check_lift

EXPECTED_RESULTS = '06a34ca92aec14a7f17f62b6ab6c1a81aa0a1fa1648f18d5d777fb9acd25dfbe'
EXPECTED_HOSTILES = '4ed179bdaf8abef534028e1b61d8c20880d020d3ffc133a478dc47c27b95c0d4'
EXPECTED_LIFT = '174fd29e5f4093240461a6e5c249813d215bcf1f8450931ec9365a033eafa854'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--baseline-root')
    args = parser.parse_args()
    old, block = pairs.load_baseline(args.baseline_root)
    results = [pairs.run(n, old, block) for n in range(6)]
    actual = pairs.semantic_hash(results)
    pairs.need(actual == EXPECTED_RESULTS, 'paired-bank semantic replay mismatch')
    hostiles = pairs.hostile_controls(old, block)
    pairs.need(pairs.semantic_hash(hostiles) == EXPECTED_HOSTILES, 'hostile replay mismatch')
    lift = check_lift.audit(old, block)
    pairs.need(lift['semantic_sha256'] == EXPECTED_LIFT, 'self-converse lift replay mismatch')
    return {'status': 'PASS', 'python_optimization_level': sys.flags.optimize,
            'semantic_sha256': actual, 'matches_expected': True,
            'pair_counts': [r['pair_classes'] for r in results],
            'independent_literal_universe_through_n': 4,
            'three_methods_equal_full_output_through_n': 5,
            'hostile_controls': hostiles,
            'lift_semantic_sha256': lift['semantic_sha256'],
            'full_self_converse_image_audit_through_order': 7,
            'source_sha256': {p.name: sha256(p.read_bytes()).hexdigest()
                              for p in sorted(Path(__file__).resolve().parent.glob('*.py'))}}

if __name__ == '__main__':
    print(json.dumps(main(), indent=2, sort_keys=True))
