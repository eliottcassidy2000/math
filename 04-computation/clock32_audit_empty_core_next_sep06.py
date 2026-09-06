#!/usr/bin/env python3
"""Independent exhaustive finite-mask audit; no primary producer imported."""
from hashlib import sha256
from itertools import combinations_with_replacement, product
import json
from math import gcd
from pathlib import Path
import subprocess

GATES = 0


def check(predicate, message):
    global GATES
    GATES += 1
    if not predicate:
        raise RuntimeError(message)


def blocks(q):
    k = (q + 6) // 7
    answer = set()
    for step in range(q):
        if gcd(step, q) != 1:
            continue
        for start in range(q):
            row = frozenset((start + step * j) % q for j in range(k))
            check(len(row) == k, 'unit orbit repeats early')
            answer.add(row)
    return sorted(answer, key=lambda s: tuple(sorted(s)))


def lift(row, q, c):
    return frozenset(j for j in range(c) if j % q in row)


def main():
    universe = {q: blocks(q) for q in (8, 10, 15, 16, 32)}
    lifted = {q: [lift(row, q, 32) for row in universe[q]]
              for q in (8, 16, 32)}
    for q in lifted:
        for row in lifted[q]:
            check(len(row) == (32 // q) * ((q + 6) // 7),
                  'wrong padded size')
    for row in universe[32]:
        check(len({x % 8 for x in row}) == 5,
              'five distinct residue classes')

    pair_minima = {}
    for q, r in ((8, 8), (8, 16), (8, 32), (16, 32), (32, 32)):
        minimum = min(len(a & b) for a, b in product(lifted[q], lifted[r]))
        check(minimum == 0, 'unexpected pair overlap')
        pair_minima[f'{q},{r}'] = minimum

    maximum = -1
    disjoint_pairs = 0
    triples = 0
    for a, b in combinations_with_replacement(lifted[8], 2):
        if not a & b:
            disjoint_pairs += 1
            check(len({j % 8 for j in a | b}) == 4, 'four classes')
        for z in lifted[32]:
            triples += 1
            overlap_credit = len(a) + len(b) + len(z) - len(a | b | z)
            check(overlap_credit >= 1, 'three-mask obstruction')
            if not a & b:
                check(bool(z & (a | b)), 'pigeonhole intersection')
            maximum = max(maximum, len(a | b | z))
    check(maximum == 20, 'exact triple union maximum')
    check(maximum + 5 + 6 == 31, 'full signature union bound')

    partition_rows = [(15, frozenset(s)) for s in ((2, 3, 4), (7, 8, 9),
                                                 (12, 13, 14))]
    partition_rows += [(10, frozenset(s)) for s in ((0, 1), (5, 6))]
    partition = []
    for q, row in partition_rows:
        check(row in universe[q], 'partition is an affine capacity block')
        partition.append(lift(row, q, 30))
    for j in range(30):
        check(sum(j in row for row in partition) == 1,
              'literal partition of Z30')

    # Freeze the already independently audited arithmetic reduction as data.
    revision = 'b77d7fb5b7cf3d5e2ea7831ae8f720e632576f0f'
    inherited_path = '05-knowledge/results/gcd_pair_hunter_audit_empty_core_next_sep06.json'
    repository = Path(__file__).resolve().parents[1]
    raw = subprocess.check_output(['git', 'show', revision + ':' + inherited_path],
                                  cwd=repository)
    check(sha256(raw).hexdigest() ==
          'dde323ce07917e530c562c2e9994386e005b5e1e42fa3fb7a9a93f93d0544a49',
          'inherited certificate hash')
    layers = json.loads(raw)['layers']
    five = layers[4]
    check(five['tails'] == 5 and len(five['profiles']) == 110,
          'five-tail universe')
    check([r for r in five['profiles'] if r[0] >= 30] ==
          [[30, [2, 2, 2, 3, 3]], [32, [1, 1, 2, 4, 4]]],
          'complete boundary profiles')
    remaining_clocks = [c for c in five['clocks'] if c != 32]
    check(max(remaining_clocks) == 30, 'new eight-body bound')
    six_removed = [row for row in layers[5]['profiles'] if 32 in row[1]]
    check(six_removed == [[32, [1, 1, 2, 4, 4, 32]],
                          [64, [1, 1, 2, 4, 4, 32]]],
          'immediate inherited child exclusions')
    check(max(c for c, g in layers[5]['profiles'] if 32 not in g) == 96,
          'no claimed six-tail cap improvement from this filter alone')

    manifest = {
        'block_counts': {q: len(rows) for q, rows in universe.items()},
        'pair_minima': pair_minima,
        'disjoint_order8_pairs': disjoint_pairs,
        'triple_configurations': triples,
        'maximum_triple_union': maximum,
        'full_clock32_union_bound': 31,
        'clock30_partition': [sorted(row) for row in partition],
        'new_five_tail_clocks': remaining_clocks,
        'six_tail_immediate_exclusions': six_removed,
    }
    payload = json.dumps(manifest, sort_keys=True, separators=(',', ':'))
    print('STATUS: independent exact audit PASS; actual LRC(14) remains OPEN')
    print('UNIVERSE: every affine unit capacity block at q=8,10,15,16,32;')
    print('          every two q8 masks with one q32 mask; one q30 partition')
    print(json.dumps(manifest, sort_keys=True, indent=2))
    print('gates:', GATES)
    print('semantic_sha256:', sha256(payload.encode()).hexdigest())


if __name__ == '__main__':
    main()
