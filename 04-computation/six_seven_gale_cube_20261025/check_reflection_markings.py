#!/usr/bin/env python3
"""Audit arbitrary binary relabelling versus the geometric hourglass half-turn.

Standard library only; no producer imports. The complete universe consists of
all 5,040 binary labelings, their 30 Fano kernels, and all 35 branch triples.
Explicit counted gates remain active under Python -O.
"""
from collections import Counter
from itertools import combinations, permutations
from pathlib import Path
import argparse
import hashlib
import json

GATES = 0
LINES = ((1, 0, 0), (1, 1, 0), (1, -1, 0), (1, 0, 1),
         (1, 0, -1), (1, 1, 1), (1, -1, -1))
TRIPLES = tuple(combinations(range(7), 3))
HALF_TURN = (0, 2, 1, 4, 3, 6, 5)


def check(condition, message):
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(message)


def xor(values):
    answer = 0
    for value in values:
        answer ^= value
    return answer


def det(a, b, c):
    return (a[0] * (b[1] * c[2] - b[2] * c[1])
            - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0]))


def j(value):
    return value ^ 7 if value != 7 else 7


def blocks(labels):
    return frozenset(t for t in TRIPLES if xor(labels[i] for i in t) == 0)


def permute_word(word, permutation):
    return sum(1 << permutation[i] for i in range(7) if word >> i & 1)


def permute_blocks(fano, permutation):
    return frozenset(tuple(sorted(permutation[i] for i in t)) for t in fano)


def induced_permutation(labels):
    return tuple(labels.index(j(value)) for value in labels)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    actual = frozenset(t for t in TRIPLES if det(*(LINES[i] for i in t)) == 0)
    check(actual == frozenset(((0, 1, 2), (0, 3, 4), (0, 5, 6))), 'actual incidence')
    for i, line in enumerate(LINES):
        check((line[0], -line[1], -line[2]) == LINES[HALF_TURN[i]], 'linear half-turn')
        check(HALF_TURN[HALF_TURN[i]] == i, 'half-turn involution')
    for t in TRIPLES:
        moved = tuple(sorted(HALF_TURN[i] for i in t))
        check((t in actual) == (moved in actual), 'geometric incidence preservation')

    good = (1, 2, 4, 3, 6, 5, 7)
    reflected_good = tuple(j(value) for value in good)
    arbitrary_permutation = induced_permutation(good)
    check(reflected_good == (6, 5, 3, 4, 1, 2, 7), 'arbitrary reflected labels')
    check(arbitrary_permutation == (4, 5, 3, 2, 0, 1, 6), 'arbitrary branch permutation')
    check(not blocks(good) & actual, 'original marking is free')
    check(blocks(reflected_good) & actual == frozenset(((0, 1, 2),)), 'one bad triple after J')
    check(permute_blocks(actual, arbitrary_permutation) != actual, 'arbitrary permutation hostile')

    adapted = (7, 1, 6, 2, 5, 4, 3)
    check(induced_permutation(adapted) == HALF_TURN, 'adapted marking gives geometric action')
    check(tuple(j(value) for value in adapted) == tuple(adapted[HALF_TURN[i]] for i in range(7)),
          'adapted marked equivariance')
    check(len(blocks(adapted) & actual) == 3, 'adapted marking stratum')
    check(len(blocks(tuple(j(value) for value in adapted)) & actual) == 3,
          'adapted reflection preserves stratum')

    standard = blocks(tuple(range(1, 8)))
    reflected_standard = blocks(tuple(j(value) for value in range(1, 8)))
    check(len(standard & reflected_standard) == 3, 'J preserves exactly three Fano lines')

    kernels, fibers = {}, Counter()
    even_words = [w for w in range(128) if w.bit_count() % 2 == 0]
    for labels in permutations(range(1, 8)):
        fano = blocks(labels)
        kernel = frozenset(w for w in even_words
                           if xor(labels[i] for i in range(7) if w >> i & 1) == 0)
        simplex = frozenset(sum(((a & labels[i]).bit_count() % 2) << i for i in range(7))
                            for a in range(8))
        check(len(kernel) == 8 and kernel == simplex, 'even Hamming kernel/simplex identity')
        check(len(fano) == 7, 'Fano block count')
        check(kernel not in kernels or kernels[kernel] == fano, 'marking fiber consistency')
        kernels[kernel] = fano
        fibers[kernel] += 1
    check(len(kernels) == 30, 'all thirty Fano kernels')
    check(set(fibers.values()) == {168}, 'GL3 marking fibers')

    rows, pair_histogram, stratum_histogram = [], Counter(), Counter()
    keys = sorted(kernels, key=lambda h: tuple(sorted(h)))
    index = {h: i for i, h in enumerate(keys)}
    for kernel in keys:
        moved = frozenset(permute_word(w, HALF_TURN) for w in kernel)
        check(moved in kernels, 'half-turn maps kernel universe to itself')
        check(moved != kernel, 'no Fano kernel is half-turn invariant')
        check(kernels[moved] == permute_blocks(kernels[kernel], HALF_TURN), 'code/incidence equivariance')
        bad = len(kernels[kernel] & actual)
        check(len(kernels[moved] & actual) == bad, 'half-turn preserves bad-triple stratum')
        check(frozenset(permute_word(w, HALF_TURN) for w in moved) == kernel, 'kernel pairing involution')
        stratum_histogram[bad] += 1
        if index[kernel] < index[moved]:
            pair_histogram[bad] += 1
        rows.append({'kernel_index': index[kernel], 'kernel_even_words': sorted(kernel),
                     'bad_triples': bad, 'half_turn_partner': index[moved]})
    check(dict(stratum_histogram) == {0: 16, 1: 12, 3: 2}, 'stratum census')
    check(dict(pair_histogram) == {0: 8, 1: 6, 3: 1}, 'pair census')

    result = {'status': 'FINITE-EXACT reflection integration audit', 'gates': GATES,
              'line_order': LINES, 'actual_triples': sorted(actual),
              'arbitrary_marking': {'labels': good, 'J_labels': reflected_good,
                                    'induced_branch_permutation': arbitrary_permutation,
                                    'bad_before': 0, 'bad_after': 1},
              'adapted_marking': {'labels': adapted, 'J_labels': tuple(j(v) for v in adapted),
                                  'induced_branch_permutation': HALF_TURN,
                                  'bad_before': 3, 'bad_after': 3},
              'labeling_universe': 5040, 'Fano_kernel_universe': 30,
              'strata': dict(sorted(stratum_histogram.items())),
              'half_turn_pairs_by_stratum': dict(sorted(pair_histogram.items())),
              'kernels': rows}
    payload = json.dumps(result, indent=2, sort_keys=True) + '\n'
    args.output.write_text(payload)
    print(json.dumps({'gates': GATES, 'Fano_kernels': 30, 'half_turn_pairs': 15,
                      'sha256': hashlib.sha256(payload.encode()).hexdigest()}, sort_keys=True))


if __name__ == '__main__':
    main()
