#!/usr/bin/env python3
"""Exact six/seven root incidence, free deck groups, and Paley-Fano bridges.

Standard library only. Two independent seven-root paths: every rank-three
subgroup of F2^6 versus all Fano-labelled branch maps. All comparisons are
explicit exceptions and survive -O. Geometry is proved in the companion note.
"""
from collections import Counter
from itertools import combinations, permutations
from pathlib import Path
import argparse
import hashlib
import json

LINES = ((1, 0, 0), (1, 1, 0), (1, -1, 0), (1, 0, 1), (1, 0, -1),
         (1, 1, 1), (1, -1, -1), (1, 1, -1), (1, -1, 1))
HOURGLASS_TRIPLES = frozenset(((0, 1, 2), (0, 3, 4), (0, 5, 6)))


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def determinant(a, b, c):
    return (a[0] * (b[1] * c[2] - b[2] * c[1])
            - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0]))


def span(generators):
    result = {0}
    for g in generators:
        result |= {x ^ g for x in tuple(result)}
    return frozenset(result)


def inertia_generators(n):
    return tuple(1 << i for i in range(n - 1)) + ((1 << (n - 1)) - 1,)


def forbidden_inertia(n, actual):
    generators = inertia_generators(n)
    # Every single and pair of distinct branch lines occurs in P2.
    supports = list(combinations(range(n), 1)) + list(combinations(range(n), 2)) + list(actual)
    return frozenset().union(*(span(generators[i] for i in support) for support in supports)) - {0}


def subgroup_census(dimension, rank):
    return {span(basis) for basis in combinations(range(1, 1 << dimension), rank)
            if len(span(basis)) == 1 << rank}


def fano_structures():
    structures = {}
    label_histogram = Counter()
    for labeling in permutations(range(1, 8)):
        blocks = frozenset(t for t in combinations(range(7), 3)
                           if labeling[t[0]] ^ labeling[t[1]] ^ labeling[t[2]] == 0)
        check(len(blocks) == 7, 'Fano block count')
        pair_counts = Counter(p for block in blocks for p in combinations(block, 2))
        check(len(pair_counts) == 21 and set(pair_counts.values()) == {1}, 'not STS(7)')
        labels = labeling[:6]
        kernel = frozenset(v for v in range(64)
                           if xor_all(labels[i] for i in range(6) if v >> i & 1) == 0)
        generators = inertia_generators(7)
        line_kernel = frozenset({0} | {xor_all(generators[i] for i in block) for block in blocks})
        check(kernel == line_kernel and len(kernel) == 8, 'Fano/Hamming kernel mismatch')
        if blocks in structures:
            check(structures[blocks]['kernel'] == kernel, 'post-label kernel changed')
            structures[blocks]['labelings'] += 1
        else:
            structures[blocks] = {'labeling': labeling, 'kernel': kernel, 'labelings': 1}
        label_histogram[len(blocks & HOURGLASS_TRIPLES)] += 1
    check(len(structures) == 30, 'wrong Fano structure count')
    check({v['labelings'] for v in structures.values()} == {168}, 'GL3 orbit size')
    check(dict(label_histogram) == {0: 2688, 1: 2016, 3: 336}, 'label histogram')
    return structures, label_histogram


def xor_all(values):
    output = 0
    for value in values:
        output ^= value
    return output


def seven_checks(structures):
    all_subgroups = subgroup_census(6, 3)
    check(len(all_subgroups) == 1395, 'rank-three subgroup universe')
    forbidden = forbidden_inertia(7, HOURGLASS_TRIPLES)
    direct = {h for h in all_subgroups if not h & forbidden}
    predicted = {v['kernel'] for f, v in structures.items() if not f & HOURGLASS_TRIPLES}
    check(direct == predicted and len(direct) == 16, 'independent free subgroup disagreement')
    # Canonical characters are even subsets omitting one of seven labels.
    characters = [63 ^ (1 << i) for i in range(6)] + [63]
    for info in structures.values():
        group = info['kernel']
        check(not any(all((c & h).bit_count() % 2 == 0 for h in group) for c in characters),
              'unexpected invariant canonical character')
        for h in group - {0}:
            check(sum((-1) ** ((c & h).bit_count() % 2) for c in characters) == -1,
                  'canonical augmentation trace')
    return {'all_rank_three_subgroups': len(all_subgroups), 'free_subgroups': len(direct),
            'canonical_invariants': 0, 'nonidentity_canonical_trace': -1,
            'good_labeling': [1, 2, 4, 3, 6, 5, 7],
            'bad_pencil_labeling': [1, 2, 3, 4, 5, 6, 7]}


def subset_checks(structures):
    actual_global = frozenset(t for t in combinations(range(9), 3)
                              if determinant(*(LINES[i] for i in t)) == 0)
    check(len(actual_global) == 8, 'nine-line triple count')
    for four in combinations(range(9), 4):
        check(any(determinant(*(LINES[i] for i in triple)) != 0 for triple in combinations(four, 3)),
              'fourfold concurrence')
    six, seven = [], []
    six_histogram, seven_histogram = Counter(), Counter()
    for subset in combinations(range(9), 6):
        positions = {v: i for i, v in enumerate(subset)}
        actual = frozenset(tuple(positions[v] for v in t) for t in actual_global if set(t) <= set(subset))
        forbidden = forbidden_inertia(6, actual)
        free_elements = set(range(1, 32)) - forbidden
        partitions = {min(sum(1 << i for i in t), 63 ^ sum(1 << i for i in t))
                      for t in combinations(range(6), 3)
                      if tuple(t) not in actual and tuple(i for i in range(6) if i not in t) not in actual}
        # Gauge last coordinate to zero; min(S, complement S) does just that.
        check(free_elements == partitions, 'six-root involution decoder')
        check(not any(not h & forbidden for h in subgroup_census(5, 2)), 'free order-four group on K3')
        six.append({'subset': subset, 'triples': sorted(actual), 'free_involutions': sorted(free_elements)})
        six_histogram[(len(actual), len(free_elements))] += 1
    for subset in combinations(range(9), 7):
        positions = {v: i for i, v in enumerate(subset)}
        actual = frozenset(tuple(positions[v] for v in t) for t in actual_global if set(t) <= set(subset))
        histogram = Counter(len(f & actual) for f in structures)
        forbidden = forbidden_inertia(7, actual)
        check(sum(not info['kernel'] & forbidden for info in structures.values()) == histogram[0],
              'seven-subset inertia/control mismatch')
        seven.append({'subset': subset, 'triples': sorted(actual), 'bad_triple_histogram': dict(sorted(histogram.items())),
                      'free_order_eight': histogram[0]})
        seven_histogram[(len(actual), histogram[0])] += 1
    check(dict(six_histogram) == {(0, 10): 2, (1, 9): 20, (2, 8): 40, (2, 9): 6, (3, 7): 16}, 'six census')
    check(dict(seven_histogram) == {(2, 20): 4, (3, 16): 10, (3, 17): 8, (4, 14): 12, (5, 12): 2}, 'seven census')
    return {'nine_line_triples': sorted(actual_global), 'six_subsets': six, 'seven_subsets': seven,
            'six_summary': [{'triples': t, 'free_involutions': f, 'subsets': count}
                            for (t, f), count in sorted(six_histogram.items())],
            'seven_summary': [{'triples': t, 'free_order_eight': f, 'subsets': count}
                              for (t, f), count in sorted(seven_histogram.items())]}


def paley_checks(structures):
    all_fanos = list(structures)
    pairs = {frozenset((a, b)) for a, b in combinations(all_fanos, 2) if not a & b}
    check(len(pairs) == 120, 'disjoint Fano pair count')
    for pair in pairs:
        a, b = tuple(pair)
        check(structures[a]['kernel'] & structures[b]['kernel'] == {0}, 'not complementary deck groups')
        check(len(span(tuple(structures[a]['kernel']) + tuple(structures[b]['kernel']))) == 64,
              'complementary groups do not span deck group')
    edge_positions = tuple(combinations(range(7), 2))
    paley_masks = set()
    for p in permutations(range(7)):
        paley_masks.add(sum(1 << index for index, (i, j) in enumerate(edge_positions)
                            if (p[j] - p[i]) % 7 in (1, 2, 4)))
    check(len(paley_masks) == 240, 'labelled Paley orbit count')
    pair_multiplicity = Counter()
    pair_masks = {}
    for mask in paley_masks:
        arcs = {(i, j) if mask >> index & 1 else (j, i) for index, (i, j) in enumerate(edge_positions)}
        triangles = frozenset(t for t in combinations(range(7), 3)
                              if all(sum((i, j) in arcs for j in t if j != i) == 1 for i in t))
        supported = frozenset(f for f in all_fanos if f <= triangles)
        check(len(triangles) == 14 and supported in pairs, 'Paley two-Fano decomposition')
        check(frozenset().union(*supported) == triangles, 'Paley triangle union')
        pair_multiplicity[supported] += 1
        pair_masks.setdefault(supported, []).append(mask)
    check(set(pair_multiplicity) == pairs and set(pair_multiplicity.values()) == {2}, 'Paley converse pair bijection')
    check(all(a ^ b == (1 << 21) - 1 for a, b in pair_masks.values()),
          'two Paley lifts are not orientation reversals')
    histogram = Counter(tuple(sorted(len(f & HOURGLASS_TRIPLES) for f in pair)) for pair in pairs)
    check(dict(histogram) == {(0, 0): 32, (0, 1): 48, (0, 3): 16, (1, 1): 24}, 'safe paired quotient census')
    return {'disjoint_fano_pairs': len(pairs), 'labelled_paley_tournaments': len(paley_masks),
            'paley_tournaments_per_fano_pair': 2, 'both_free_pairs': histogram[(0, 0)],
            'joint_bad_triples': [{'bad_counts': key, 'pairs': count} for key, count in sorted(histogram.items())]}


def reflection_check():
    """A support involution is not a binary-linear incidence automorphism."""
    original = (1, 2, 4, 3, 6, 5, 7)
    reflected = tuple(x ^ 7 if x != 7 else 7 for x in original)
    standard_lines = frozenset(t for t in combinations(range(1, 8), 3) if xor_all(t) == 0)
    j = lambda x: x ^ 7 if x != 7 else 7
    preserved = frozenset(t for t in standard_lines if xor_all(j(x) for x in t) == 0)
    check(len(preserved) == 3, 'modified complement preserves wrong Fano line count')
    before = sorted(t for t in HOURGLASS_TRIPLES if xor_all(original[i] for i in t) == 0)
    after = sorted(t for t in HOURGLASS_TRIPLES if xor_all(reflected[i] for i in t) == 0)
    check(not before and after == [(0, 1, 2)], 'modified complement free quotient hostile')
    return {'original_labels': original, 'reflected_labels': reflected,
            'preserved_standard_fano_lines': sorted(preserved),
            'bad_actual_triples_before': before, 'bad_actual_triples_after': after}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    structures, labels = fano_structures()
    seven = seven_checks(structures)
    subsets = subset_checks(structures)
    paley = paley_checks(structures)
    output = {'status': 'FINITE-EXACT exhaustive declared combinatorial universes',
              'line_order': LINES, 'fano_labeling_histogram': dict(sorted(labels.items())),
              'hourglass': seven, 'subsets': subsets, 'paley': paley,
              'modified_complement_hostile': reflection_check()}
    payload = json.dumps(output, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.write_text(payload)
    print(json.dumps({'fano_structures': len(structures), 'hourglass_free': seven['free_subgroups'],
                      'six_subsets': 84, 'seven_subsets': 36, 'both_free_paley_pairs': 32,
                      'sha256': hashlib.sha256(payload.encode()).hexdigest()}, sort_keys=True))


if __name__ == '__main__':
    main()
