#!/usr/bin/env python3
"""Norm-seven face systems select two complementary free hourglass kernel pairs.

Standalone exact integration check. No producer imports or optional packages.
Print the full JSON receipt by default; --output additionally saves that JSON.
"""
from collections import Counter
from itertools import combinations, permutations
from pathlib import Path
import argparse
import json

GATES = 0
BRANCH_VECTORS = ((0, 0), (1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (-1, -1))
TRIPLES = tuple(combinations(range(7), 3))
ACTUAL = frozenset(((0, 1, 2), (0, 3, 4), (0, 5, 6)))


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


def design(labels):
    return frozenset(t for t in TRIPLES if xor(labels[i] for i in t) == 0)


def kernel(labels):
    return frozenset(w for w in range(128) if w.bit_count() % 2 == 0
                     and xor(labels[i] for i in range(7) if w >> i & 1) == 0)


def face_system(k, face_type, residues):
    inverse = {value: i for i, value in enumerate(residues)}
    middle = 1 if face_type == 'U' else k
    return frozenset(tuple(sorted(inverse[value] for value in (t, (t + middle) % 7,
                                                              (t + 1 + k) % 7)))
                     for t in range(7))


def key(fano):
    return tuple(sorted(fano))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()

    representatives, fibers = {}, Counter()
    for labels in permutations(range(1, 8)):
        fano = design(labels)
        representatives.setdefault(fano, labels)
        fibers[fano] += 1
    check(len(representatives) == 30, 'all thirty labelled Fano structures')
    check(set(fibers.values()) == {168}, 'GL3 fibers of all 5040 markings')
    all_fanos = sorted(representatives, key=key)
    indices = {fano: i for i, fano in enumerate(all_fanos)}
    kernels = {}
    for fano in all_fanos:
        h = kernel(representatives[fano])
        # Independent decoder: nonzero even representatives are complements of Fano triples.
        decoded = frozenset({0} | {127 ^ sum(1 << i for i in t) for t in fano})
        check(h == decoded and len(h) == 8, 'Hamming/simplex kernel decoder')
        check(all(a ^ b in h for a in h for b in h), 'kernel is an additive subgroup')
        kernels[fano] = h
    free = [f for f in all_fanos if not f & ACTUAL]
    check(len(free) == 16, 'hourglass free Fano census')
    free_pairs = []
    for a, b in combinations(free, 2):
        check((not a & b) == (kernels[a] & kernels[b] == {0}), 'disjoint designs iff complementary kernels')
        if not a & b:
            free_pairs.append((indices[a], indices[b]))
    check(len(free_pairs) == 32, 'all thirty-two complementary free pairs')

    rows, selected_pairs, selected_fanos = [], [], []
    for k, ideal in ((2, (2, -1)), (4, (3, 1))):
        a, b = ideal
        check(a * a - a * b + b * b == 7, 'norm-seven ideal index')
        check((k * k + k + 1) % 7 == 0, 'residue map respects Eisenstein relation')
        check((a + k * b) % 7 == 0, 'ideal is killed by its residue map')
        residues = tuple((x + k * z) % 7 for x, z in BRANCH_VECTORS)
        check(set(residues) == set(range(7)), 'common branch labels biject onto residues')
        faces = [face_system(k, face_type, residues) for face_type in ('U', 'D')]
        check(not faces[0] & faces[1], 'each ideal has disjoint U/D face systems')
        for face_type, fano in zip(('U', 'D'), faces):
            pair_counts = Counter(p for t in fano for p in combinations(t, 2))
            check(len(fano) == 7 and len(pair_counts) == 21 and set(pair_counts.values()) == {1},
                  'face system is an STS(7)')
            check(fano in representatives, 'face system in independent Fano universe')
            check(not fano & ACTUAL, 'face system avoids all actual branch triples')
            h = kernels[fano]
            # Direct inertia audit in the even-representative model G7=E7.
            inertia = [127 ^ (1 << i) for i in range(7)]
            forbidden = set(inertia) | {inertia[i] ^ inertia[j] for i, j in combinations(range(7), 2)}
            forbidden |= {xor(inertia[i] for i in t) for t in ACTUAL}
            check(not h & forbidden, 'direct single/pair/triple inertia freeness')
            rows.append({'ideal_root': k, 'ideal_generator': ideal, 'face_type': face_type,
                         'branch_residues': residues, 'Fano_index': indices[fano],
                         'branch_triples': sorted(fano), 'binary_marking': representatives[fano],
                         'kernel_even_words': sorted(h)})
            selected_fanos.append(fano)
        h_u, h_d = (kernels[f] for f in faces)
        check(h_u & h_d == {0}, 'selected kernels have trivial intersection')
        check(len({u ^ d for u in h_u for d in h_d}) == 64, 'selected kernels sum to full deck group')
        pair = tuple(sorted(indices[f] for f in faces))
        check(pair in free_pairs, 'selected pair belongs to the thirty-two pair census')
        selected_pairs.append(pair)
    check(len(set(selected_fanos)) == 4, 'four distinct Fano structures on common branch labels')
    check(len({kernels[f] for f in selected_fanos}) == 4, 'four distinct embedded deck kernels')
    check(len(set(selected_pairs)) == 2, 'two distinct selected unordered kernel pairs')

    result = {'status': 'FINITE-EXACT norm-seven complementary-pair integration',
              'gates': GATES, 'branch_order': ['E', 'E+x', 'E-x', 'E+z', 'E-z', 'E+x+z', 'E-x-z'],
              'branch_vectors': BRANCH_VECTORS, 'actual_triples': sorted(ACTUAL),
              'labeling_universe': 5040, 'Fano_universe': 30, 'free_kernels': 16,
              'complementary_free_pair_count': 32, 'selected_pair_count': 2,
              'distinct_selected_kernel_count': 4, 'selected_pairs': selected_pairs,
              'all_complementary_free_pairs': free_pairs, 'constructions': rows,
              'geometric_consequence': 'Each selected pair recovers the normal root cover by normalization of the finite quotient fiber product.',
              'scope': 'No tournament orientation or chromatic assertion; no arithmetic point lifting.'}
    payload = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.write_text(payload)
    print(payload, end='')


if __name__ == '__main__':
    main()
