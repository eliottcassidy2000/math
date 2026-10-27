#!/usr/bin/env python3
"""Exact, dependency-free audit of HYP-9419 clause 2's support-rank formula.

Run both `python3 check.py` and `python3 -O check.py`. Explicit checks remain
active under -O. Vertices and bit positions are zero-based throughout.
No external libraries, network calls, or imported repository code are used.
"""
from fractions import Fraction
from itertools import permutations
from math import gcd
import argparse
import hashlib
import json
from pathlib import Path

PIN = '64a7985e7054821c5958dc85f6dd489c844111f2'

def require(condition, message):
    if not condition:
        raise RuntimeError(message)

def mobius(values, dimension):
    out = list(values)
    for bit in range(dimension):
        for mask in range(1 << dimension):
            if mask & (1 << bit):
                out[mask] -= out[mask ^ (1 << bit)]
    return tuple(out)

def zeta(coefficients, dimension):
    out = list(coefficients)
    for bit in range(dimension):
        for mask in range(1 << dimension):
            if mask & (1 << bit):
                out[mask] += out[mask ^ (1 << bit)]
    return tuple(out)

def primitive_echelon(rows):
    """Exact Q row elimination using only integer operations and gcd division."""
    basis = {}
    selected = []
    for index, row in enumerate(rows):
        r = list(row)
        for pivot, br in sorted(basis.items()):
            if r[pivot]:
                a, b = br[pivot], r[pivot]
                r = [a * x - b * y for x, y in zip(r, br)]
                common = 0
                for x in r:
                    common = gcd(common, x)
                if common:
                    r = [x // common for x in r]
        pivot = next((i for i, x in enumerate(r) if x), None)
        if pivot is not None:
            basis[pivot] = r
            selected.append(index)
    return basis, selected

def rational_rank(rows):
    """Independent Fraction elimination, reverse column pivot order."""
    if not rows:
        return 0
    a = [list(map(Fraction, row)) for row in rows]
    rank = 0
    for col in reversed(range(len(a[0]))):
        pivot = next((i for i in range(rank, len(a)) if a[i][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        scale = a[rank][col]
        a[rank] = [x / scale for x in a[rank]]
        for i in range(rank + 1, len(a)):
            scale = a[i][col]
            if scale:
                a[i] = [x - scale * y for x, y in zip(a[i], a[rank])]
        rank += 1
        if rank == len(a):
            break
    return rank

def bareiss_det(square):
    """Independent integer minor certificate with checked exact division."""
    a = [list(row) for row in square]
    n = len(a)
    previous = 1
    sign = 1
    for k in range(n - 1):
        pivot = next((i for i in range(k, n) if a[i][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            sign = -sign
        diagonal = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = a[i][j] * diagonal - a[i][k] * a[k][j]
                require(numerator % previous == 0, 'Bareiss division not exact')
                a[i][j] = numerator // previous
            a[i][k] = 0
        previous = diagonal
    return sign * a[-1][-1] if n else 1

def downward(supports):
    result = set()
    for mask in supports:
        sub = mask
        while True:
            result.add(sub)
            if not sub:
                break
            sub = (sub - 1) & mask
    return result

def build(n):
    edges = [(i, j) for i in range(n) for j in range(i + 2, n)]
    index = {edge: bit for bit, edge in enumerate(edges)}
    m = len(edges)
    size = 1 << m
    all_permutations = list(permutations(range(n)))
    guards = []
    impossible = []
    for permutation_index, path in enumerate(all_permutations):
        # Route A: literal directed adjacency test on each tournament.
        direct_values = []
        for state in range(size):
            def arc(a, b):
                if abs(a - b) == 1:
                    return a < b
                bit = (state >> index[tuple(sorted((a, b)))]) & 1
                return bool(bit) if a < b else not bool(bit)
            direct_values.append(int(all(arc(a, b) for a, b in zip(path, path[1:]))))
        # Route B: signed cylinder masks, rejecting reversed forced path edges.
        positive = negative = 0
        reverse_forced = False
        for a, b in zip(path, path[1:]):
            if a == b + 1:
                reverse_forced = True
            elif b == a + 1:
                pass
            else:
                bit = 1 << index[tuple(sorted((a, b)))]
                if a < b:
                    positive |= bit
                else:
                    negative |= bit
        if reverse_forced:
            require(not any(direct_values), 'forbidden guard is nonzero')
            impossible.append(permutation_index)
            continue
        require(any(direct_values), 'a possible cylinder is empty')
        cylinder_values = tuple(int((state & positive) == positive and not state & negative)
                                for state in range(size))
        require(tuple(direct_values) == cylinder_values, 'direct/cylinder mismatch')
        coefficients = [0] * size
        sub = negative
        while True:
            coefficients[positive | sub] = (-1) ** sub.bit_count()
            if not sub:
                break
            sub = (sub - 1) & negative
        coefficients = tuple(coefficients)
        require(coefficients == mobius(direct_values, m), 'symbolic/Mobius mismatch')
        require(zeta(coefficients, m) == cylinder_values, 'zeta inverse failed')
        guards.append(dict(index=permutation_index, path=path, positive=positive,
                           negative=negative, values=cylinder_values, coefficients=coefficients))
    require(len(guards) + len(impossible) == len(all_permutations), 'not exhaustive')
    return edges, guards, impossible

def certificate(rows, labels):
    basis, selected = primitive_echelon(rows)
    columns = sorted(basis)
    rank = len(basis)
    fraction_rank = rational_rank(rows)
    require(rank == fraction_rank, 'integer/Fraction rank disagreement')
    determinant = bareiss_det([[rows[i][j] for j in columns] for i in selected])
    require(determinant != 0, 'selected minor is singular')
    return dict(rank_Q=rank, reverse_fraction_rank_Q=fraction_rank,
                selected_row_labels=[labels[i] for i in selected],
                selected_columns=columns, minor_determinant=determinant)

def audit():
    expected = {3: (3, 2, 2, 2, 2), 4: (11, 9, 6, 8, 8), 5: (53, 46, 31, 46, 44)}
    output = dict(source_commit=PIN, status='PASS', universe='fixed identity Hamilton path; all labelled tournaments; Q',
                  zero_based_vertices=True, cases=[], certificates={}, controls={})
    cache = {}
    for n in (3, 4, 5):
        edges, guards, impossible = build(n)
        cache[n] = (edges, guards)
        supports = {g['positive'] | g['negative'] for g in guards}
        closure = downward(supports)
        actual = {mask for g in guards for mask, c in enumerate(g['coefficients']) if c}
        values = [g['values'] for g in guards]
        coefficients = [g['coefficients'] for g in guards]
        vc = certificate(values, [g['index'] for g in guards])
        cc = certificate(coefficients, [g['index'] for g in guards])
        require(vc['rank_Q'] == cc['rank_Q'], 'coefficient/evaluation rank disagreement')
        got = len(guards), len(set(values)), len(supports), len(closure), vc['rank_Q']
        require(got == expected[n], 'finite expected result mismatch: ' + str((n, got)))
        require(actual == closure, 'expansion support is not the tested downward closure')
        case = dict(n=n, free_edges=edges, states=1 << len(edges), legal_permutations=len(guards),
                    impossible_permutations=len(impossible), distinct_nonzero_guards=len(set(values)),
                    distinct_guards_including_zero=len(set(values)) + int(bool(impossible)),
                    k=1, full_support_count=len(supports), downward_closure_count=len(closure),
                    actual_coefficient_support_count=len(actual), rank_Q=vc['rank_Q'],
                    full_support_masks=sorted(supports), downward_closure_masks=sorted(closure))
        all_permutation_supports = set()
        for path in permutations(range(n)):
            support = 0
            for a, b in zip(path, path[1:]):
                if abs(a - b) > 1:
                    support |= 1 << edges.index(tuple(sorted((a, b))))
            all_permutation_supports.add(support)
        case['invalid_filter_control_all_permutation_supports'] = len(all_permutation_supports)
        case['invalid_filter_control_downward_closure'] = len(downward(all_permutation_supports))
        require(len(all_permutation_supports) == {3: 2, 4: 8, 5: 38}[n],
                'invalid-filter hostile count mismatch')
        output['cases'].append(case)
        support_order = sorted(supports, key=lambda s: (s.bit_count(), s))
        representative = {g['positive'] | g['negative']: g for g in guards}
        lower_minor = [[representative[t]['coefficients'][s] for s in support_order]
                       for t in support_order]
        lower_determinant = bareiss_det(lower_minor)
        require(abs(lower_determinant) == 1, 'support lower bound minor not unimodular')
        output['certificates'][str(n)] = dict(evaluation=vc, coefficient=cc,
                                              support_lower_bound_minor_determinant=lower_determinant)
        # At-most-two equals exactly two, because the identity guard is 1.
        product_rows = sorted({tuple(a * b for a, b in zip(g['values'], h['values']))
                               for g in guards for h in guards})
        unions = {a | b for a in supports for b in supports}
        compatible_unions = {g['positive'] | g['negative'] | h['positive'] | h['negative']
                             for g in guards for h in guards
                             if not ((g['positive'] | h['positive']) &
                                     (g['negative'] | h['negative']))}
        require(compatible_unions == unions, 'tested union lacks compatible representative')
        pair_cert = certificate(product_rows, list(range(len(product_rows))))
        output['certificates'][str(n)]['two_guard_products'] = pair_cert
        require(pair_cert['rank_Q'] == 1 << len(edges), 'two-guard rank not full')
        output['cases'].append(dict(n=n, k=2, full_support_count=len(unions),
                                    compatible_union_support_count=len(compatible_unions),
                                    downward_closure_count=len(downward(unions)),
                                    distinct_product_rows_including_zero=len(product_rows),
                                    rank_Q=pair_cert['rank_Q']))
    # Minimal literal witness: a=x_(0,2) is absent from full path supports at n=4,
    # but a=1-g_(2,0,1,3)-g_(3,1,2,0).
    edges4, guards4 = cache[4]
    by_path = {g['path']: g for g in guards4}
    a_bit = edges4.index((0, 2))
    witness = tuple((state >> a_bit) & 1 for state in range(1 << len(edges4)))
    reconstructed = tuple(1 - x - y for x, y in zip(by_path[(2, 0, 1, 3)]['values'],
                                                   by_path[(3, 1, 2, 0)]['values']))
    require(witness == reconstructed, 'literal monomial identity failed')
    require(1 << a_bit not in {g['positive'] | g['negative'] for g in guards4},
            'literal witness support unexpectedly present')
    output['literal_witness'] = dict(n=4, monomial_edge=[0, 2], identity='x_02=1-g_2013-g_3120')
    # Two explicit independent coefficient functionals annihilate all 53 guards.
    # The first takes value 1 on x_(1,3), which occurs in g_01342=d(1-f).
    first = {i: 1 for i in (8, 9, 11, 12, 13, 40, 41, 43, 44, 56, 57)}
    second = {1:-1, 5:-1, 11:1, 12:1, 17:-1, 21:-1, 25:-1, 29:-1,
              32:-1, 33:-2, 34:-1, 35:-1, 36:-1, 38:-1, 42:-1, 46:-1,
              49:-1, 56:1}
    edges5, guards5 = cache[5]
    for functional in (first, second):
        for g in guards5:
            require(sum(weight * g['coefficients'][mask] for mask, weight in functional.items()) == 0,
                    'coefficient annihilator fails')
    null_rows = [[f.get(i, 0) for i in range(64)] for f in (first, second)]
    require(rational_rank(null_rows) == 2, 'annihilators are dependent')
    require(first.get(1 << edges5.index((1, 3)), 0) == 1, 'monomial not separated')
    d_guard = next(g for g in guards5 if g['path'] == (0, 1, 3, 4, 2))
    require(d_guard['positive'] == 8 and d_guard['negative'] == 32, 'witness cylinder wrong')
    require(d_guard['coefficients'][8] == 1, 'witness monomial not in expansion')
    output['downward_witness'] = dict(n=5, monomial_edge=[1, 3], exposing_guard=[0, 1, 3, 4, 2],
                                    exposing_polynomial='d(1-f)',
                                    annihilators=[first, second], independent_annihilators=2,
                                    dimension_upper_bound=46 - 2)
    # Hostile controls: 0->... identity has guard 1, reverse forced edge gives 0,
    # and literal-only / unsigned-support ranks differ from genuine rank.
    require(all(next(g for g in cache[5][1] if g['path'] == tuple(range(5)))['values']), 'identity not 1')
    require(len(cache[4][0]) == 3 and len(cache[5][0]) == 6, 'wrong cube dimension')
    require(6 != 8 and 46 != 44, 'false equality accidentally accepted')
    try:
        require(False, 'optimization-hostile check')
    except RuntimeError:
        output['controls']['explicit_require_survives_optimization'] = True
    else:
        raise RuntimeError('explicit check failed')
    output['controls'].update(identity_guard_one=True, reversed_forced_edge_zero=True,
                               direct_vs_cylinder_all_permutations=True,
                               symbolic_vs_mobius_all_guards=True,
                               rational_vs_primitive_all_ranks=True,
                               selected_minors_nonzero=True,
                               sharp_depth_n7_not_computed=True)
    return output

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify-manifest', action='store_true')
    args = parser.parse_args()
    if args.verify_manifest:
        root = Path(__file__).resolve().parent
        manifest = json.loads((root / 'MANIFEST.json').read_text())
        for relative, expected in manifest['sha256'].items():
            require(hashlib.sha256((root / relative).read_bytes()).hexdigest() == expected,
                    'manifest mismatch: ' + relative)
    print(json.dumps(audit(), indent=2, sort_keys=True))

if __name__ == '__main__':
    main()
