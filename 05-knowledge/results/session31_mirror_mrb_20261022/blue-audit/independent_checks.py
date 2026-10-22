#!/usr/bin/env python3
"""Independent small exact audit. No imports from the submission.

Determinants: signed subset assignment DP on every principal minor.
Forests: undirected union-find rejection while assigning outgoing arcs.
All-presentation enumeration stops at n=6 (1,024 presentations).
"""
from collections import defaultdict
from fractions import Fraction as Q
from itertools import permutations, product, combinations
from math import comb
from pathlib import Path
import argparse
import hashlib
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def pairs(n):
    return [(i, j) for i in range(n) for j in range(i + 2, n)]


def orbits(n):
    grouped = defaultdict(list)
    for i, j in pairs(n):
        grouped[min((i, j), (n - 1 - j, n - 1 - i))].append((i, j))
    return list(grouped.values())


def tournament(n, groups, mask):
    a = [[0] * n for _ in range(n)]
    for i in range(n - 1):
        a[i][i + 1] = 1
    for k, group in enumerate(groups):
        for i, j in group:
            a[i][j] = (mask >> k) & 1
            a[j][i] = 1 - a[i][j]
    return a


def transpose(a):
    return [list(row) for row in zip(*a)]


def complement(a):
    n = len(a)
    return [[a[i][j] if abs(i - j) <= 1 else a[j][i]
             for j in range(n)] for i in range(n)]


def reflected(a):
    n = len(a)
    return [[a[n - 1 - j][n - 1 - i] for j in range(n)]
            for i in range(n)]


def laplacian(a):
    return [[sum(row) if i == j else -x for j, x in enumerate(row)]
            for i, row in enumerate(a)]


def principal_det(a, ids):
    """Leibniz sum by a column-subset DP, no elimination or trace identities."""
    m = len(ids)
    dp = [0] * (1 << m)
    dp[0] = 1
    for mask in range((1 << m) - 1):
        row = mask.bit_count()
        if not dp[mask]:
            continue
        for col in range(m):
            if mask & (1 << col):
                continue
            sign = -1 if (mask >> (col + 1)).bit_count() % 2 else 1
            dp[mask | (1 << col)] += sign * dp[mask] * a[ids[row]][ids[col]]
    return dp[-1]


def forest_polynomial(a):
    l = laplacian(a)
    n = len(l)
    out = [sum(principal_det(l, ids) for ids in combinations(range(n), k))
           for k in range(n + 1)]
    assert out[-1] == 0
    return out[:-1]


def energy(a):
    return sum((2 * sum(row) - len(a) + 1) ** 2 for row in a)


def triangles(a):
    return sum(a[i][j] * a[j][k] * a[k][i] +
               a[j][i] * a[k][j] * a[i][k]
               for i, j, k in combinations(range(len(a)), 3))


def all_mean_formula(n):
    p = [Q(0)] * n
    for ell in range((n - 1) // 2 + 1):
        for j in range(n - 2 * ell):
            p[j + 2 * ell] += Q(comb(n - 1 - ell, ell) *
                comb(n - 1 - 2 * ell, j) * n ** j, 4 ** ell * 2 ** j)
    return p


def forest_occurrences(n):
    """Assign root or one outgoing neighbor; union-find excludes any cycle.

    Count occurrence totals across cubes, with incompatible orbit pins dropped
    only for blue. This does not call a determinant or enumerate tournaments.
    """
    gs = orbits(n)
    oid = {e: k for k, g in enumerate(gs) for e in g}
    D, o = len(pairs(n)), len(gs)
    totals_all, totals_blue = [0] * n, [0] * n
    count = conflicts = 0

    def visit(i, components, pins, consistent, degree, free_edges):
        nonlocal count, conflicts
        if i == n:
            count += 1
            totals_all[degree] += 1 << (D - free_edges)
            if consistent:
                totals_blue[degree] += 1 << (o - len(pins))
            else:
                conflicts += 1
            return
        visit(i + 1, components, pins, consistent, degree, free_edges)
        for j in range(n):
            if i == j or j == i - 1 or components[i] == components[j]:
                continue
            src, dst = components[i], components[j]
            joined = tuple(dst if c == src else c for c in components)
            if j == i + 1:
                visit(i + 1, joined, pins, consistent, degree + 1, free_edges)
            else:
                key = oid[tuple(sorted((i, j)))]
                bit = int(i < j)
                okay = consistent and (key not in pins or pins[key] == bit)
                new_pins = dict(pins)
                new_pins[key] = bit
                visit(i + 1, joined, new_pins, okay, degree + 1, free_edges + 1)

    visit(0, tuple(range(n)), {}, True, 0, 0)
    return totals_all, totals_blue, count, conflicts


def isomorphisms(a, b):
    n = len(a)
    sa, sb = list(map(sum, a)), list(map(sum, b))
    return [p for p in permutations(range(n))
            if all(sa[i] == sb[p[i]] for i in range(n)) and
            all(a[i][j] == b[p[i]][p[j]] for i in range(n) for j in range(n))]


def code(a):
    n = len(a)
    return min(tuple(a[p[i]][p[j]] for i in range(n) for j in range(i + 1, n))
               for p in permutations(range(n)))


def main():
    if not __debug__:
        raise SystemExit('Run without -O: this audit deliberately uses checked assertions.')
    parser = argparse.ArgumentParser(description=__doc__)
    default_source = next((p for p in (ROOT / 'session31-blue-forests', ROOT / 'blue-forests')
                           if p.is_dir()), ROOT / 'session31-blue-forests')
    parser.add_argument('--source-dir', type=Path, default=default_source)
    source = parser.parse_args().source_dir
    receipt = {'status': 'PASS', 'method': __doc__, 'cube_checks': [],
               'orbit_gram_checks': [], 'forest_occurrence_checks': [],
               'merged_inventory_checks': [], 'zero_flux_checks': []}
    # Algebraic incidence-vector evaluation, without presentation enumeration.
    for n in range(2, 32):
        gs = orbits(n)
        fixed = sum(len(g) == 1 for g in gs)
        assert fixed == (n - 1) // 2
        assert len(gs) == (n - 1) ** 2 // 4
        gram = 2
        shared = 0
        for group in gs:
            v = [0] * n
            for i, j in group:
                v[i] += 1
                v[j] -= 1
            gram += sum(x * x for x in v)
            shared += int(len(group) == 2 and len({i for e in group for i in e}) == 3)
        delta = n - 3 if n % 2 else 0
        assert gram == n * n - 3 * n + 4 - delta
        assert shared == ((n - 3) // 2 if n % 2 else 0)
        receipt['orbit_gram_checks'].append({'n': n, 'blue_energy_mean': gram,
            'center_sharing_orbits': shared})

    for n in range(2, 7):
        edges, gs = pairs(n), orbits(n)
        D, o = len(edges), len(gs)
        M, B = 1 << D, 1 << o
        all_sum, blue_sum = [0] * n, [0] * n
        ea = eb = 0
        mats, blues = [], []
        # Independent free-edge cube followed by a matrix-level reflection test.
        for mask in range(M):
            a = tournament(n, [[e] for e in edges], mask)
            f, w, c = forest_polynomial(a), energy(a), triangles(a)
            ac = complement(a)
            isblue = a == reflected(a)
            mats.append(a)
            blues.append(isblue)
            ea += w
            all_sum = [x + y for x, y in zip(all_sum, f)]
            if isblue:
                eb += w
                blue_sum = [x + y for x, y in zip(blue_sum, f)]
            if n >= 3:
                N = comb(n, 2)
                assert 8 * f[2] == 4 * N * N - n * (n - 1) ** 2 - w
                assert f[2] == c + comb(N, 2) - comb(n, 3)
            flux = sum(a[0]) - sum(a[-1]) - 1
            assert w - energy(ac) == 8 * flux
            assert triangles(ac) - c == flux
            assert w == n * (n * n - 1) // 3 - 8 * c
            residue = n if n % 2 == 0 else 0
            assert (w - residue) % 8 == 0
            if isblue:
                assert flux == 2 * sum(a[0]) - n
                assert flux % 2 == n % 2
                K = [[2 * laplacian(a)[i][j] - n * int(i == j) + 1
                      for j in range(n)] for i in range(n)]
                assert all(sum(row) == 0 for row in K)
                assert all(K[n - 1 - i][n - 1 - j] == -K[i][j]
                           for i in range(n) for j in range(n))
        assert sum(blues) == B
        mu, delta = n * n - 3 * n + 4, (n - 3 if n % 2 else 0)
        assert ea == M * mu and eb == B * (mu - delta)
        p, q = [Q(x, M) for x in all_sum], [Q(x, B) for x in blue_sum]
        assert p == all_mean_formula(n)
        diff = [x - y for x, y in zip(q, p)]
        if n >= 3:
            assert diff[2] == Q(delta, 8)
        if n >= 4:
            assert diff[3] == Q(delta * n * (n - 3), 16)
            assert Q(ea - eb, M - B) == mu + Q(B * delta, M - B)
        if n >= 3:
            position = {e: k for k, e in enumerate(edges)}
            sigma = lambda mask: sum(((mask >> k) & 1) << position[(n - 1 - j, n - 1 - i)]
                                     for k, (i, j) in enumerate(edges))
            seen, sizes = set(), defaultdict(int)
            for mask in range(M):
                if mask in seen:
                    continue
                orbit = {mask, mask ^ (M - 1), sigma(mask), sigma(mask) ^ (M - 1)}
                sizes[len(orbit)] += 1
                seen |= orbit
            assert dict(sizes) == ({2: B // 2} if M == B else {2: B // 2, 4: (M - B) // 4})
        row = {'n': n, 'all_presentations': M, 'blue_presentations': B,
            'all_energy_sum': ea, 'blue_energy_sum': eb,
            'all_polynomial_sum': all_sum, 'blue_polynomial_sum': blue_sum,
            'blue_minus_all': list(map(str, diff))}
        receipt['cube_checks'].append(row)
        af, bf, cnt, bad = forest_occurrences(n)
        assert af == all_sum and bf == blue_sum
        receipt['forest_occurrence_checks'].append({'n': n,
            'path_compatible_forests': cnt, 'blue_orbit_conflicts': bad,
            'all_occurrence_total': af, 'blue_occurrence_total': bf})

        if n in (4, 6):
            candidates = [a for a, blue in zip(mats, blues)
                          if blue and sum(a[0]) == n // 2]
            loops = [a for a in candidates if isomorphisms(a, complement(a))]
            assert (len(candidates), len(loops)) == ((2, 2) if n == 4 else (24, 4))
            receipt['zero_flux_checks'].append({'n': n,
                'candidate_endpoints': len(candidates), 'loop_endpoints': len(loops)})

        if n <= 5:
            classes = {}
            for a, blue in zip(mats, blues):
                key = code(a)
                if key not in classes:
                    classes[key] = {'adjacency': a, 'h': 0, 'b': 0,
                        'converse': code(transpose(a)), 'F': forest_polynomial(a)}
                classes[key]['h'] += 1
                classes[key]['b'] += int(blue)
            seen, bs, ks = set(), [Q(0)] * n, [Q(0)] * n
            nonconverse_polynomial_count = 0
            for key, c in classes.items():
                if key in seen:
                    continue
                mate = classes[c['converse']]
                assert c['h'] == mate['h']
                sc = key == c['converse']
                seen |= {key, c['converse']}
                w = c['h'] if sc else 2 * c['h']
                b = c['b']
                assert sc or b == 0
                assert (w - b) % 2 == 0
                fbar = c['F'] if sc else [Q(x + y, 2) for x, y in zip(c['F'], mate['F'])]
                nonconverse_polynomial_count += int(c['F'] != mate['F'])
                k = (w - b) // 2
                bs = [x + b * y for x, y in zip(bs, fbar)]
                ks = [x + k * y for x, y in zip(ks, fbar)]
            assert bs == blue_sum
            assert ks == [Q(x - y, 2) for x, y in zip(all_sum, blue_sum)]
            receipt['merged_inventory_checks'].append({'n': n,
                'ordinary_classes': len(classes),
                'distinct_converse_pairs_with_different_polynomials': nonconverse_polynomial_count,
                'blue_polynomial_mass': list(map(str, bs)),
                'fundamental_black_incidence_polynomial_mass': list(map(str, ks))})

    diffs = {r['n']: r['blue_minus_all'] for r in receipt['cube_checks']}
    assert diffs[5] == ['0', '0', '1/4', '5/4', '27/16']
    assert diffs[6] == ['0', '0', '0', '0', '5/16', '15/16']

    # Literal adjacency in the note, independent of orbit-mask construction.
    outgoing = ((2, 3, 4), (3,), (4, 6), (2, 5, 6), (1, 2, 3, 6), (1, 2))
    a = [[int(j + 1 in outgoing[i]) for j in range(6)] for i in range(6)]
    ac = complement(a)
    assert all(a[i][j] + a[j][i] == 1 for i in range(6) for j in range(i + 1, 6))
    assert all(a[i][i + 1] for i in range(5))
    assert a == reflected(a) and ac == reflected(ac)
    fa, fc = forest_polynomial(a), forest_polynomial(ac)
    assert fa == [1, 15, 91, 279, 433, 273]
    assert fc == [1, 15, 91, 279, 432, 270]
    assert energy(a) == energy(ac) == 22
    assert triangles(a) == triangles(ac) == 6
    assert sorted(map(sum, a)) == sorted(map(sum, ac)) == [1, 2, 2, 3, 3, 4]
    assert not isomorphisms(a, ac)
    assert not isomorphisms(a, transpose(ac))
    receipt['witness'] = {'adjacency': a, 'complement_adjacency': ac,
        'forest_polynomial': fa, 'complement_forest_polynomial': fc,
        'outdegrees': list(map(sum, a)), 'complement_outdegrees': list(map(sum, ac)),
        'energy': 22, 'cyclic_triangles': 6, 'isomorphism_count': 0,
        'converse_isomorphism_count': 0,
        'self_converse_permutation': [5, 4, 3, 2, 1, 0]}

    source_files = ['BLUE_FOREST_CONSTRAINTS.md', 'check_blue_forests.py',
                    'check_wiring.py', 'probe.py', 'verification.json', 'wiring_verification.json']
    receipt['source_sha256'] = {name: hashlib.sha256(
        (source / name).read_bytes()).hexdigest() for name in source_files}
    submission = json.loads((source / 'verification.json').read_text())
    known = {row['n']: row for row in submission['full_polynomials']}
    for row in receipt['cube_checks']:
        assert row['blue_polynomial_sum'] == known[row['n']]['blue_total_coefficients']
        assert row['blue_minus_all'] == known[row['n']]['blue_minus_all_coefficients']
    (HERE / 'independent_receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print('PASS: exact cubes n=2..6; Gram formulas n=2..31; independent forest occurrence sums n=2..6; merged inventory n=2..5; minimal n=6 witness.')


if __name__ == '__main__':
    main()
