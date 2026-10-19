#!/usr/bin/env python3
"""Finite-exact p=29 certificates; Python standard library only.

This enumerates 1,134 supports of 12 points, not sections or colorings.
Every elimination uses integer arithmetic in F_29. No solver is used.
"""
from collections import Counter
from itertools import combinations
from pathlib import Path
import json

P = 29
DIRECTIONS = ((1, 0), (0, 1), (3, 9))
SIGMA = (-3, -9, 1)
PAIRS = tuple(tuple(DIRECTIONS[j] for j in range(3) if j != i) for i in range(3))


def add(x, y):
    return ((x[0] + y[0]) % P, (x[1] + y[1]) % P)


def scale(a, x):
    return ((a * x[0]) % P, (a * x[1]) % P)


def cross(x, y):
    return (x[0] * y[1] - x[1] * y[0]) % P


def derivative(f, d):
    """Sparse INTEGER derivative D_d f(x)=f(x+d)-f(x)."""
    out = {}
    for x, value in f.items():
        for y, sign in ((add(x, scale(-1, d)), 1), (x, -1)):
            out[y] = out.get(y, 0) + sign * value
    return {x: value for x, value in out.items() if value}


def triple(f):
    for d in DIRECTIONS:
        f = derivative(f, d)
    return f


def translate(f, t):
    return {add(x, t): value for x, value in f.items()}


def dipole_cube(e, q):
    f = {(0, 0): 1, scale(q, DIRECTIONS[e]): -1}
    for d in PAIRS[e]:
        f = derivative(f, d)
    return f


def rref(matrix, ncols=None):
    """Gauss-Jordan elimination, with all pivot inverses in F_29."""
    a = [[value % P for value in row] for row in matrix]
    if ncols is None:
        ncols = len(a[0]) if a else 0
    assert all(len(row) == ncols for row in a)
    pivots = []
    for col in range(ncols):
        rank = len(pivots)
        pivot = next((i for i in range(rank, len(a)) if a[i][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        inv = pow(a[rank][col], -1, P)
        a[rank] = [value * inv % P for value in a[rank]]
        for i in range(len(a)):
            if i != rank and a[i][col]:
                factor = a[i][col]
                a[i] = [(x - factor * y) % P for x, y in zip(a[i], a[rank])]
        pivots.append(col)
        if len(pivots) == len(a):
            break
    return a[:len(pivots)], pivots


def nullspace(matrix, ncols):
    reduced, pivots = rref(matrix, ncols)
    basis = []
    for free in range(ncols):
        if free in pivots:
            continue
        vector = [0] * ncols
        vector[free] = 1
        for row, pivot in zip(reduced, pivots):
            vector[pivot] = -row[free] % P
        assert all(sum(x * y for x, y in zip(row, vector)) % P == 0 for row in matrix)
        basis.append(vector)
    return basis


def line_constraints(support, missing):
    rows = []
    for d in PAIRS[missing]:
        for label in sorted({cross(d, x) for x in support}):
            rows.append([int(cross(d, x) == label) for x in support])
    return rows


def incidence_certificate(cubes, expected_dimension, direct=True):
    support = sorted(set().union(*(set(c) for c in cubes)))
    n = len(support)
    lines = [line_constraints(support, i) for i in range(3)]
    spaces = [nullspace(rows, n) for rows in lines]
    matrix = [[SIGMA[i] * vector[r] % P
               for i in range(3) for vector in spaces[i]] for r in range(n)]
    coupled_rank = len(rref(matrix)[1])
    dimension = sum(map(len, spaces)) - coupled_rank
    assert dimension == expected_dimension

    # Explicit model vectors: (C,0,3C) and (0,C,9C), one pair per cube.
    model = []
    for cube in cubes:
        for projection in ((1, 0, 3), (0, 1, 9)):
            vector = [projection[i] * cube.get(x, 0) % P
                      for i in range(3) for x in support]
            for i in range(3):
                assert all(sum(row[r] * vector[i*n+r] for r in range(n)) % P == 0
                           for row in lines[i])
            assert all(sum(SIGMA[i] * vector[i*n+r] for i in range(3)) % P == 0
                       for r in range(n))
            model.append(vector)
    model_rank = len(rref(model)[1])
    assert model_rank == 2 * len(cubes)

    # A second presentation: all 3n entries at once, impose line constraints
    # and pointwise coupling. Its nullity independently agrees with the
    # restricted-domain calculation above.
    full_rank = None
    if direct:
        full = []
        for i in range(3):
            for row in lines[i]:
                full.append([0] * (i*n) + row + [0] * ((2-i)*n))
        for r in range(n):
            full.append([SIGMA[i] if j == r else 0 for i in range(3) for j in range(n)])
        full_rank = len(rref(full)[1])
        assert full_rank == 3*n - dimension
    return {'support_size': n, 'pair_kernel_dimensions': list(map(len, spaces)),
            'coupled_rank': coupled_rank, 'homogeneous_dimension': dimension,
            'direct_constraint_rank': full_rank, 'model_rank': model_rank}


def largest_gap_paths(e, q, r):
    points = (0, q, r)
    gaps = [((points[(i+1) % 3] - points[i]) % P, points[(i+1) % 3])
            for i in range(3)]
    gap, origin = max(gaps)
    lengths = sorted((x-origin) % P for x in points if x != origin)
    assert gap >= 10 and lengths[-1] == P-gap <= 19
    paths = [{scale(origin+j, DIRECTIONS[e]): 1 for j in range(1, k+1)}
             for k in lengths]
    bounds = []
    for i, d in enumerate(DIRECTIONS):
        changes = [derivative(path, d) for path in paths]
        for change in changes:
            assert sum(change.values()) == 0
            assert all(sum(value for x, value in change.items() if cross(d, x) == t) == 0
                       for t in range(P))
        support = set().union(*(set(change) for change in changes))
        column_direction = DIRECTIONS[e] if i != e else PAIRS[i][0]
        counts = Counter(cross(column_direction, x) for x in support)
        assert len(counts) <= (2 if i != e else 3)
        assert max(counts.values()) <= (19 if i != e else 1)
        assert P > 2*len(counts) and max(counts.values()) <= P-2
        bounds.append({'columns': len(counts), 'maximum_changes': max(counts.values())})
    return {'origin_coordinate': origin, 'largest_gap': gap, 'path_lengths': lengths,
            'normalization_support_bounds': bounds}


def sparse_record(f):
    return [[x[0], x[1], f[x]] for x in sorted(f)]


single_cube_certificates = []
for e in range(3):
    for k in range(1, 15):
        certificate = incidence_certificate([dipole_cube(e, k)], 2)
        assert certificate['support_size'] == 8
        assert certificate['pair_kernel_dimensions'] == [2, 2, 2]
        single_cube_certificates.append({'exceptional_direction': e, 'short_separation': k,
                                         **certificate})

records = []
histogram = Counter()
for e in range(3):
    for q, r in combinations(range(1, P), 2):
        cubes = [dipole_cube(e, q), dipole_cube(e, r)]
        assert all(len(cube) == 8 and set(cube.values()) == {-1, 1} for cube in cubes)
        assert len(set(cubes[0]) & set(cubes[1])) == 4
        assert all(cubes[0][x] == cubes[1][x] for x in set(cubes[0]) & set(cubes[1]))
        certificate = incidence_certificate(cubes, 4)
        assert certificate['support_size'] == 12 and certificate['model_rank'] == 4
        key = (e, tuple(certificate['pair_kernel_dimensions']), certificate['coupled_rank'])
        histogram[key] += 1
        records.append({'exceptional_direction': e, 'q': q, 'r': r,
                        **certificate, **largest_gap_paths(e, q, r)})
assert len(records) == 3 * 28 * 27 // 2 == 1134

# Disjoint cubes do NOT by themselves force the four-dimensional model.
# C_j = D_u D_v D_w delta_(j v), j=0,1,2.
cubes = [triple({scale(j, DIRECTIONS[1]): 1}) for j in range(3)]
assert set(cubes[0]).isdisjoint(cubes[2])
assert set(cubes[1]) <= set(cubes[0]) | set(cubes[2])
outer_support = set(cubes[0]) | set(cubes[2])
assert len(outer_support) == 16
hostile = incidence_certificate([cubes[0], cubes[2]], 6)
assert hostile['pair_kernel_dimensions'] == [6, 4, 6]
assert hostile['coupled_rank'] == 10 and hostile['direct_constraint_rank'] == 42
# The intermediate cube supplies precisely the two missing modes.
complete = incidence_certificate(cubes, 6)
assert complete['model_rank'] == 6
hostile['cubes_C0_C1_C2'] = [sparse_record(cube) for cube in cubes]
hostile['full_six_mode_basis_verified'] = True
hostile['is_binary_error_section'] = False
hostile['interpretation'] = 'Linear incidence hostile only; intermediate cube adds two local-motion modes.'

# A nonvacuous separated configuration satisfying the four-mode criterion.
passing_cubes = [dipole_cube(0, 1), translate(dipole_cube(1, 1), (7, 5))]
assert set(passing_cubes[0]).isdisjoint(passing_cubes[1])
passing = incidence_certificate(passing_cubes, 4)
assert passing['pair_kernel_dimensions'] == [4, 4, 4]
passing['dipoles'] = [{'direction': 0, 'P': [0, 0], 'Q': [1, 0]},
                     {'direction': 1, 'P': [7, 5], 'Q': [7, 6]}]

# Retain the sharp support warning for the normalization lemma: at p-1
# changes/column, weight and total-sum preservation no longer suffice.
V = [(s, t) for s in range(P) for t in range(P)]
J = {(s, t): int(s == 0) + int(t == 0) - int(t == 1) for s, t in V}
I = dict(J)
for s in range(P):
    if s != 1:
        I[s, 0] -= 1
    if s != 0:
        I[s, 1] += 1
assert set(I.values()) == {0, 1} and sum(I.values()) == sum(J.values()) == P
assert not set(J.values()) <= {0, 1}
assert [sum(I[s, t] != J[s, t] for s in range(P)) for t in (0, 1)] == [28, 28]

result = {
    'p': P, 'directions': DIRECTIONS, 'stress_coefficients': SIGMA,
    'scope': 'First-carry section defect in the unramified anisotropic quadratic norm plane.',
    'theorem': 'A section with 89 failures on this triple cannot have one three-failure cycle.',
    'unchanged_lower_bounds': {'triple': 89, 'full_graph': 445},
    'single_cube_certificates': single_cube_certificates,
    'three_failure_cycle_geometries_checked': len(records),
    'histogram': [{'exceptional_direction': e, 'pair_kernel_dimensions': dims,
                   'coupled_rank': rank, 'count': count}
                  for (e, dims, rank), count in sorted(histogram.items())],
    'three_failure_cycle_certificates': records,
    'disjoint_cube_six_dimensional_hostile': hostile,
    'four_dimensional_separated_example': passing,
    'binary_normalization_boundary_hostile': {'weight': 29, 'changes_per_column': 28,
                                             'J_binary': False},
    'unresolved': '89 failures in two distinct double-failure cycles outside the correction criterion.'
}
Path(__file__).with_name('two-extra-exact-results.json').write_text(json.dumps(result, indent=2) + '\n')
print('42 single-cube controls: pair dimensions [2, 2, 2]; coupled homogeneous dimension 2.')
print('1134 three-failure-cycle geometries: support size 12; coupled homogeneous dimension 4 in every case.')
print('All 1134 direct 36-variable constraint systems have rank 32; explicit model rank 4.')
print('Largest-gap path normalization: maximum path length 19; all support bounds pass.')
print('Disjoint-cube hostile: 16 points; pair dimensions [6, 4, 6]; coupled rank 10; homogeneous dimension 6.')
print('The intermediate cube supplies exactly the two additional local-motion modes.')
print('Separated example: pair dimensions [4, 4, 4]; coupled rank 8; homogeneous dimension 4.')
print('Normalization boundary hostile: 28 changes per column can hide a nonbinary additive array.')
print('Scoped theorem: any 89-failure section must have two distinct double-failure cycles.')
print('Numerical bounds remain 89 on the triple and 445 on the full graph; attainment is unresolved.')
