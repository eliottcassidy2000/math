#!/usr/bin/env python3
"""Dependency-free exact referee for THM-3048."""

from fractions import Fraction
import json


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def jdump(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def eye(n):
    return [[int(i == j) for j in range(n)] for i in range(n)]


def zero(n):
    return [[0 for _ in range(n)] for _ in range(n)]


def add(a, b):
    return [[a[i][j] + b[i][j] for j in range(len(a[0]))] for i in range(len(a))]


def matmul(a, b):
    return [
        [sum(a[i][k] * b[k][j] for k in range(len(b))) for j in range(len(b[0]))]
        for i in range(len(a))
    ]


def translation(n, shift):
    shift %= n
    out = zero(n)
    for j in range(n):
        out[(j + shift) % n][j] = 1
    return out


def trace(a):
    return sum(a[i][i] for i in range(len(a)))


def total(a):
    return sum(sum(row) for row in a)


def det_i_plus_t(a):
    """Coefficients of det(I+tA), by Newton identities."""
    n = len(a)
    powers = [eye(n)]
    traces = [None]
    for _ in range(1, n + 1):
        powers.append(matmul(powers[-1], a))
        traces.append(trace(powers[-1]))
    e = [Fraction(1)]
    for k in range(1, n + 1):
        value = sum(((-1) ** (j - 1)) * e[k - j] * traces[j] for j in range(1, k + 1))
        e.append(value / k)
    require(all(x.denominator == 1 for x in e), "characteristic coefficient lost integrality")
    return [int(x) for x in e]


def poly_mul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def poly_pow(a, exponent):
    out = [1]
    base = list(a)
    while exponent:
        if exponent & 1:
            out = poly_mul(out, base)
        base = poly_mul(base, base)
        exponent //= 2
    return out


def orbit_sum(seed, u, v):
    n = len(seed)
    out = zero(n)
    for j in range(n):
        out = add(out, matmul(matmul(translation(n, j * u), seed), translation(n, -j * v)))
    return out


def seed_matrix(n, seed):
    return [
        [((3 * i + 5 * j + seed) % 7 if (i + 2 * j + seed) % 3 == 0 else 0) for j in range(n)]
        for i in range(n)
    ]


print("THM-3048 PRIME-HOLONOMY PARALLEL CORRESPONDENCE ARRIVAL")

parallel_cells = 0
parallel_digest = []
for n, pairs in (
    (3, ((0, 1), (1, 0))),
    (5, ((0, 1), (2, 4))),
    (13, ((0, 1), (2, 7), (5, 12))),
):
    for u, v in pairs:
        g = (v - u) % n
        require(g != 0, "test holonomy vanished")
        for seed_id in (1, 2):
            seed = seed_matrix(n, seed_id)
            cmat = orbit_sum(seed, u, v)
            require(cmat == matmul(matmul(translation(n, u), cmat), translation(n, -v)), "parallelism failed")
            determinants = []
            ghosts = []
            for c in range(n):
                pointed = matmul(cmat, translation(n, c))
                determinants.append(det_i_plus_t(pointed))
                ghosts.append(trace(pointed))
                conjugate = matmul(matmul(translation(n, u), pointed), translation(n, -u))
                require(matmul(cmat, translation(n, c + g)) == conjugate, "pointing covariance failed")
            require(all(d == determinants[0] for d in determinants), "prime holonomy did not flatten determinants")
            require(all(x == total(cmat) // n for x in ghosts), "first ghost did not equal total/n")
            require(total(cmat) % n == 0 and total(cmat) > 0, "invalid positive root mass")
            norm = [1]
            for dpoly in determinants:
                norm = poly_mul(norm, dpoly)
            require(norm == poly_pow(determinants[0], n), "pointing norm is not a prime power")
            parallel_cells += 1
            if len(parallel_digest) < 8:
                parallel_digest.append({"g": g, "ghost": ghosts[0], "n": n, "seed": seed_id})
print(jdump({"parallel_prime_cells": parallel_cells, "digest": parallel_digest}))

# The unpointed norm always has first ghost equal to total root mass.
norm_cells = 0
for n in (3, 5, 13):
    for seed_id in (1, 2, 3):
        cmat = seed_matrix(n, seed_id)
        determinants = [det_i_plus_t(matmul(cmat, translation(n, c))) for c in range(n)]
        require(sum(poly[1] for poly in determinants) == total(cmat), "norm first ghost lost total mass")
        norm_cells += 1
print(jdump({"unpointed_norm_first_ghost_cells": norm_cells}))

# Pure shifts show that the norm forgets semantic displacement.
pure_shift_cells = 0
for n in (3, 5, 13):
    expected = poly_mul(poly_pow([1, 1], n), poly_pow([1] + [0] * (n - 1) + [1], n - 1))
    for d in range(n):
        cmat = translation(n, d)
        norm = [1]
        for c in range(n):
            norm = poly_mul(norm, det_i_plus_t(matmul(cmat, translation(n, c))))
        require(norm == expected, "pure-shift norm depends on displacement")
        pure_shift_cells += 1
print(jdump({"pure_shift_norm_cells": pure_shift_cells, "p13_formula":"(1+t)^13(1+t^13)^12"}))

# Sharp failures: zero holonomy, nonparallel transport, and composite nontransitivity.
p = 13
zero_holonomy_c = translation(p, 1)
require(zero_holonomy_c == matmul(matmul(translation(p, 1), zero_holonomy_c), translation(p, -1)), "zero-holonomy control not parallel")
require(trace(zero_holonomy_c) == 0 and total(zero_holonomy_c) == p, "zero-holonomy hostile changed")

nonparallel_c = translation(p, 1)
require(nonparallel_c != matmul(nonparallel_c, translation(p, -1)), "nonparallel hostile became parallel")
require(trace(nonparallel_c) == 0 and total(nonparallel_c) == p, "nonparallel hostile changed")

n = 6
composite_c = [[int(((b - h) % n) in {1, 3, 5}) for b in range(n)] for h in range(n)]
require(composite_c == matmul(composite_c, translation(n, -2)), "composite hostile not parallel")
composite_ghosts = [trace(matmul(composite_c, translation(n, c))) for c in range(n)]
require(composite_ghosts == [0, 6, 0, 6, 0, 6], "composite orbit hostile changed")

holonomy_orbits = []
for a in range(1, 13):
    g = (7 * a) % 13
    orbit = sorted({(j * g) % 13 for j in range(13)})
    require(orbit == list(range(13)), "THM-2591 nonzero holonomy lost transitivity")
    holonomy_orbits.append(g)
print(jdump({"boundaries":{"cemetery_only_root_mass":0,"composite_C6_g2":composite_ghosts,"nonparallel_arrival":0,"zero_holonomy_arrival":0},"thm2591_holonomy_permutation":holonomy_orbits}))

print("all_exact_checks=PASS")
