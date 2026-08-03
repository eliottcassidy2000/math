#!/usr/bin/env python3
"""Exact hostile audit for the THM-3211 boundary-cocycle reframe.

The paper-level cocycle and aggregate bounds are proved in THM-3211.  This
companion checks their finite controls, hostiles, and declared census.
"""

from fractions import Fraction as F
from math import gcd

L = 168
THETA = F(1, 14)
H = (1, 2, 3, 4, 6, 12)
EDGES = ((0, 1), (0, 2), (0, 3), (0, 4), (0, 5),
         (1, 2), (1, 3), (1, 4), (1, 5))
LANES = tuple((H[a], H[b]) for a, b in EDGES for (a, b) in ((a, b), (b, a)))


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def frac(x):
    return x - (x.numerator // x.denominator)


def b2bar(x):
    t = frac(x)
    return t * t - t + F(1, 6)


def b3bar(x):
    t = frac(x)
    return t ** 3 - F(3, 2) * t ** 2 + F(1, 2) * t


def comb(P, phase):
    """Disjoint clipped representatives of chi(P*x-phase) on [0,1]."""
    phase = frac(phase)
    ans = []
    for k in range(-2, P + 3):
        lo = (F(k) + phase - THETA) / P
        hi = (F(k) + phase + THETA) / P
        lo, hi = max(F(0), lo), min(F(1), hi)
        if lo < hi:
            ans.append((lo, hi))
    ans.sort()
    return ans


def overlap_components(P, alpha, Q, beta):
    left, right = comb(P, alpha), comb(Q, beta)
    i = j = 0
    ans = []
    while i < len(left) and j < len(right):
        lo = max(left[i][0], right[j][0])
        hi = min(left[i][1], right[j][1])
        if lo < hi:
            ans.append((lo, hi))
        if left[i][1] < right[j][1]:
            i += 1
        elif right[j][1] < left[i][1]:
            j += 1
        else:
            i += 1
            j += 1
    return ans


def barycenter_direct(P, alpha, Q, beta):
    return sum(((hi - F(1, 2)) ** 2 - (lo - F(1, 2)) ** 2) / 2
               for lo, hi in overlap_components(P, alpha, Q, beta))


def barycenter_b2(P, alpha, Q, beta):
    # sigma=+1 at a left endpoint and -1 at a right endpoint.
    return -sum(b2bar(lo) - b2bar(hi)
                for lo, hi in overlap_components(P, alpha, Q, beta)) / 2


def ccell(P, Q, e, f, j):
    a0, b0 = F(j * e, L), F(j * f, L)
    a1, b1 = F((j + 1) * e, L), F((j + 1) * f, L)
    return barycenter_direct(P, a1, Q, b1) - barycenter_direct(P, a0, Q, b0)


def limit_cell(P, Q, e, f, j):
    C = Q * e - P * f
    require(C != 0, ("rank-one limit", P, Q, e, f))
    R, S = j * e % L, j * f % L
    D = Q * R - P * S
    a, b = F(P, 14), F(Q, 14)
    u, v = F(D + C, L), F(-D, L)
    psi = (b3bar(u + a - b) + b3bar(u - a + b)
           + b3bar(v + a - b) + b3bar(v - a + b)
           - b3bar(u + a + b) - b3bar(u - a - b)
           - b3bar(v + a + b) - b3bar(v - a - b))
    return F(1, 49) + F(28, P * Q * C) * psi


def located_comb(p, e, j):
    z, r = L * p - e, j * e % L
    ans = []
    for k in range(-2, p + 3):
        lo, hi = F(r + L * k - 12, z), F(r + L * k + 12, z)
        lo, hi = max(F(0), lo), min(F(1), hi)
        if lo < hi:
            ans.append((lo, hi))
    ans.sort()
    return ans


def interval_overlap(left, right, moment=False):
    i = j = 0
    total = F(0)
    while i < len(left) and j < len(right):
        lo, hi = max(left[i][0], right[j][0]), min(left[i][1], right[j][1])
        if lo < hi:
            total += (((hi - F(1, 2)) ** 2 - (lo - F(1, 2)) ** 2) / 2
                      if moment else hi - lo)
        if left[i][1] < right[j][1]:
            i += 1
        elif right[j][1] < left[i][1]:
            j += 1
        else:
            i += 1
            j += 1
    return total


def cell_mass(P, Q, e, f, j, g):
    return interval_overlap(located_comb(g * P, e, j), located_comb(g * Q, f, j))


def orbit_order(e, f):
    return L // gcd(L, gcd(e, f))


# Independent B2 endpoint identity, oddness, exact path additivity, and the
# finite cyclic-cell law on a nontrivial but broad exact universe.
states = 0
cell_signs = {"positive": 0, "negative": 0, "zero": 0}
all_zero_orbits = 0
zero_records = []
for Q in range(2, 26):
    for P in range((Q + 1) // 2, Q):
        if gcd(P, Q) != 1 or P + Q < 8:
            continue
        for e, f in LANES:
            T = orbit_order(e, f)
            potentials = []
            bulk_sum = F(0)
            for j in range(T + 1):
                a, b = F(j * e, L), F(j * f, L)
                components = overlap_components(P, a, Q, b)
                direct = sum(((hi - F(1, 2)) ** 2 - (lo - F(1, 2)) ** 2) / 2
                             for lo, hi in components)
                endpoint = -sum(b2bar(lo) - b2bar(hi) for lo, hi in components) / 2
                require(direct == endpoint, ("B2 endpoint mismatch", P, Q, e, f, j))
                potentials.append(direct)
                if j < T:
                    bulk_sum += limit_cell(P, Q, e, f, j)
                states += 1
            require(bulk_sum == F(T, 49), ("universal bulk mean", P, Q, e, f, bulk_sum))
            require(all(potentials[T - j] == -potentials[j] for j in range(T + 1)),
                    ("oddness mismatch", P, Q, e, f))
            corrections = tuple(potentials[j + 1] - potentials[j] for j in range(T))
            require(sum(corrections) == 0, ("cycle sum", P, Q, e, f))
            require(potentials[0] == potentials[-1] == 0, ("closed potential", P, Q, e, f))
            require(all(corrections[T - 1 - j] == corrections[j] for j in range(T)),
                    ("reflection palindrome", P, Q, e, f))
            # Three-segment hostile additivity, including wrapped phases.
            cuts = (0, T // 3, (2 * T) // 3, T)
            pieces = tuple(potentials[cuts[k + 1]] - potentials[cuts[k]] for k in range(3))
            require(sum(pieces) == 0, ("path additivity", P, Q, e, f, pieces))
            for j, value in enumerate(corrections):
                sign = "positive" if value > 0 else "negative" if value < 0 else "zero"
                cell_signs[sign] += 1
                if value == 0:
                    zero_records.append((P, Q, e, f, j))
            if not any(corrections):
                all_zero_orbits += 1


# Explicit THM-3211 values at cell 90.
controls = {
    (3, 5, 1, 2): F(71, 264600),
    (3, 5, 3, 1): F(23, 12600),
    (3, 5, 2, 3): F(1, 5400),
    (3, 5, 6, 1): -F(8213, 1411200),
}
for key, expected in controls.items():
    P, Q, e, f = key
    require(ccell(P, Q, e, f, 90) == expected, ("cell90 control", key))


# Direct finite-cell-to-closed-geodesic identity and its explicit O(g^-2)
# complete-orbit bound on the four canonical rays.
aggregate_controls = []
for P, Q, e, f in controls:
    d = gcd(L, gcd(e, f))
    T, a, b = L // d, e // d, f // d
    K = abs(Q * e - P * f) // d
    for g in (1, 2):
        n, m = g * P * T - a, g * Q * T - b
        h = gcd(n, m)
        require(K % h == 0, ("gcd determinant divisor", P, Q, e, f, g, h, K))
        direct = sum(cell_mass(P, Q, e, f, j, g) for j in range(T))
        closed = T * interval_overlap(comb(n, F(0)), comb(m, F(0)))
        require(direct == closed, ("closed geodesic identity", P, Q, e, f, g))
        error = direct - F(T, 49)
        bound = F(T * h * h, 3 * n * m)
        require(abs(error) <= bound, ("complete orbit bound", P, Q, e, f, g, error, bound))
        aggregate_controls.append((P, Q, e, f, g, direct, error, bound, h))

# The complete-orbit bound is two-sided, not a hidden lower bound.
negative_aggregate_key = (3, 5, 1, 12, 2)
P, Q, e, f, g = negative_aggregate_key
d = gcd(L, gcd(e, f))
T, a, b = L // d, e // d, f // d
n, m = g * P * T - a, g * Q * T - b
negative_aggregate = T * interval_overlap(comb(n, F(0)), comb(m, F(0))) - F(T, 49)
require(negative_aggregate == -F(10, 979811), ("negative aggregate hostile", negative_aggregate))


# A zero first correction does not kill the periodic second-order term or
# settle a finite head.
zero_hostile = (5, 6, 1, 12, 56)
P, Q, e, f, j = zero_hostile
require(ccell(P, Q, e, f, j) == 0, ("zero hostile correction", zero_hostile))
zero_limit = limit_cell(P, Q, e, f, j)
zero_heads = tuple((g, cell_mass(P, Q, e, f, j, g) - zero_limit) for g in (1, 2, 6, 9))
require(zero_limit == F(4, 189), ("zero hostile limit", zero_limit))
require(zero_heads == ((1, -F(340, 13161393)),
                       (2, F(11483, 52994277)),
                       (6, -F(697, 479042613)),
                       (9, -F(4, 1078631505))), ("zero hostile heads", zero_heads))


# Sign/zero anatomy on the full cell orbit for the four canonical rays.
print("LRC BOUNDARY COCYCLE EXACT AUDIT")
print(f"exact_phase_states={states};Q<=25;lanes={len(LANES)}")
print(f"cell_signs={cell_signs};all_zero_orbits={all_zero_orbits};zero_records={zero_records}")
print(f"aggregate_controls={aggregate_controls}")
print(f"negative_aggregate_hostile={negative_aggregate_key};error={negative_aggregate}")
print(f"zero_correction_hostile={zero_hostile};limit={zero_limit};head_differences={zero_heads}")
for P, Q, e, f in controls:
    T = orbit_order(e, f)
    vals = tuple(ccell(P, Q, e, f, j) for j in range(T))
    pos = sum(v > 0 for v in vals)
    neg = sum(v < 0 for v in vals)
    zero = sum(v == 0 for v in vals)
    print(f"ray={(P,Q,e,f)};T={T};signs=+{pos}/-{neg}/0:{zero};cell90={vals[90 % T]};maxabs={max(map(abs, vals))}")
print("B2_ENDPOINT_IDENTITY=PASS")
print("EXACT_COBoundary_CYCLE_SUM=PASS")
print("REFLECTION_PALINDROME=PASS")
print("UNIVERSAL_COMPLETE_ORBIT_BULK=T/49=PASS")
print("FINITE_CLOSED_GEODESIC_BOUND=PASS")
print("FAILED_CHECKS=NONE")
