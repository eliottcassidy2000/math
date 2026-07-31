#!/usr/bin/env python3
"""Exact referee for THM-2647.

All rational operations use ``Fraction`` and every logical check uses
``require`` so optimized Python executes the same certificate.
"""

from fractions import Fraction
from itertools import combinations
from math import gcd


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def convolution(left, right):
    n = len(left)
    require(len(right) == n, "convolution size mismatch")
    return tuple(
        sum(left[a] * right[(x - a) % n] for a in range(n))
        for x in range(n)
    )


def delta(n, position):
    return tuple(Fraction(int(x == position % n)) for x in range(n))


def mask(n, support):
    support = set(support)
    return tuple(int(x in support) for x in range(n))


def pair_convolution_profile(left_pair, right_pair, n):
    out = [0] * n
    for left in left_pair:
        for right in right_pair:
            out[(left + right) % n] += 1
    return tuple(out)


def translate(function, shift):
    n = len(function)
    return tuple(function[(x - shift) % n] for x in range(n))


def additive_order(d, n):
    require(d % n != 0, "nonzero step required")
    return n // gcd(d, n)


def two_point_inverse(n, d):
    q = additive_order(d, n)
    require(q % 2 == 1, "odd step order required")
    out = [Fraction(0) for _ in range(n)]
    for j in range(q):
        out[(j * d) % n] += Fraction(1 if j % 2 == 0 else -1, 2)
    return tuple(out)


def convolution_matrix(kernel):
    n = len(kernel)
    return [[kernel[(row - column) % n] for column in range(n)] for row in range(n)]


def bareiss_determinant(matrix):
    a = [[int(value) for value in row] for row in matrix]
    n = len(a)
    if n == 0:
        return 1
    sign = 1
    previous = 1
    for k in range(n - 1):
        pivot_row = next((row for row in range(k, n) if a[row][k] != 0), None)
        if pivot_row is None:
            return 0
        if pivot_row != k:
            a[k], a[pivot_row] = a[pivot_row], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = a[i][j] * pivot - a[i][k] * a[k][j]
                require(numerator % previous == 0, "Bareiss exact division")
                a[i][j] = numerator // previous
        for i in range(k + 1, n):
            a[i][k] = 0
        previous = pivot
    return sign * a[n - 1][n - 1]


def centered_energy(profile):
    n = len(profile)
    mean = Fraction(sum(profile), n)
    return sum((value - mean) ** 2 for value in profile) / n


def main():
    odd_inverse_checks = 0
    determinant_checks = 0
    for n in (3, 5, 7, 9, 11, 13):
        for d in range(1, n):
            q = additive_order(d, n)
            inverse = two_point_inverse(n, d)
            kernel = tuple(delta(n, 0)[x] + delta(n, d)[x] for x in range(n))
            require(convolution(kernel, inverse) == delta(n, 0),
                    "odd-cycle inverse failure")
            require(sum(abs(value) for value in inverse) == Fraction(q, 2),
                    "odd-cycle l1 invoice")
            require(sum(inverse) == Fraction(1, 2), "inverse augmentation")
            determinant = bareiss_determinant(convolution_matrix(kernel))
            require(abs(determinant) == 2 ** (n // q), "odd-cycle determinant")
            odd_inverse_checks += 1
            determinant_checks += 1

    even_hostiles = 0
    for n in (2, 4, 6, 8, 10, 12):
        kernel = tuple(delta(n, 0)[x] + delta(n, 1)[x] for x in range(n))
        alternating = tuple(Fraction(1 if x % 2 == 0 else -1) for x in range(n))
        require(convolution(kernel, alternating) == (Fraction(0),) * n,
                "even-cycle alternating hostile")
        require(bareiss_determinant(convolution_matrix(kernel)) == 0,
                "even-cycle singular determinant")
        even_hostiles += 1

    p = 13
    cayley_checks = 0
    for d in range(1, p):
        inverse = two_point_inverse(p, d)
        identity = delta(p, 0)
        h = tuple(identity[x] - 2 * inverse[x] for x in range(p))
        kernel = tuple(identity[x] + delta(p, d)[x] for x in range(p))
        rhs = tuple(delta(p, d)[x] - identity[x] for x in range(p))
        require(convolution(kernel, h) == rhs, "Cayley telescoping identity")
        require(tuple((identity[x] - h[x]) / 2 for x in range(p)) == inverse,
                "Cayley inverse formula")
        require(sum(value > 0 for value in inverse) == 7 and
                sum(value < 0 for value in inverse) == 6,
                "C13 signed inverse census")
        require(sum(abs(value) for value in inverse) == Fraction(13, 2),
                "C13 thirteen-halves tax")
        require(abs(bareiss_determinant(convolution_matrix(kernel))) == 2,
                "C13 determinant two")
        cayley_checks += 1

    pairs = tuple(combinations(range(p), 2))
    require(len(pairs) == 78, "C13 two-set count")

    # For each fixed left endpoint, exhaust all 78 right endpoints once and
    # verify that their multiplicity profiles are pairwise distinct.  This is
    # the finite-bank form of uniqueness against every candidate B.
    fixed_endpoint_profiles = 0
    for a_pair in pairs:
        profiles = {
            pair_convolution_profile(a_pair, b_pair, p)
            for b_pair in pairs
        }
        require(len(profiles) == 78, "fixed endpoint convolution is not injective")
        fixed_endpoint_profiles += len(profiles)

    recovery_checks = 0
    gauge_checks = 0
    one_point_anchor_checks = 0

    for a_pair in pairs:
        a_mask = mask(p, a_pair)
        a_complement = tuple(1 - value for value in a_mask)
        for b_pair in pairs:
            b_mask = mask(p, b_pair)
            b_complement = tuple(1 - value for value in b_mask)
            multiplicity = pair_convolution_profile(a_pair, b_pair, p)
            dense_profile = convolution(a_complement, b_complement)
            recovered_multiplicity = tuple(value - 9 for value in dense_profile)
            require(recovered_multiplicity == multiplicity, "r=9+m recovery")

            for a0, a1 in (a_pair, tuple(reversed(a_pair))):
                d = (a1 - a0) % p
                inverse = two_point_inverse(p, d)
                twice_inverse = tuple(int(2 * value) for value in inverse)
                recovered_twice = convolution(twice_inverse,
                                                translate(multiplicity, -a0))
                require(recovered_twice == tuple(2 * value for value in b_mask),
                        "endpoint deconvolution failure")
                require(all(value in (0, 2) for value in recovered_twice),
                        "deconvolution did not return a Boolean endpoint")
                recovery_checks += 1

            orbit = set()
            one_point_representatives = set()
            for shift in range(p):
                moved_a = tuple(sorted((a + shift) % p for a in a_pair))
                moved_b = tuple(sorted((b - shift) % p for b in b_pair))
                require(pair_convolution_profile(moved_a, moved_b, p) == multiplicity,
                        "common-origin gauge changed multiplicity")
                state = (moved_a, moved_b)
                orbit.add(state)
                if 0 in moved_a:
                    one_point_representatives.add(state)
                gauge_checks += 1
            require(len(orbit) == p, "common-origin gauge is not free")
            require(len(one_point_representatives) == 2,
                    "one-point anchor did not leave exactly two orbit representatives")
            one_point_anchor_checks += 1

    require(recovery_checks == 2 * 78 * 78, "orientation recovery count")
    require(fixed_endpoint_profiles == 78 * 78, "candidate uniqueness count")
    require(gauge_checks == 78 * 78 * 13, "gauge check count")
    require(one_point_anchor_checks == 78 * 78, "one-point anchor count")

    hostile_a = mask(p, (0, 1))
    hostile_b2 = mask(p, (0, 2))
    hostile_b3 = mask(p, (0, 3))
    hostile_r2 = convolution(tuple(1 - x for x in hostile_a),
                             tuple(1 - x for x in hostile_b2))
    hostile_r3 = convolution(tuple(1 - x for x in hostile_a),
                             tuple(1 - x for x in hostile_b3))
    require(hostile_r2 != hostile_r3, "energy hostile profiles coincided")
    require(centered_energy(hostile_r2) == Fraction(36, 169) and
            centered_energy(hostile_r3) == Fraction(36, 169),
            "energy-only hostile")

    print("THM2647 endpoint-anchored two-point deconvolution exact referee")
    print(f"odd_inverse_checks={odd_inverse_checks} determinant_checks={determinant_checks}")
    print(f"even_alternating_hostiles={even_hostiles}")
    print(f"C13_cayley_checks={cayley_checks} determinant=2 l1_inverse=13/2 signs=7plus6minus")
    print(f"ordered_endpoint_pairs={len(pairs) ** 2} oriented_recoveries={recovery_checks}")
    print(f"fixed_A_candidate_profiles={fixed_endpoint_profiles} solutions_per_fixed_A=1")
    print(f"common_origin_gauge_checks={gauge_checks} orbit_size=13")
    print(f"one_point_anchor_checks={one_point_anchor_checks} representatives_per_orbit=2")
    print("energy_only_hostile=distinct_profiles_same_36/169")
    print("PASS")


if __name__ == "__main__":
    main()
