#!/usr/bin/env python3
"""Exact companion for THM-2413.

Checks operation-shadow and diagonal-deletion laws, prime valuation
coordinates, the twin-center sum/product weld, the prime-index slope-two
drift identity, a finite exact A373813 line-cover prefix, and the reciprocal
twin-center identity.  Standard-library integer/Fraction arithmetic only.
"""

from fractions import Fraction
from functools import lru_cache
from math import isqrt


def is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    p = 3
    while p * p <= n:
        if n % p == 0:
            return False
        p += 2
    return True


def primes_first(count):
    out = []
    n = 2
    while len(out) < count:
        if is_prime(n):
            out.append(n)
        n += 1
    return out


def factorization(n):
    out = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0) + 1
            n //= p
        p = 3 if p == 2 else p + 2
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def divisors(n):
    out = [1]
    for p, e in factorization(n).items():
        out = [d * p**j for d in out for j in range(e + 1)]
    return sorted(out)


def tau(n):
    out = 1
    for e in factorization(n).values():
        out *= e + 1
    return out


def collinear(points, i, j, k):
    xi, yi = points[i]
    xj, yj = points[j]
    xk, yk = points[k]
    return (xj - xi) * (yk - yi) == (xk - xi) * (yj - yi)


def affine_line_masks(points):
    n = len(points)
    masks = {1 << i for i in range(n)}
    for i in range(n):
        for j in range(i + 1, n):
            mask = 0
            for k in range(n):
                if collinear(points, i, j, k):
                    mask |= 1 << k
            masks.add(mask)

    # A proper subset of another available collinear mask is never useful in
    # an unweighted minimum cover.
    ordered = sorted(masks, key=int.bit_count, reverse=True)
    maximal = []
    for mask in ordered:
        if not any(mask | other == other for other in maximal):
            maximal.append(mask)
    return tuple(maximal)


def minimum_line_cover(points):
    n = len(points)
    full = (1 << n) - 1
    lines = affine_line_masks(points)
    containing = [[] for _ in range(n)]
    for mask in lines:
        for i in range(n):
            if mask >> i & 1:
                containing[i].append(mask)

    @lru_cache(None)
    def solve(covered):
        if covered == full:
            return 0
        uncovered = full ^ covered
        i = (uncovered & -uncovered).bit_length() - 1
        return 1 + min(solve(covered | mask) for mask in containing[i])

    return solve(0)


def main():
    limit = 80
    additive_edges = [
        (x, z) for x in range(1, limit + 1) for z in range(x + 1, limit + 1)
    ]
    multiplicative_edges_with_unit = [
        (x, z)
        for x in range(1, limit + 1)
        for z in range(x + 1, limit + 1)
        if z % x == 0
    ]
    multiplicative_edges_nonunit = [
        (x, z) for x, z in multiplicative_edges_with_unit if x >= 2
    ]
    assert len(additive_edges) == limit * (limit - 1) // 2 == 3160
    assert len(multiplicative_edges_with_unit) == 288
    assert len(multiplicative_edges_nonunit) == 209

    checked_defects = 0
    for x in range(2, limit + 1):
        additive_defects = [(c, x - c) for c in range(1, x)]
        multiplicative_defects = [
            (c, x // c)
            for c in range(2, x)
            if x % c == 0 and x // c >= 2
        ]
        assert len(additive_defects) == x - 1
        assert len(multiplicative_defects) == tau(x) - 2
        assert (len(multiplicative_defects) == 0) == is_prime(x)

        # Each pair is a genuine two-step path whose composite direct edge is
        # precisely the deleted equal-parent edge.
        for c, d in additive_defects:
            assert c != x  # first equal-parent-deleted step exists
            assert d != x + c  # second such step exists
            assert (x + c) + d == 2 * x
            assert c + d == x  # direct edge is deleted
        for c, d in multiplicative_defects:
            assert c != x
            assert d != x * c
            assert (x * c) * d == x * x
            assert c * d == x  # direct edge is deleted
        checked_defects += 1

    checked_valuations = 0
    for a in range(1, 51):
        va = factorization(a)
        for b in range(1, 51):
            vb = factorization(b)
            vab = factorization(a * b)
            primes = set(va) | set(vb) | set(vab)
            assert all(
                vab.get(p, 0) == va.get(p, 0) + vb.get(p, 0)
                for p in primes
            )
            checked_valuations += 1

    prime_list = primes_first(40)
    twin_centers = [
        prime_list[k] + 1
        for k in range(len(prime_list) - 1)
        if prime_list[k + 1] - prime_list[k] == 2
    ]
    assert twin_centers[:12] == [4, 6, 12, 18, 30, 42, 60, 72, 102, 108, 138, 150]

    checked_twin_equivalences = 0
    checked_drift = 0
    reciprocal_midpoints = Fraction(0)
    reciprocal_pairs_half = Fraction(0)
    reciprocal_correction = Fraction(0)
    drift = [p - 2 * (k + 1) for k, p in enumerate(prime_list)]
    for k in range(len(prime_list) - 1):
        p = prime_list[k]
        q = prime_list[k + 1]
        plateau = drift[k + 1] == drift[k]
        twin = q - p == 2
        assert plateau == twin
        checked_drift += 1
        if not twin:
            continue

        m = p + 1
        n = m * m - 1
        assert q == m + 1
        assert n == p * q
        assert divisors(n) == [1, p, q, n]
        assert tau(n) == 4
        s = p + q
        product = p * q
        assert s == 2 * m
        assert product == m * m - 1
        assert s * s - 4 * product == 4
        assert product + 1 == m * m
        assert m == 4 or m % 6 == 0

        reciprocal_midpoints += Fraction(1, m)
        reciprocal_pairs_half += Fraction(1, 2) * (
            Fraction(1, p) + Fraction(1, q)
        )
        reciprocal_correction += Fraction(1, m * (m * m - 1))
        checked_twin_equivalences += 1

    assert (
        reciprocal_midpoints
        == reciprocal_pairs_half - reciprocal_correction
    )

    # Converse finite hostile scan for the discriminant/factor-diamond weld.
    for m in range(4, 500):
        twin = is_prime(m - 1) and is_prime(m + 1)
        n = m * m - 1
        diamond_gap_two = (
            len(divisors(n)) == 4
            and divisors(n)[2] - divisors(n)[1] == 2
            and divisors(n)[1] + divisors(n)[2] == 2 * m
        )
        assert twin == diamond_gap_two

    # Exact A373813 prefix from the prime-index point hypergraph.
    line_prefix = []
    for n in range(1, 17):
        points = [(k + 1, prime_list[k]) for k in range(n)]
        line_prefix.append(minimum_line_cover(points))
    assert line_prefix == [1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4]
    assert all(
        line_prefix[i] <= line_prefix[i + 1] <= line_prefix[i] + 1
        for i in range(len(line_prefix) - 1)
    )
    awkward_indices = [
        i + 1
        for i, value in enumerate(line_prefix)
        if value > (line_prefix[i - 1] if i else 0)
    ]
    assert len(awkward_indices) == line_prefix[-1]

    print("theorem=THM-2413")
    print("arithmetic=integer/Fraction-only")
    print(f"additive-shadow-edges-N80={len(additive_edges)}")
    print(
        "multiplicative-shadow-edges-N80="
        f"with-unit:{len(multiplicative_edges_with_unit)},"
        f"nonunit-source:{len(multiplicative_edges_nonunit)}"
    )
    print(f"diagonal-defect-vertices-checked={checked_defects}")
    print(f"valuation-products-checked={checked_valuations}")
    print(f"twin-centers-checked={checked_twin_equivalences}")
    print("A014574-first=" + ",".join(map(str, twin_centers[:12])))
    print(f"prime-drift-adjacencies-checked={checked_drift}")
    print("twin-drift-law=d_(k+1)=d_k iff prime-gap=2")
    print(
        "reciprocal-midpoint-identity="
        "sum(1/m)=1/2 sum(1/(m-1)+1/(m+1))-sum(1/(m(m^2-1)))"
    )
    print("A373813-prefix=" + ",".join(map(str, line_prefix)))
    print("awkward-indices-prefix=" + ",".join(map(str, awkward_indices)))
    print("A373813-guard=no equivalence with its plateaux or optimal-line slopes")
    print("status=PASS")


if __name__ == "__main__":
    main()
