#!/usr/bin/env python3
"""Independent hostile audit for the H cell-90 high-channel floor scout.

This deliberately does not import the primary C++ implementation.  It uses
three paths:

* literal Fraction intervals and a two-pointer intersection sweep;
* the O(min(p,q)) integer digital strip;
* brute floor-moment definitions for the Euclidean identities.

Scope is finite-exact audit plus an algebraic audit of the exceptional ray;
it does not turn the finite level ceiling into an infinite theorem.
"""

from fractions import Fraction as F
from hashlib import sha256
from math import gcd
from pathlib import Path
from random import Random

L = 168
CELL = 90
H = (1, 2, 3, 4, 6, 12)
EDGES = ((0, 1), (0, 2), (0, 3), (0, 4), (0, 5),
         (1, 2), (1, 3), (1, 4), (1, 5))
PRIMARY = Path(__file__).with_name("lrc_high_channel_floor_referee_thm3171.cpp")
EXPECTED_PRIMARY_SHA256 = "48b635746237b15c3a0ab2874b6d163c42b511019aaa8ae5aa3ba036ecaed1a6"


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def arcs(e, p):
    """Literal canonical reflected arcs, retaining exact endpoints."""
    z = L * p - e
    r = e * CELL % L
    rows = tuple((F(r + L * k - 12, z), F(r + L * k + 12, z))
                 for k in range(p))
    require(all(F(0) <= a < b <= 1 for a, b in rows),
            ("arc endpoint", e, p, rows))
    require(all(rows[k][1] < rows[k + 1][0] for k in range(p - 1)),
            ("arc disjointness", e, p, rows))
    return rows


def literal_mass(e, p, f, q):
    first, second = arcs(e, p), arcs(f, q)
    i = j = 0
    total = F(0)
    while i < len(first) and j < len(second):
        a, b = first[i]
        c, d = second[j]
        total += max(F(0), min(b, d) - max(a, c))
        if b < d:
            i += 1
        else:
            j += 1
    return total


def digital_numerator(e, p, f, q):
    """Independent linear scan of every possible nearest interval."""
    if p > q:
        return digital_numerator(f, q, e, p)
    z, w = L * p - e, L * q - f
    r, s = e * CELL % L, f * CELL % L
    determinant = r * w - s * z
    require(determinant % L == 0, ("nonintegral phase", e, p, f, q))
    base = determinant // L
    total = 0
    for k in range(p):
        x = base + k * w
        quotient = x // z
        contributing = 0
        for ell in (quotient, quotient + 1):
            if not 0 <= ell < q:
                continue
            delta = abs(x - ell * z)
            outer = max(0, 12 * (z + w) - L * delta)
            inner = max(0, 12 * abs(w - z) - L * delta)
            if outer - inner:
                contributing += 1
                total += outer - inner
        require(contributing <= 1, ("cap-two collision", e, p, f, q, k))
    return total, z * w


def brute_moments(n, m, a, b):
    values = tuple((a * k + b) // m for k in range(n))
    return (sum(values), sum(k * y for k, y in enumerate(values)),
            sum(y * y for y in values))


def recursive_moments(n, m, a, b):
    """Euclidean three-moment recursion, written independently in Python."""
    if n == 0:
        return 0, 0, 0
    s1 = n * (n - 1) // 2
    s2 = n * (n - 1) * (2 * n - 1) // 6
    qa, a0 = divmod(a, m)
    qb, b0 = divmod(b, m)
    base = (
        qa * s1 + qb * n,
        qa * s2 + qb * s1,
        qa * qa * s2 + 2 * qa * qb * s1 + qb * qb * n,
    )
    if a0 == 0:
        return base
    height = (a0 * (n - 1) + b0) // m
    if height == 0:
        return base
    u0, u1, u2 = recursive_moments(height, a0, m, m - b0 + a0 - 1)
    r0 = n * height - u0
    r1 = height * s1 - (u2 - u0) // 2
    r2 = n * height * height - 2 * u1 - u0
    return (base[0] + r0, base[1] + r1,
            base[2] + 2 * qa * r1 + 2 * qb * r0 + r2)


def exceptional_numerator(g):
    n = 2 * g // 3
    return 3276 * g * g - 1848 * g + n * (3024 * g - 1008 - 2268 * n)


def main():
    require(sha256(PRIMARY.read_bytes()).hexdigest() == EXPECTED_PRIMARY_SHA256,
            ("primary hash", sha256(PRIMARY.read_bytes()).hexdigest(),
             EXPECTED_PRIMARY_SHA256))

    # Literal arc geometry versus the integer digital strip.
    literal_controls = 0
    for i, j in EDGES:
        for p in range(6, 41):
            for q in range(p + 1, min(2 * p, 55) + 1):
                for e, f in ((H[i], H[j]), (H[j], H[i])):
                    numerator, denominator = digital_numerator(e, p, f, q)
                    require(literal_mass(e, p, f, q) == F(numerator, denominator),
                            ("literal/digital mismatch", e, p, f, q))
                    literal_controls += 1

    # Complete direct high-channel audit through level 300.
    pair_count = edge_checks = exceptional_checks = 0
    minima = [None] * 18
    for p in range(6, 301):
        for q in range(p + 1, min(2 * p, 300) + 1):
            common = gcd(p, q)
            P, Q = p // common, q // common
            if P + Q < 8:
                continue
            pair_count += 1
            for edge, (i, j) in enumerate(EDGES):
                for orientation, (e, f) in enumerate(((H[i], H[j]), (H[j], H[i]))):
                    numerator, denominator = digital_numerator(e, p, f, q)
                    row = (F(numerator, denominator), p, q)
                    index = 2 * edge + orientation
                    if minima[index] is None or row < minima[index]:
                        minima[index] = row
                    edge_checks += 1
                    exception = edge == 3 and (P, Q) == (3, 5)
                    if exception:
                        exceptional_checks += 1
                        require(numerator * 280393 >= 2030 * denominator,
                                ("exceptional floor", e, p, f, q))
                    else:
                        require(numerator * 105 >= denominator,
                                ("ordinary floor", e, p, f, q))

    # Hostile exact checks of the Euclidean recursion itself.
    rng = Random(20260802)
    moment_controls = 0
    for _ in range(50000):
        n = rng.randrange(0, 80)
        m = rng.randrange(1, 80)
        a = rng.randrange(0, 240)
        b = rng.randrange(0, 240)
        require(recursive_moments(n, m, a, b) == brute_moments(n, m, a, b),
                ("moment recursion", n, m, a, b))
        moment_controls += 1

    # Algebraic exceptional ray.  The three residue classes give the same
    # floor formula and every consecutive cleared difference is positive.
    constants = (0, -336, 84)
    ray_controls = 0
    for g in range(1, 1001):
        polynomial = 4284 * g * g - 2520 * g + constants[g % 3]
        require(exceptional_numerator(g) == polynomial,
                ("exceptional polynomial", g))
        direct, denominator = digital_numerator(1, 5 * g, 6, 3 * g)
        require(direct == polynomial, ("exceptional direct", g))
        require(denominator == (840 * g - 1) * (504 * g - 6),
                ("exceptional denominator", g))
        ray_controls += 1
    require(F(exceptional_numerator(2), (1680 - 1) * (1008 - 6))
            == F(2030, 280393), "sharp exceptional value")
    require(all(c > 0 for c in (
        # Positive polynomial coefficients after clearing I_(g+1)-I_g.
        1787436, 2073474, 17,
        1211238, 1314819, 139285,
        999558, 964791, -34808 + 999558 + 964791,
    )), "exceptional monotonicity coefficients")

    semantic = repr((literal_controls, pair_count, edge_checks, exceptional_checks,
                     tuple(minima), moment_controls, ray_controls)).encode()
    print("LRC H HIGH-CHANNEL FLOOR -- INDEPENDENT HOSTILE AUDIT")
    print(f"literal_arc_controls={literal_controls}")
    print(f"finite_universe=6<=p<q<=300;high_level_pairs={pair_count};"
          f"edge_orientation_checks={edge_checks};exception_checks={exceptional_checks}")
    print(f"moment_recursion_hostile_controls={moment_controls}")
    print(f"exceptional_ray_formula_controls={ray_controls};sharp_g2=2030/280393")
    print(f"semantic_sha256={sha256(semantic).hexdigest()}")


if __name__ == "__main__":
    main()
