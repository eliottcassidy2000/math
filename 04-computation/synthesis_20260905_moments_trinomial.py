"""Exact audit of trinomial first-return classification and carry profile.

Reproduce: python 04-computation/synthesis_20260905_moments_trinomial.py
Universe: primitive (-a,b,c), 1<=a,b<c<=60; all coefficients nonzero.
The gcd census is FINITE-EXACT, not an all-exponent noncancellation proof.
"""
from collections import Counter
from math import factorial, gcd
import sympy as sp

BOUND = 60
T = sp.Symbol("T")


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def direct_rows(a, b, c, m):
    """Independent charge-balance implementation, scanning negative count."""
    result = []
    for x in range(m + 1):
        numerator = (a + c) * x - c * m
        if numerator % (b - c):
            continue
        y = numerator // (b - c)
        z = m - x - y
        if y >= 0 and z >= 0:
            result.append((x, y, z))
    return sorted(result)


def coordinates(a, b, c):
    g = gcd(a + b, a + c)
    A, B = (a + b) // g, (a + c) // g
    z0 = 0 if A == 1 else (a * pow(B, -1, A)) % A
    y0 = (a - B * z0) // A
    return g, A, B, y0, z0


def in_semigroup(n, A, B):
    if n < 0:
        return False
    z = 0 if A == 1 else (n * pow(B, -1, A)) % A
    return n >= B * z


def compressed_moment(rows, B):
    """Remove a torus monomial; T is the middle-coefficient B-th power."""
    m = sum(rows[0])
    ymin = min(y for _, y, _ in rows)
    coeffs = {}
    for x, y, z in rows:
        require((y - ymin) % B == 0, '(y - ymin) % B == 0')
        coeffs[((y - ymin) // B,)] = factorial(m) // (
            factorial(x) * factorial(y) * factorial(z))
    return sp.Poly.from_dict(coeffs, T, domain=sp.QQ).monic()


def raw_moment(a, b, c, m):
    return sp.Poly.from_dict({(y,): factorial(m) // (
        factorial(x) * factorial(y) * factorial(z))
        for x, y, z in direct_rows(a, b, c, m)}, T, domain=sp.QQ)


def main():
    counts = Counter()
    maximum_channels = (0, None)
    for a in range(1, BOUND + 1):
        for b in range(1, BOUND + 1):
            for c in range(b + 1, BOUND + 1):
                if gcd(a, gcd(b, c)) != 1:
                    continue
                counts["primitive_patterns"] += 1
                g, A, B, y0, z0 = coordinates(a, b, c)
                pair_bound = min((a + b) // gcd(a, b),
                                 (a + c) // gcd(a, c))
                k0 = next(k for k in range(1, pair_bound // g + 1)
                          if in_semigroup(a * k, A, B))
                predicted_m = g * k0
                # Direct scan checks every earlier m, including nonmultiples of g.
                for m in range(1, predicted_m + 1):
                    rr = direct_rows(a, b, c, m)
                    require(bool(rr) == (m == predicted_m), 'bool(rr) == (m == predicted_m)')
                collision = len(rr) > 1
                require(collision == in_semigroup(a - A * B, A, B), 'collision == in_semigroup(a - A * B, A, B)')
                require(collision == (y0 >= B), 'collision == (y0 >= B)')
                counts["collided" if collision else "singleton"] += 1
                if not collision:
                    require(predicted_m <= a + c, 'predicted_m <= a + c')
                    continue
                require(predicted_m == g and g > B and (2 * g <= a + c), 'predicted_m == g and g > B and (2 * g <= a + c)')
                q, residue = divmod(y0, B)
                require(len(rr) == q + 1, 'len(rr) == q + 1')
                if len(rr) > maximum_channels[0]:
                    maximum_channels = (len(rr), (a, b, c, g))
                first = set(rr)
                for k in range(1, 6):
                    rows_k = direct_rows(a, b, c, k * g)
                    expected = k * q + (k * residue) // B + (k * z0) // A + 1
                    require(len(rows_k) == expected, 'len(rows_k) == expected')
                    counts["carry_profile_checks"] += 1
                doubled = {tuple(v[i] + w[i] for i in range(3))
                           for v in first for w in first}
                actual_second = set(direct_rows(a, b, c, 2 * g))
                require(doubled <= actual_second, 'doubled <= actual_second')
                require(len(doubled) == 2 * q + 1, 'len(doubled) == 2 * q + 1')
                excess = (2 * residue) // B + (2 * z0) // A
                require(len(actual_second - doubled) == excess, 'len(actual_second - doubled) == excess')
                counts[f"second_rung_excess_{excess}"] += 1
                P = compressed_moment(rr, B)
                Q = compressed_moment(sorted(actual_second), B)
                require(sp.gcd(P, Q).degree() == 0, 'sp.gcd(P, Q).degree() == 0')
                counts["two_rung_gcd_one"] += 1

    print("UNIVERSE: primitive (-a,b,c), 1<=a<=60, 1<=b<c<=60")
    for key in sorted(counts):
        print(f"{key}: {counts[key]}")
    print("maximum_first_channels:", maximum_channels)

    print("CONTROL: excluded old c=1 classification, charges (-3,1,9)")
    P, Q = raw_moment(3, 1, 9, 4), raw_moment(3, 1, 9, 8)
    require(sp.rem(Q, P) == -224, 'sp.rem(Q, P) == -224')
    print("moments:", P.as_expr(), ";", Q.as_expr(), "; remainder:", sp.rem(Q, P).as_expr())

    print("HOSTILE: two first channels do not imply three second channels")
    first = set(direct_rows(13, 1, 8, 7))
    second = set(direct_rows(13, 1, 8, 14))
    doubled = {tuple(v[i] + w[i] for i in range(3)) for v in first for w in first}
    require(len(first) == 2 and len(second) == 5, 'len(first) == 2 and len(second) == 5')
    print("(-13,1,8): first=", sorted(first), "; new second channels=", sorted(second - doubled))

    print("HOSTILE FAMILY: unbounded shortest signed-relation height")
    for A, B in [(2, 3), (3, 4), (4, 5), (7, 8)]:
        a, g = A * B, A * B + 1
        require(g > 2 * B, 'g > 2 * B')
        b, c = g * A - a, g * B - a
        relation = (B - A, -B, A)
        require(-a * relation[0] + b * relation[1] + c * relation[2] == 0, '-a * relation[0] + b * relation[1] + c * relation[2] == 0')
        require(len(direct_rows(a, b, c, g)) == 2, 'len(direct_rows(a, b, c, g)) == 2')
        # Independent finite shortest-relation audit of displayed witnesses.
        found = []
        for r in range(-2 * B, 2 * B + 1):
            for s in range(-2 * B, 2 * B + 1):
                numerator = a * r - b * s
                if numerator % c:
                    continue
                t = numerator // c
                norm = abs(r) + abs(s) + abs(t)
                if norm and norm <= 2 * B:
                    found.append((norm, (r, s, t)))
        require(min((n for n, _ in found)) == 2 * B, 'min((n for n, _ in found)) == 2 * B')
        print(f"charges={(-a,b,c)}, m0={g}, min_relation_l1={2*B}, relation={relation}")

    print("BOUNDARY: neutral coefficient detects at 1; zero coefficients reduce support;")
    print("reflection permutes rows; common charge scaling preserves every row.")
    for a, b, c in [(2, 1, 4), (3, 1, 9), (13, 1, 8)]:
        for m in range(1, 16):
            require(direct_rows(2 * a, 2 * b, 2 * c, m) == direct_rows(a, b, c, m), 'direct_rows(2 * a, 2 * b, 2 * c, m) == direct_rows(a, b, c, m)')
    print("ALL EXACT ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
