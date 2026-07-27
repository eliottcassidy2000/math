#!/usr/bin/env python3
"""Exact/arithmetic audit for the T(2,7) mirror-mixture selector packet."""

from fractions import Fraction


def add(a, b):
    out = [0] * max(len(a), len(b))
    for i, x in enumerate(a):
        out[i] += x
    for i, x in enumerate(b):
        out[i] += x
    return out


def mul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def det_recurrence(n):
    """D_n=(t-1)D_(n-1)+tD_(n-2), in ascending powers of t."""
    d0 = [1]
    if n == 0:
        return d0
    d1 = [-1, 1]
    if n == 1:
        return d1
    for _ in range(2, n + 1):
        d0, d1 = d1, add(mul([-1, 1], d1), mul([0, 1], d0))
    return d1


def inertia_at_r(r):
    """Exact inertia at theta=r*pi/7 for 0<=r<=14.

    Since sin(r*pi/14)=cos(|7-r|*pi/14) and cosine is strictly decreasing
    on [0,pi], the sign of sin(r*pi/14)-cos(j*pi/7) is the sign of
    2*j-|7-r|.
    """
    r = Fraction(r)
    assert 0 <= r <= 14
    angle_index = abs(Fraction(7) - r)
    signs = [Fraction(2 * j) - angle_index for j in range(1, 7)]
    pos = sum(value > 0 for value in signs)
    neg = sum(value < 0 for value in signs)
    zero = 6 - pos - neg
    return pos - neg, zero


def mixture_check():
    root_sig = [1, 3, 5, 5, 3, 1]
    off_sig = [0, 2, 4, 6, 4, 2, 0]
    for p in range(9):
        for q in range(9):
            n = p + q
            signed = p - q
            plus = [n + signed * s for s in root_sig]
            minus = [n - signed * s for s in root_sig]
            plus.extend(signed * s for s in off_sig)
            minus.extend(-signed * s for s in off_sig)
            max_plus = max(plus + [0])
            max_minus = max(minus + [0])
            mu = Fraction(max_plus + max_minus, 2)
            expected = n + 2 * abs(signed)
            assert mu == expected, (p, q, mu, expected)


def curvature_check():
    rows = []
    for delta in (Fraction(1), Fraction(2), Fraction(3), Fraction(7, 2)):
        t0 = delta / (8 - delta)
        m0 = (8 - delta) / 2
        a = m0 * (1 - t0)
        b = m0 * t0
        assert a == 4 - delta
        assert b == delta / 2
        rows.append((delta, t0, m0, a, b))
    return rows


def main():
    phi14 = [1, -1, 1, -1, 1, -1, 1]
    d6 = det_recurrence(6)
    assert d6 == phi14

    roots = [1, 3, 5, 9, 11, 13]
    root_table = [inertia_at_r(r) for r in roots]
    assert root_table == [(1, 1), (3, 1), (5, 1), (5, 1), (3, 1), (1, 1)]

    off_samples = [Fraction(1, 2), 2, 4, 7, 10, 12, Fraction(27, 2)]
    off_table = [inertia_at_r(r) for r in off_samples]
    assert off_table == [(0, 0), (2, 0), (4, 0), (6, 0), (4, 0), (2, 0), (0, 0)]

    mixture_check()

    # Phi_(zeta5,zeta1) and its mirror-swapped companion on (K,Kbar).
    gauge = ((3, -1), (-1, 3))
    determinant = gauge[0][0] * gauge[1][1] - gauge[0][1] * gauge[1][0]
    assert determinant == 8
    for p in range(9):
        for q in range(9):
            values = (3 * p - q, -p + 3 * q)
            lower = max(abs(values[0]), abs(values[1]))
            assert lower == 3 * max(p, q) - min(p, q)

    curvature_rows = curvature_check()

    print("T(2,7) BLANCHFIELD / GORDIAN SELECTOR AUDIT")
    print("det(tV-V^T) ascending coefficients:", d6)
    print("root r values:", roots)
    print("root (signature,nullity):", root_table)
    print("off-root sample (signature,nullity):", off_table)
    print("mixture BF formula checked for 0<=P,Q<=8: N+2|P-Q|")
    print("Gordian gauge matrix:", gauge, "determinant:", determinant)
    print("gauge norm checked for 0<=P,Q<=8: 3 max(P,Q)-min(P,Q)")
    print("curvature extremal rows (delta,t0,m0,A,B):")
    for row in curvature_rows:
        print(" ", tuple(str(x) for x in row))
    print("delta=4 endpoint: A<=0 forces open-interval curvature measure zero")
    print("PASS")


if __name__ == "__main__":
    main()
