#!/usr/bin/env python3
"""Independent hull/bitset audit for the finite Newton-barcode claim."""

from math import gcd, isqrt

import gmpy2
from flint import fmpz_poly


V = fmpz_poly([0, 1])


def prime(n):
    if n < 2:
        return False
    q = 2
    while q * q <= n:
        if n % q == 0:
            return n == q
        q += 1
    return True


def factor_count(n):
    count = 0
    q = 2
    while q * q <= n:
        if n % q == 0:
            count += 1
            while n % q == 0:
                n //= q
        q += 1
    return count + (n > 1)


def is_residual(d):
    return not prime(d) and not prime(d - 2) and factor_count(d - 1) != 1


def moments(d):
    previous = fmpz_poly([1])
    current = fmpz_poly([d - 1, 2])
    d_to_k = d
    for k in range(1, d - 1):
        nxt = (
            fmpz_poly([d_to_k * (d - k - 1)])
            + 2 * (k + 1) * (2 * k + 1) * V * current
            + k * (k + 1) * (previous - 4 * d * V * previous)
        )
        previous, current = current, nxt
        d_to_k *= d
    return previous, current


def first_remainder(p, q, d):
    n = d - 2
    return (2 * n - 1) * (q - 2 * (n + 1) * (2 * n + 1) * V * p) + 2 * d * (n + 1) * p


def valuation(a, p):
    a = gmpy2.mpz(a)
    if not a:
        return None
    return int(gmpy2.remove(a, p)[1])


def profile(poly, p):
    points = []
    zero_order = None
    for x in range(poly.degree() + 1):
        y = valuation(poly[x], p)
        if y is None:
            continue
        if zero_order is None:
            zero_order = x
        while len(points) >= 2:
            x0, y0 = points[-2]
            x1, y1 = points[-1]
            # Pop when the new slope is no larger than the preceding slope.
            if (y1 - y0) * (x - x1) >= (y - y1) * (x1 - x0):
                points.pop()
            else:
                break
        points.append((x, y))
    caps = {}
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        dy, dx = y1 - y0, x1 - x0
        g = gcd(abs(dy), dx)
        slope = (dy // g, dx // g)
        assert slope not in caps
        caps[slope] = dx
    return zero_order, caps


def allowed(rows, p):
    profiles = [profile(row, p) for row in rows]
    possibilities = set(range(min(z for z, _ in profiles) + 1))
    common = set(profiles[0][1])
    for _, caps in profiles[1:]:
        common &= set(caps)
    blocks = []
    for slope in sorted(common):
        denominator = slope[1]
        capacity = min(caps[slope] for _, caps in profiles)
        assert capacity % denominator == 0
        possibilities = {
            old + use
            for old in possibilities
            for use in range(0, capacity + 1, denominator)
        }
        blocks.append((slope, capacity, denominator))
    return possibilities, blocks


def main():
    bank = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47)
    ds = [d for d in range(1001, 1101) if is_residual(d)]
    assert len(ds) == 56
    closed = []
    witness_primes = {}
    for d in ds:
        p, q = moments(d)
        r = first_remainder(p, q, d)
        assert (p.degree(), q.degree(), r.degree()) == (d - 2, d - 1, d - 3)
        possible = set(range(1, d - 2))
        used = []
        for ell in bank:
            before = set(possible)
            possible &= allowed((p, q, r), ell)[0]
            if possible != before:
                used.append(ell)
            if not possible:
                break
        assert not possible, (d, sorted(possible))
        closed.append(d)
        witness_primes[d] = tuple(used)

    p, q = moments(1001)
    r = first_remainder(p, q, 1001)
    pair2, pair_blocks2 = allowed((p, q), 2)
    pair3, pair_blocks3 = allowed((p, q), 3)
    flag2, flag_blocks2 = allowed((p, q, r), 2)
    flag3, flag_blocks3 = allowed((p, q, r), 3)
    assert not ((pair2 & pair3) - {0})
    assert not ((flag2 & flag3) - {0})
    print("independent residual rows", len(closed), closed)
    print("d=1001 pair blocks", pair_blocks2, pair_blocks3)
    print("d=1001 flag blocks", flag_blocks2, flag_blocks3)
    print("witness prime traces", witness_primes)
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
