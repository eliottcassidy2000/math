#!/usr/bin/env python3
"""Recount, by a subset dynamic programme (a routine independent of audit_new_items.py's depth-first tables), of the
second-order-gap witnesses found by audit_new_items.py's random search, and of the packet's order-8 witness.
For each: c = H(R), d_a, d_b, mu and the ratio mu / max(d).  The ratios found at orders 7..10 are
(2^j + 1)/(3 2^j + 1), j = 0..3, i.e. 1/2, 3/7, 5/13, 9/25 (the packet reports 14/34 as the sampled minimum at order 10).
"""
from fractions import Fraction


def H(out, mask):
    n = len(out)
    if mask == 0:
        return 1
    cnt = {}
    vs = [v for v in range(n) if mask >> v & 1]
    for v in vs:
        cnt[(1 << v, v)] = 1
    total = 0
    # process masks in increasing order of popcount
    layer = dict(cnt)
    for _ in range(len(vs) - 1):
        nxt = {}
        for (m, v), c in layer.items():
            cand = out[v] & mask & ~m
            while cand:
                low = cand & -cand
                cand ^= low
                key = (m | low, low.bit_length() - 1)
                nxt[key] = nxt.get(key, 0) + c
        layer = nxt
    return sum(layer.values())


def report(name, out, a, b):
    n = len(out)
    assert all(((out[i] >> j) & 1) + ((out[j] >> i) & 1) == 1 for i in range(n) for j in range(i + 1, n))
    full = (1 << n) - 1
    R = full ^ (1 << a) ^ (1 << b)
    c = H(out, R)
    da = H(out, R | (1 << a)) - c
    db = H(out, R | (1 << b)) - c
    mu = H(out, full) - c - da - db
    print(f"{name}: order {n}, masks {out}, a = {a}, b = {b}: (c, d_a, d_b, mu) = {(c, da, db, mu)}; "
          f"mu/max(d) = {Fraction(mu, max(da, db))}; 2 mu >= max(d): {2 * mu >= max(da, db)}")
    return Fraction(mu, max(da, db))


def main():
    r8p = report("packet order-8 witness", [246, 252, 248, 177, 160, 64, 152, 32], 1, 2)
    r8 = report("audit order-8 witness ", [250, 204, 249, 224, 234, 66, 128, 32], 3, 6)
    r9 = report("audit order-9 witness ", [494, 508, 504, 496, 449, 80, 256, 352, 32], 6, 7)
    r10 = report("audit order-10 witness", [894, 1020, 1016, 1008, 992, 448, 768, 833, 512, 32], 6, 8)
    ok = (r8p == Fraction(3, 7) and r8 == Fraction(3, 7) and r9 == Fraction(5, 13) and r10 == Fraction(9, 25))
    print("ratios 3/7, 3/7, 5/13, 9/25 = (2^j + 1)/(3 2^j + 1) for j = 1, 1, 2, 3:", "OK" if ok else "FAIL")
    print("9/25 < 14/34:", Fraction(9, 25) < Fraction(14, 34))
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
