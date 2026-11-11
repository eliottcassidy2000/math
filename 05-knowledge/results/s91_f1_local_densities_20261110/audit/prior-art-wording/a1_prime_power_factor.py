#!/usr/bin/env python3
"""Auditor prior-art-wording, check a1.

Question: is the local factor of F1 at an odd prime p a function of [p | N] only (README section 2, THM-5128 (D3)
"singular series"), or does it depend on the p-adic valuation of N?

(1) Uniform-state model: u, v independent uniform on (Z/p^2)^2, form u1 v1 + c u2 v2 (c a unit).  Exact class law
    modulo p^2 by brute force, compared with the closed form
        p^2 P(B = c) = (1 + 1/p)(1 - p^-(v+1))         for v = v_p(c) in {0, 1},
        p^2 P(B = 0) = (1 + 1/p)(1 - p^-2) + p^-2.
(2) F1 itself (own enumeration, ordered triples (h1, h2, s), N = x1 x2 + 2^s y1 y2 - [s = 0]) on a window of odd N:
    class means of the representation count modulo 25, 9 and 49, split by sheet and by valuation.
"""
import sys
import numpy as np
from numba import njit
from fractions import Fraction


def model(p, c):
    q = p * p
    cnt = [0] * q
    for u1 in range(q):
        for u2 in range(q):
            for v1 in range(q):
                a = u1 * v1
                for v2 in range(q):
                    cnt[(a + c * u2 * v2) % q] += 1
    tot = q ** 4
    return [Fraction(x * q, tot) for x in cnt]


def closed(p, v):
    P = Fraction(p)
    if v == 2:
        return (1 + 1 / P) * (1 - P ** -2) + P ** -2
    return (1 + 1 / P) * (1 - P ** -(v + 1))


@njit(cache=False)
def table(nmax):
    psi = np.zeros(nmax, dtype=np.int64)
    for h in range(1, nmax):
        if h % 2 == 0:
            psi[h] = psi[h // 2] + h // 2
        else:
            k = h // 2
            psi[h] = psi[k] + k // 2 + 1
    return psi


@njit(cache=False)
def enumerate_window(x, y, A, B, r0, r1):
    n = len(x)
    # sheet 0: N = x1 x2 + y1 y2 - 1 ; sheet 1: N = x1 x2 + 2^s y1 y2
    for h1 in range(n):
        y1 = y[h1]
        x1 = x[h1]
        if y1 * 1 >= B + 1:
            break
        for h2 in range(n):
            y2 = y[h2]
            if y1 * y2 >= B + 1:
                break
            x2 = x[h2]
            xx = x1 * x2
            N = xx + y1 * y2 - 1
            if A <= N < B:
                r0[(N - A) // 2] += 1
            m = 2 * y1 * y2
            while m < B:
                N = xx + m
                if A <= N < B:
                    r1[(N - A) // 2] += 1
                m *= 2


def main():
    ok = True
    print("(1) uniform-state model modulo p^2 (exact rationals)")
    for p, c in ((3, 1), (3, 2), (5, 1), (5, 2), (5, 4)):
        law = model(p, c)
        q = p * p
        v0 = {law[r] for r in range(q) if r % p != 0}
        v1 = {law[r] for r in range(q) if r % p == 0 and r != 0}
        good = v0 == {closed(p, 0)} and v1 == {closed(p, 1)} and law[0] == closed(p, 2)
        ok &= good
        modp = Fraction(sum(law[r] for r in range(0, q, p)), p)
        print(f"    p = {p}, c = {c}: factor at units {float(closed(p, 0)):.5f}, at v = 1 {float(closed(p, 1)):.5f}, "
              f"at 0 mod p^2 {float(closed(p, 2)):.5f}; mean over the class 0 mod p = {float(modp):.5f} "
              f"(1 - 1/p^2 + 1/p = {1 - 1 / p**2 + 1 / p:.5f});  closed form agrees: {good}")
    kA, kB = 22, 24
    if len(sys.argv) > 2:
        kA, kB = int(sys.argv[1]), int(sys.argv[2])
    A, B = (1 << kA) + 1, (1 << kB)
    nmax = (B // 2) + 2
    psi = table(nmax)
    h = np.arange(nmax, dtype=np.int64)
    x = 2 * psi + 1
    y = 2 * h + 1
    assert x[1] == 3 and x[2] == 5 and x[3] == 5 and x[10] == 19 and x[13] == 19 and x[17] == 27 and x[256] == 513
    nodd = (B - A + 1) // 2
    r0 = np.zeros(nodd, dtype=np.int32)
    r1 = np.zeros(nodd, dtype=np.int32)
    enumerate_window(x, y, A, B, r0, r1)
    Ns = A + 2 * np.arange(nodd, dtype=np.int64)
    print(f"(2) F1, odd N in [2^{kA}, 2^{kB}): {nodd} integers; triples per odd N: sheet 0 {r0.mean():.4f}, "
          f"sheet 1 {r1.mean():.4f}, total {(r0 + r1).mean():.4f}; misses {(r0 + r1 == 0).mean() * 100:.3f} %")
    for p in (5, 3, 7):
        q = p * p
        for name, r, shift in (("sheet 1 (class of N)", r1, 0), ("sheet 0 (class of N + 1)", r0, 1)):
            M = (Ns + shift)
            mean = r.mean()
            f_unit = r[M % p != 0].mean() / mean
            f_v1 = r[(M % p == 0) & (M % q != 0)].mean() / mean
            f_v2 = r[M % q == 0].mean() / mean
            f_p = r[M % p == 0].mean() / mean
            print(f"    p = {p}, {name}: factor at p-units {f_unit:.4f}; p | . exactly {f_v1:.4f}; p^2 | . {f_v2:.4f}; "
                  f"(all of the class 0 mod p: {f_p:.4f});  ratio (p^2 | .)/(p || .) = {f_v2 / f_v1:.4f}; "
                  f"limit model {float(closed(p, 2) / closed(p, 1)):.4f}")
    print("ALL CHECKS PASSED" if ok else "CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
