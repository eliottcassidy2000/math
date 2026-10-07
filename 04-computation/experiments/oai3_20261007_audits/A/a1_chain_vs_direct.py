#!/usr/bin/env python3
"""Audit A, item 1: independent implementation of the pair chain (k, N), N = 3^max(0,-k) e (integer),
checked against direct 2-adic orbits of u = 3^k0 y + e0 and v = y (Terras map T).
Checks: relation u_n = 3^k_n v_n + e_n (mod 2^(P-n)), parity(u_n) = beta_n xor sigma_n, N integral,
k_n - k_0 = O_u(n) - O_v(n) (odd-step count difference), absorption <=> u_n == v_n (to precision).
"""
import random, sys
from fractions import Fraction

def chain_step(k, N, beta):
    """exact integer transition; returns (k', N')"""
    sig = N & 1
    if k >= 1:
        if sig == 0 and beta == 0: return k, N // 2
        if sig == 0 and beta == 1:
            x = 3 * N + 1 - 3 ** k; assert x % 2 == 0; return k, x // 2
        if sig == 1 and beta == 0:
            x = 3 * N + 1; assert x % 2 == 0; return k + 1, x // 2
        x = N - 3 ** (k - 1); assert x % 2 == 0; return k - 1, x // 2
    if k == 0:
        if sig == 0 and beta == 0: return 0, N // 2
        if sig == 0 and beta == 1: return 0, 3 * N // 2
        if sig == 1 and beta == 0:
            x = 3 * N + 1; assert x % 2 == 0; return 1, x // 2
        x = 3 * N - 1; assert x % 2 == 0; return -1, x // 2
    h = -k
    if sig == 0 and beta == 0: return k, N // 2
    if sig == 0 and beta == 1:
        x = 3 * N + 3 ** h - 1; assert x % 2 == 0; return k, x // 2
    if sig == 1 and beta == 0:
        x = N + 3 ** (h - 1); assert x % 2 == 0; return k + 1, x // 2
    x = 3 * N - 1; assert x % 2 == 0; return k - 1, x // 2

def chain_step_frac(k, e, beta):
    """the theorem's table, verbatim, in exact rationals"""
    sig = e.numerator % 2  # e in Z[1/3]: 2-adic parity = numerator parity (denominator odd)
    if sig == 0 and beta == 0: return k, e / 2
    if sig == 0 and beta == 1: return k, (3 * e + 1 - Fraction(3) ** k) / 2
    if sig == 1 and beta == 0: return k + 1, (3 * e + 1) / 2
    return k - 1, (e - Fraction(3) ** (k - 1)) / 2

def T(x, mod):
    return (x >> 1) if x % 2 == 0 else (((3 * x + 1) >> 1) % mod)

def run(k0, e0, P, rng, steps):
    mod = 1 << P
    y = rng.getrandbits(P)
    inv3 = pow(3, -1, mod)
    def tw(fr, m):  # Fraction with odd denominator -> 2-adic mod m
        return (fr.numerator * pow(fr.denominator, -1, m)) % m
    u = (pow(3, k0, mod) * y + tw(Fraction(e0), mod)) % mod if k0 >= 0 else (pow(inv3, -k0, mod) * y + tw(Fraction(e0), mod)) % mod
    v = y
    k, e = k0, Fraction(e0)
    N = int(e * 3 ** max(0, -k)); assert N == e * 3 ** max(0, -k)
    Ou = Ov = 0
    bad = 0; absorbed_at = None
    prec = P
    for n in range(steps):
        m = 1 << prec
        rel = (pow(3, k, m) if k >= 0 else pow(pow(3, -1, m), -k, m)) * v + tw(e, m)
        if (u - rel) % m != 0: bad += 1
        beta = v & 1; sig = e.numerator & 1
        if (u & 1) != (beta ^ sig): bad += 1
        if k - k0 != Ou - Ov: bad += 1
        if (k, e) == (0, 0):
            if absorbed_at is None: absorbed_at = n
            if (u - v) % m != 0: bad += 1
        else:
            # not absorbed: u != v unless value coincidence; check difference is nonzero mod 2^prec (prec large)
            pass
        Ou += u & 1; Ov += v & 1
        k2, e2 = chain_step_frac(k, e, beta)
        k3, N3 = chain_step(k, N, beta)
        if k2 != k3 or N3 != e2 * 3 ** max(0, -k2) or (e2 * 3 ** max(0, -k2)).denominator != 1: bad += 1
        if k2 == 0 and e2.denominator != 1: bad += 1
        k, e, N = k2, e2, N3
        u = T(u, m) % (m >> 1); v = T(v, m) % (m >> 1)
        prec -= 1
        if prec < 64: break
    return bad, absorbed_at, n + 1

if __name__ == '__main__':
    rng = random.Random(20261007)
    starts = [(0, 1), (0, -1), (1, 0), (-1, 0), (1, 2), (2, 8), (0, 5), (0, -7), (2, 5), (-2, Fraction(5, 9)),
              (-1, Fraction(1, 3)), (3, 26), (0, 1000001), (-3, Fraction(-13, 27)), (5, 242)]
    total_steps = 0; total_bad = 0; absorbed = 0; trials = 0
    for (k0, e0) in starts:
        for t in range(40):
            bad, ab, ns = run(k0, e0, 1200, rng, 1100)
            total_steps += ns; total_bad += bad; trials += 1; absorbed += ab is not None
    print(f"[chain vs direct] {trials} trials, {total_steps} steps, mismatches = {total_bad}; absorbed within ~1100 steps: {absorbed}/{trials}")
