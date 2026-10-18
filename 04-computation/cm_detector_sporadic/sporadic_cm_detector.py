#!/usr/bin/env python3
"""S43's p-adic CM detector (second-order Lucas-block congruence) applied to SPORADIC Apery-like families.
kind-pasteur-2026-10-18-S2.  See PREREGISTRATION.md (written and pushed before this script was run).

Test (B1'), t_k = u_k m^{-k}:   sum_{k=p}^{2p-1} t_k == t_1 sum_{k<p} t_k   (mod p^2),
filter primes 7..43, confirm primes 47..200; m ranges over 0 < |m| <= M (m with p | m count as passing that prime).
Control: supersingular density = fraction of primes 7 <= p <= 500 (p not dividing m) with sum_{k<p} u_k m^{-k} == 0 mod p.
Usage: python sporadic_cm_detector.py [M]
"""
import sys, json
from fractions import Fraction
from math import comb
import numpy as np

NMAX = 520   # need u_k for k < 2p <= 400 (scan) and k < 500 (density control)


def primes_between(a, b):
    return [p for p in range(a, b + 1) if p > 1 and all(p % d for d in range(2, int(p ** 0.5) + 1))]


def zagier_rec(a, b, c, N=NMAX):
    """(n+1)^2 u_{n+1} = (a n^2 + a n + b) u_n - c n^2 u_{n-1}, u_0 = 1, u_1 = b (exact; integrality checked)."""
    u = [Fraction(1), Fraction(b)]
    for n in range(1, N):
        nxt = ((a * n * n + a * n + b) * u[n] - c * n * n * u[n - 1]) / Fraction((n + 1) ** 2)
        u.append(nxt)
    assert all(x.denominator == 1 for x in u), "non-integral sequence: parameters wrong"
    return [int(x) for x in u]


def from_sum(f, N=NMAX):
    return [f(n) for n in range(N + 1)]


FAMILIES = {
    # zeta(2)-type sporadic (Zagier 2009)
    'Zagier A (Franel)': lambda: zagier_rec(7, 2, -8),
    'Zagier B': lambda: zagier_rec(9, 3, 27),
    'Zagier C': lambda: zagier_rec(10, 3, 9),
    'Zagier D (Apery zeta2)': lambda: zagier_rec(11, 3, -1),
    'Zagier E': lambda: zagier_rec(12, 4, 32),
    'Zagier F': lambda: zagier_rec(17, 6, 72),
    # zeta(3)-type
    'Apery zeta3': lambda: from_sum(lambda n: sum(comb(n, k) ** 2 * comb(n + k, k) ** 2 for k in range(n + 1))),
    'Domb': lambda: from_sum(lambda n: sum(comb(n, k) ** 2 * comb(2 * k, k) * comb(2 * n - 2 * k, n - k) for k in range(n + 1))),
    'Cooper s7': lambda: from_sum(lambda n: sum(comb(n, k) ** 2 * comb(n + k, k) * comb(2 * k, n) for k in range(n + 1))),
    'Cooper s10': lambda: from_sum(lambda n: sum(comb(n, k) ** 4 for k in range(n + 1))),
}


def lucas_test(u, ms, p):
    pp = p * p
    A = [u[k] % pp for k in range(2 * p)]
    r = np.mod(ms, pp)
    invtab = np.zeros(pp, dtype=np.int64)
    for x in range(pp):
        if x % p:
            invtab[x] = pow(x, -1, pp)
    w = invtab[r]
    S0 = np.zeros_like(ms); S1 = np.zeros_like(ms); pw = np.ones_like(ms)
    for k in range(2 * p):
        if k < p:
            S0 = (S0 + A[k] * pw) % pp
        else:
            S1 = (S1 + A[k] * pw) % pp
        pw = (pw * w) % pp
    ok = (S1 - (A[1] * w % pp) * S0) % pp == 0
    return ok | (r % p == 0)


def ss_density(u, m, primes):
    z = tot = 0
    for p in primes:
        if m % p == 0:
            continue
        im = pow(m % p, -1, p)
        s = 0; w = 1
        for k in range(p):
            s = (s + (u[k] % p) * w) % p
            w = w * im % p
        tot += 1; z += (s == 0)
    return z / tot if tot else float('nan')


def main():
    M = int(sys.argv[1]) if len(sys.argv) > 1 else 100000
    ms_all = np.array([m for m in range(-M, M + 1) if m != 0], dtype=np.int64)
    filt, conf, ssp = primes_between(7, 43), primes_between(47, 200), primes_between(7, 500)
    out = {}
    for name, gen in FAMILIES.items():
        u = gen()
        ms = ms_all.copy()
        for p in filt:
            ms = ms[lucas_test(u, ms, p)]
        survivors = ms.copy()
        for p in conf:
            survivors = survivors[lucas_test(u, survivors, p)]
        surv = sorted(int(x) for x in survivors)
        dens = {m: round(ss_density(u, m, ssp), 3) for m in surv}
        ctrl = {m: round(ss_density(u, m, ssp), 3) for m in (2, 3, 5, 7, 11, 100, -100, 999) if m not in surv}
        out[name] = {'survivors': surv, 'after_filter_count': int(len(ms)), 'ss_density': dens, 'controls': ctrl,
                     'u_head': [int(x) for x in u[:8]]}
        print(f"{name}: u = {u[:7]}...", flush=True)
        print(f"   survivors (|m| <= {M}, p <= 200): {surv}   [after p<=43 filter: {len(ms)}]")
        print(f"   supersingular density of survivors: {dens}")
        print(f"   controls: {ctrl}", flush=True)
    json.dump(out, open(f'sporadic_cm_detector_M{M}.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
