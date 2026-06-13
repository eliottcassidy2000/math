#!/usr/bin/env python3
"""EXTENSION of THM-398 headline claim to n=19,20 (monad-compute-2026-06-03-S596).

Continuation of lrc_multiple_loose_extended_monad_s595.py, which verified the
THM-398 / HYP-2102 reduction claim "every multiple-of-n speed set is LOOSE
(M(S) > 1/n), 0 tight-with-mult" exactly through n=18. The S595 handoff noted the
next extension to n>=19 is straightforward (cost grows with denominators).

This script pushes the direct EXACT looseness verification to n=19,20, with
n=16,17,18 reproduced as control. Same definitions and method as S595:

LOOSENESS TEST (exact, rational): a speed set S on n runners is LOOSE iff
M(S) = sup_t min_i ||v_i t|| > 1/n, equivalently iff the OPEN safe set
  A_open = { t in [0,1) : min_i ||v_i t|| > 1/n }
has positive measure. Computed exactly with Fraction arithmetic. Breakpoints of
t |-> 1[min_i ||v_i t|| > 1/n] are exactly t = (k*n +- 1)/(n*v_i); each elementary
interval is decided by its midpoint via a STRICT inequality.

If mu(A_open) > 0 for every multiple-of-n config tested, C' has no tight witness
at that n. HARDEST slice = v = n exactly (THM-398) with small companions; plus a
random slice v = n*w, w in {1,2,3}, wider companions.

No proof; pure computation.
"""
from fractions import Fraction as F
from math import gcd, comb
import itertools, random


def dist(x):
    x %= 1
    return min(x, 1 - x)


def open_safe_measure(V, n):
    """Exact measure of { t : min_i ||v_i t|| > 1/n } (STRICT)."""
    THR = F(1, n)
    eps = {F(0)}
    for v in V:
        for k in range(v + 1):
            for s in (-1, 1):
                eps.add(F(k * n + s, n * v) % 1)
    pts = sorted(eps)
    meas = F(0)
    L = len(pts)
    for i in range(L):
        a = pts[i]
        b = pts[(i + 1) % L]
        ln = (b - a) if b > a else (b - a + 1)
        mid = (a + ln / 2) % 1
        if all(dist(v * mid) > THR for v in V):  # STRICT => open safe set
            meas += ln
    return meas


def prim(V):
    g = 0
    for v in V:
        g = gcd(g, v)
    return tuple(sorted(v // g for v in V))


def check_n(n, n_random=800, seed=1234):
    """Return (loose, total, n_tight, min_measure) over sampled+hardest
    multiple-of-n configs on n runners."""
    m = n - 1  # number of speeds (runners excluding the stationary observer 0)
    rng = random.Random(seed + n)
    configs = set()

    # --- HARDEST systematic slice: v = n exactly, small companions in [1, n+6] ---
    pool = [x for x in range(1, n + 7) if x % n != 0]
    if comb(len(pool), m - 1) <= 4000:
        for comb_others in itertools.combinations(pool, m - 1):
            V = prim(tuple(sorted(set(comb_others) | {n})))
            if len(V) == m and any(x % n == 0 for x in V):
                configs.add(V)
    else:
        tries = 0
        while len(configs) < 2000 and tries < 60000:
            tries += 1
            others = rng.sample(pool, m - 1)
            V = prim(tuple(sorted(set(others) | {n})))
            if len(V) == m and any(x % n == 0 for x in V):
                configs.add(V)

    # --- RANDOM general slice: v = n*w, w in {1,2,3}, wider companions ---
    target = len(configs) + n_random
    tries = 0
    wide = [x for x in range(1, n + 12) if x % n != 0]
    while len(configs) < target and tries < 120000:
        tries += 1
        w = rng.randint(1, 3)
        v = n * w
        others = rng.sample(wide, m - 1)
        if v in others:
            continue
        V = prim(tuple(sorted(set(others) | {v})))
        if len(V) == m and any(x % n == 0 for x in V):
            configs.add(V)

    loose = 0
    tight = 0
    min_meas = None
    for V in configs:
        mu = open_safe_measure(V, n)
        if mu > 0:
            loose += 1
            if min_meas is None or mu < min_meas:
                min_meas = mu
        else:
            tight += 1
    return loose, len(configs), tight, min_meas


def main():
    print("THM-398 headline: every multiple-of-n speed set is LOOSE (M(S) > 1/n).")
    print("Direct EXACT looseness test via open-safe-set measure.")
    print("Control n=16,17,18 (reproduces S595); NEW extension n=19,20.")
    print()
    print(f"{'n':>3} {'parity':>6} {'#configs':>9} {'loose':>7} {'tight':>6} "
          f"{'all_loose':>10} {'min_margin(measure)':>22}")
    grand_loose = grand_tot = grand_tight = 0
    for n in [16, 17, 18, 19, 20]:
        loose, tot, tight, min_meas = check_n(n)
        parity = "odd" if n % 2 else "even"
        all_loose = (tight == 0 and loose == tot)
        mm = f"{float(min_meas):.5f} (={min_meas})" if min_meas is not None else "n/a"
        print(f"{n:>3} {parity:>6} {tot:>9} {loose:>7} {tight:>6} "
              f"{str(all_loose):>10} {mm:>22}", flush=True)
        grand_loose += loose
        grand_tot += tot
        grand_tight += tight
    print()
    print(f"TOTAL: {grand_loose}/{grand_tot} loose, {grand_tight} tight.")
    verdict = ("CONFIRMED: 0 tight multiple-of-n configs through n=20 "
               "(C' holds for tested configs n<=20)"
               if grand_tight == 0 else
               "*** TIGHT MULTIPLE-OF-N CONFIG FOUND -- would REFUTE C' ***")
    print(verdict)


if __name__ == '__main__':
    main()
