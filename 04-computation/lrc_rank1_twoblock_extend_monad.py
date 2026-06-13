#!/usr/bin/env python3
"""Extend the bounded-CRT-automaton emptiness check (lrc_rank1_twoblock_s595.py)
to larger n and more trials.

Frontier result (opus-S595): for multiple-of-n residual rows, the per-component
"allowed-w" bounded CRT automaton has EMPTY intersection (=> loose) for 400/400
random configs at n=10,12,14. The open question is whether this holds for ALL
large-owner residuals at every n. This script pushes the empirical envelope:
larger n (16,18,20,22,24) and 1000 trials per n.

monad-compute-2026-06-03-S1. Pure compute, no proof. Logic copied verbatim from
the S595 automaton (same G_components, same allowed-w transition) so the numbers
are directly comparable.
"""
from fractions import Fraction as F
from math import gcd
import random

def dist(x):
    x %= 1
    return min(x, 1 - x)

def G_components(Sp, n):
    THR = F(1, n); pts = {}
    for u in Sp:
        for k in range(u + 1):
            for eps in (1, -1):
                pts.setdefault(F(k * n + eps, n * u) % 1, []).append((u, k, eps))
    order = sorted(pts); comps = []; L = len(order)
    for i in range(L):
        a = order[i]; b = order[(i + 1) % L]
        ln = (b - a) if b > a else (b - a + 1); mid = (a + ln / 2) % 1
        if all(dist(u * mid) > THR for u in Sp):
            comps.append((a, b, ln, pts[a], pts[b]))
    return comps

def run(n, trials, rng):
    m = n - 1; B = 2 * n + 8; emptyc = 0; tot = 0
    for _ in range(trials):
        others = set(rng.sample([x for x in range(1, B + 1) if x % n != 0], m - 1))
        ww = rng.randint(1, 3); v = n * ww
        if v in others: continue
        V = tuple(sorted(others | {v}))
        g = 0
        for x in V: g = gcd(g, x)
        V = tuple(sorted(x // g for x in V))
        if len(V) != m or not any(x % n == 0 for x in V): continue
        mults = [x for x in V if x % n == 0]; vv = mults[0]
        Sp = tuple(x for x in V if x != vv)
        comps = G_components(Sp, n)
        if not comps: continue
        tot += 1
        wmax = ((n - 1) * max(Sp)) // n + 1
        allowed = set(range(1, min(wmax, 200) + 1))
        for (a, b, ln, oa, ob) in comps:
            mid = (a + ln / 2) % 1
            allowed = {w for w in allowed if dist(n * w * mid) <= F(1, n) - F(n * w, 2) * ln}
            if not allowed: break
        if not allowed: emptyc += 1
    return tot, emptyc

def main():
    print("Bounded-CRT-automaton emptiness (large-owner residual => loose), extended.")
    print("EMPTY/tot = configs where per-component allowed-w intersection is empty.")
    print("=" * 64)
    rng = random.Random(1)
    grand_tot = 0; grand_empty = 0
    for n in [10, 12, 14, 16, 18, 20, 22, 24]:
        tot, emptyc = run(n, 1000, rng)
        grand_tot += tot; grand_empty += emptyc
        flag = "ALL EMPTY" if emptyc == tot else f"*** {tot - emptyc} NON-EMPTY ***"
        print(f"   n={n:3d}: residual rows={tot:4d}; automaton EMPTY={emptyc:4d}/{tot:<4d}  {flag}")
    print("=" * 64)
    print(f"   TOTAL across n=10..24: {grand_empty}/{grand_tot} empty "
          f"({'NO COUNTEREXAMPLE' if grand_empty == grand_tot else 'COUNTEREXAMPLE FOUND'})")

if __name__ == '__main__':
    main()
