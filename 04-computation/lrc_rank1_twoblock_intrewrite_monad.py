#!/usr/bin/env python3
"""Integer-arithmetic rewrite of the bounded-CRT-automaton emptiness check
(lrc_rank1_twoblock_extend_monad.py / lrc_rank1_twoblock_s595.py), so that
n >= 22 becomes feasible.

monad-compute-2026-06-03-S2. Pure compute, no proof.

WHY: the S595 automaton uses fractions.Fraction throughout G_components and the
allowed-w transition. Fraction reduces by gcd on every op, which makes the cost
super-linear in n (~95 s / 1000 trials already at n=20). All the rationals here
share the common denominator D = n * lcm(Sp); midpoints need 2D. So every
comparison can be done in exact integers over 2D, with NO gcd reductions.

INTEGER MODEL (denominator 2D, D = n*L, L = lcm(Sp)):
  - A breakpoint at fraction (k*n+eps)/(n*u) becomes the integer
        pint = ((k*n+eps) * (L//u)) % D            (numerator over D, in [0,D))
  - A gap (a,b) (a,b integers over D) has length lnint = (b-a) % D, and
    midpoint numerator over 2D:  midnum = 2*a + lnint   (in [0,2D)).
  - dist(u*mid) > 1/n   <=>   min(r, 2D-r) > 2L, where r = (u*midnum) % (2D).
        [since 2D/n = 2L]
  - keep automaton state w iff  dist(n*w*mid) <= 1/n - (n*w/2)*ln, i.e.
        min(r2, 2D-r2) <= 2L - n*w*lnint,  r2 = (n*w*midnum) % (2D).
    (RHS may be negative -> condition fails, since LHS >= 0.)

This is an EXACT reimplementation; for the cross-check primes it must reproduce
the Fraction version's (tot, EMPTY) counts bit-for-bit.
"""
from fractions import Fraction as F
from math import gcd
import random
import time


# ---------- original Fraction reference (verbatim from S595/S1) ----------
def dist_F(x):
    x %= 1
    return min(x, 1 - x)


def G_components_F(Sp, n):
    THR = F(1, n); pts = {}
    for u in Sp:
        for k in range(u + 1):
            for eps in (1, -1):
                pts.setdefault(F(k * n + eps, n * u) % 1, []).append((u, k, eps))
    order = sorted(pts); comps = []; L = len(order)
    for i in range(L):
        a = order[i]; b = order[(i + 1) % L]
        ln = (b - a) if b > a else (b - a + 1); mid = (a + ln / 2) % 1
        if all(dist_F(u * mid) > THR for u in Sp):
            comps.append((a, b, ln, mid))
    return comps


def automaton_empty_F(Sp, n):
    """Returns (counts_as_residual_row: bool, is_empty: bool) via Fractions."""
    comps = G_components_F(Sp, n)
    if not comps:
        return False, False
    wmax = ((n - 1) * max(Sp)) // n + 1
    allowed = set(range(1, min(wmax, 200) + 1))
    for (a, b, ln, mid) in comps:
        allowed = {w for w in allowed
                   if dist_F(n * w * mid) <= F(1, n) - F(n * w, 2) * ln}
        if not allowed:
            break
    return True, (not allowed)


# ---------- integer rewrite ----------
def lcm(a, b):
    return a // gcd(a, b) * b


def automaton_empty_int(Sp, n):
    """Exact integer version. Returns (counts_as_residual_row, is_empty)."""
    if not Sp:
        return False, False
    L = 1
    for u in Sp:
        L = lcm(L, u)
    D = n * L
    D2 = 2 * D
    twoL = 2 * L  # = 2D / n

    # breakpoints as integer numerators over D, with their owners attached
    pts = {}
    for u in Sp:
        c = L // u
        for k in range(u + 1):
            for eps in (1, -1):
                pint = ((k * n + eps) * c) % D
                pts.setdefault(pint, None)
    order = sorted(pts)
    Lo = len(order)

    comps = []  # store (lnint, midnum)
    for i in range(Lo):
        a = order[i]
        b = order[(i + 1) % Lo]
        lnint = (b - a) % D
        if lnint == 0:
            lnint = D  # full wrap (single distinct point) -- matches Fraction +1
        midnum = (2 * a + lnint) % D2
        ok = True
        for u in Sp:
            r = (u * midnum) % D2
            if min(r, D2 - r) <= twoL:   # dist(u*mid) > 1/n  fails
                ok = False
                break
        if ok:
            comps.append((lnint, midnum))

    if not comps:
        return False, False

    wmax = ((n - 1) * max(Sp)) // n + 1
    allowed = set(range(1, min(wmax, 200) + 1))
    for (lnint, midnum) in comps:
        new = set()
        for w in allowed:
            r2 = (n * w * midnum) % D2
            if min(r2, D2 - r2) <= twoL - n * w * lnint:
                new.add(w)
        allowed = new
        if not allowed:
            break
    return True, (not allowed)


# ---------- driver (RNG sequence identical to S1, so counts are comparable) ----------
def run(n, trials, rng, fn):
    m = n - 1; B = 2 * n + 8; emptyc = 0; tot = 0
    for _ in range(trials):
        others = set(rng.sample([x for x in range(1, B + 1) if x % n != 0], m - 1))
        ww = rng.randint(1, 3); v = n * ww
        if v in others:
            continue
        V = tuple(sorted(others | {v}))
        g = 0
        for x in V:
            g = gcd(g, x)
        V = tuple(sorted(x // g for x in V))
        if len(V) != m or not any(x % n == 0 for x in V):
            continue
        mults = [x for x in V if x % n == 0]; vv = mults[0]
        Sp = tuple(x for x in V if x != vv)
        counts, is_empty = fn(Sp, n)
        if not counts:
            continue
        tot += 1
        if is_empty:
            emptyc += 1
    return tot, emptyc


def main():
    print("Integer-rewrite bounded-CRT-automaton emptiness (large-owner residual => loose).")
    print("monad-compute-2026-06-03-S2. EMPTY/tot over random configs.")
    print("=" * 70)

    # ---- cross-check: integer vs Fraction must agree exactly ----
    print("CROSS-CHECK (int vs Fraction, same RNG, must match):")
    ok_all = True
    for n in [10, 12, 14]:
        rng_f = random.Random(1)
        rng_i = random.Random(1)
        tf = ef = ti = ei = 0
        for _ in range(n - 1):  # advance unused? no -- just run full
            pass
        tf, ef = run(n, 300, rng_f, automaton_empty_F)
        ti, ei = run(n, 300, rng_i, automaton_empty_int)
        match = (tf == ti and ef == ei)
        ok_all &= match
        print(f"   n={n:3d}: Fraction (tot={tf},EMPTY={ef})  int (tot={ti},EMPTY={ei})  "
              f"{'MATCH' if match else '*** MISMATCH ***'}")
    print(f"   cross-check: {'ALL MATCH' if ok_all else '*** FAILURE ***'}")
    print("=" * 70)

    # ---- extended integer run ----
    print("EXTENDED INTEGER RUN (1000 trials/n):")
    rng = random.Random(1)
    grand_tot = 0; grand_empty = 0
    for n in [10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30]:
        t0 = time.time()
        tot, emptyc = run(n, 1000, rng, automaton_empty_int)
        dt = time.time() - t0
        grand_tot += tot; grand_empty += emptyc
        flag = "ALL EMPTY" if emptyc == tot else f"*** {tot - emptyc} NON-EMPTY ***"
        print(f"   n={n:3d}: residual rows={tot:4d}; EMPTY={emptyc:4d}/{tot:<4d}  "
              f"{flag}   ({dt:.1f}s)")
    print("=" * 70)
    print(f"   TOTAL across n=10..30: {grand_empty}/{grand_tot} empty "
          f"({'NO COUNTEREXAMPLE' if grand_empty == grand_tot else '*** COUNTEREXAMPLE FOUND ***'})")
    print("   monad-compute-2026-06-03-S2")


if __name__ == '__main__':
    main()
