#!/usr/bin/env python3
"""
DISPROOF HUNT for LRC@14, given the crux "AP is the unique tight 13-set".
opus-2026-06-01-S553 (remote-control).

Two ways to break the program:
  (D1) a NON-AP tight 13-set (M=1/14): kills the "AP unique-tight => LRC" route;
  (D2) a config with M < 1/14: an actual COUNTEREXAMPLE to LRC(14) (historic).

We hunt both, going BEYOND the [1,21] census (speeds up to ~26 and via hill-climb
to ~60), targeting exactly the structures that produced sporadic tight sets at
small n:

  A. ALL 2^13 antipodal TRANSVERSALS mod 27 (speeds = chosen residues in 1..26),
     i.e. oracle-S553's reduced family -- tested at the real target n=14 (they
     only did n<=8).  Classify each; report any non-AP tight or any M<1/14.
  B. NON-transversals missing a NON-UNIT antipodal pair mod 27 (the n=8 mechanism;
     2n-1=27=3^3 has non-units 3,6,9,...).  Random + structured sample.
  C. Hill-climb MINIMISING the safe-measure at threshold 1/14 toward 0, over
     speeds <= 60, from many seeds (AP, transversals, n=8-style).  Any config
     reaching measure 0 that is not the AP is a D1/D2 hit.

Exact verdicts (Fraction); float pre-scan only for speed.
"""

from fractions import Fraction
from math import gcd
from itertools import product
import random

N = 14
M = 2 * N - 1          # 27
THR = Fraction(1, N)
THRF = 1.0 / N


def primitive(V):
    g = 0
    for v in V:
        g = gcd(g, v)
    return tuple(sorted(v // g for v in V)) if g else tuple(sorted(V))


def classify_exact(V):
    """'loose' (M>1/n) / 'tight' (M=1/n) / 'COUNTEREX' (M<1/n). Exact."""
    endpoints = set()
    for v in V:
        for k in range(v + 1):
            for s in (-1, 1):
                endpoints.add(Fraction(k * N + s, v * N) % 1)
    pts = sorted(endpoints)
    L = len(pts)
    for i in range(L):                      # positive-measure open interval?
        a, b = pts[i], pts[(i + 1) % L]
        length = (b - a) if b > a else (b - a + 1)
        mid = (a + length / 2) % 1
        if all(min((v * mid) % 1, 1 - (v * mid) % 1) > THR for v in V):
            return "loose"
    for t in pts:                           # measure-zero tight point?
        if all(min((v * t) % 1, 1 - (v * t) % 1) >= THR for v in V):
            return "tight"
    return "COUNTEREX"


def float_measure(V):
    """Fast float approx of safe-set measure at threshold 1/N (for screening)."""
    pts = []
    for v in V:
        for k in range(v + 1):
            base = k / v
            pts.append((base + THRF / v) % 1.0)
            pts.append((base - THRF / v) % 1.0)
    pts.sort()
    L = len(pts)
    meas = 0.0
    for i in range(L):
        a = pts[i]
        b = pts[i + 1] if i + 1 < L else pts[0] + 1.0
        mid = (a + (b - a) / 2.0) % 1.0
        ok = True
        for v in V:
            x = (v * mid) % 1.0
            if (x if x < 1 - x else 1 - x) < THRF - 1e-12:
                ok = False
                break
        if ok:
            meas += (b - a)
    return meas


AP = tuple(range(1, N))


# ----------------------------------------------------------------------
def partA_transversals():
    print("== A. ALL 2^13 antipodal transversals mod 27 (speeds = residues) ==")
    pairs = [(a, M - a) for a in range(1, N)]      # 13 pairs, a=1..13
    tight_nonAP = []
    counterex = []
    count = 0
    seen = set()
    for choice in product(*[(lo, hi) for (lo, hi) in pairs]):
        V = tuple(sorted(choice))
        if V in seen:
            continue
        seen.add(V)
        g = 0
        for v in V:
            g = gcd(g, v)
        if g != 1:
            continue
        count += 1
        if float_measure(V) > 1e-6:           # clearly loose, skip exact
            continue
        c = classify_exact(V)
        if c == "tight" and primitive(V) != AP:
            tight_nonAP.append(V)
        elif c == "COUNTEREX":
            counterex.append(V)
    print(f"   gcd-1 transversals tested: {count}")
    print(f"   NON-AP tight transversals (D1 hits): {len(tight_nonAP)}")
    for V in tight_nonAP[:30]:
        print(f"        D1: {V}")
    print(f"   COUNTEREXAMPLES M<1/14 (D2 hits): {len(counterex)}")
    for V in counterex[:30]:
        print(f"        D2: {V}")
    # sanity: AP itself is tight
    print(f"   sanity AP {AP}: {classify_exact(AP)}")
    print()
    return tight_nonAP, counterex


def partB_nonunit(samples=120000, seed=1):
    print("== B. Non-transversals missing a NON-UNIT antipodal pair mod 27 ==")
    rng = random.Random(seed)
    nonunit_pairs = [(a, M - a) for a in range(1, N) if gcd(a, M) != 1]  # 3,6,9,12
    print(f"   non-unit pairs mod 27: {nonunit_pairs}")
    tight_nonAP = []
    counterex = []
    tested = 0
    for _ in range(samples):
        miss = rng.choice(nonunit_pairs)
        allowed = [r for r in range(1, M) if r not in miss]     # residues != missed pair
        # pick 13 residues (with repetition across speeds via +27 lifts), speeds<=~80
        V = set()
        while len(V) < N - 1:               # n-1 = 13 speeds
            r = rng.choice(allowed)
            lift = r + M * rng.randint(0, 2)
            if lift > 0:
                V.add(lift)
        V = tuple(sorted(V))
        g = 0
        for v in V:
            g = gcd(g, v)
        if g != 1:
            continue
        tested += 1
        if float_measure(V) > 1e-6:
            continue
        c = classify_exact(V)
        if c == "tight" and primitive(V) != AP:
            tight_nonAP.append(V)
        elif c == "COUNTEREX":
            counterex.append(V)
    print(f"   tested (gcd-1, residues missing a non-unit pair): {tested}")
    print(f"   NON-AP tight (D1 hits): {len(tight_nonAP)}")
    for V in tight_nonAP[:30]:
        print(f"        D1: {V}  -> classify {classify_exact(V)}")
    print(f"   COUNTEREXAMPLES (D2 hits): {len(counterex)}")
    for V in counterex[:30]:
        print(f"        D2: {V}")
    print()
    return tight_nonAP, counterex


def partC_hillclimb(restarts=250, steps=250, cap=60, seed=3):
    print(f"== C. Hill-climb MINIMISING safe-measure (cap {cap}) toward tight ==")
    rng = random.Random(seed)
    hits = []
    best_overall = (1.0, None)
    # seeds: random, AP-neighbours, transversal-like
    SP = N - 1                               # 13 speeds
    for r in range(restarts):
        if r % 4 == 0:
            V = list(AP)
            for _ in range(rng.randint(1, 4)):
                V[rng.randrange(len(V))] += rng.choice([N, -N, 2 * N, -2 * N, 1, -1])
            V = sorted({abs(x) for x in V if x != 0})
            if len(V) != SP:
                V = rng.sample(range(1, cap + 1), SP)
        else:
            V = rng.sample(range(1, cap + 1), SP)
        V = list(V)
        cur = float_measure(tuple(V))
        for _ in range(steps):
            i = rng.randrange(len(V))
            old = V[i]
            cand = rng.randint(1, cap)
            if cand in V:
                continue
            V[i] = cand
            nm = float_measure(tuple(V))
            if nm <= cur:
                cur = nm
            else:
                V[i] = old
        Vt = primitive(tuple(V))
        if cur < best_overall[0]:
            best_overall = (cur, Vt)
        if cur < 1e-7:                         # looks tight/counterex -> confirm
            c = classify_exact(Vt)
            if c in ("tight", "COUNTEREX") and Vt != AP:
                hits.append((Vt, c))
    print(f"   best float-measure reached: {best_overall[0]:.3e} at {best_overall[1]}")
    print(f"   confirmed non-AP tight/counterex hits: {len(hits)}")
    for V, c in hits[:30]:
        print(f"        {c}: {V}")
    print()
    return hits


if __name__ == "__main__":
    a1, a2 = partA_transversals()
    b1, b2 = partB_nonunit()
    c = partC_hillclimb()
    print("================ VERDICT ================")
    total_D1 = len(a1) + len(b1) + len([1 for V, cc in c if cc == 'tight'])
    total_D2 = len(a2) + len(b2) + len([1 for V, cc in c if cc == 'COUNTEREX'])
    print(f"  D1 (non-AP tight 13-set) hits: {total_D1}")
    print(f"  D2 (M<1/14 counterexample) hits: {total_D2}")
    if total_D1 == 0 and total_D2 == 0:
        print("  => NO disproof found: AP remains the only tight 13-set across all")
        print("     tested transversals (speeds<=26), non-unit-pair non-transversals,")
        print("     and minimise-measure hill-climbs (speeds<=60). Crux survives.")
