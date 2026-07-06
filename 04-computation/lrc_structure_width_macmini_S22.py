#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S22 -- HYP-4502: THE METRIC HALF of opus-S113's STRUCTURE x
WIDTH.  opus: gap member => generalized AP (structure) + Farey wall q>=3k+2.
I quantify: a generalized AP's min M-RISE vs the window width, across k, to see
WHY the gap is nonempty at small k but empty at k=12.

For LRC(k+1): k speeds, floor 1/(k+1), second value 2/(2k+1), window
(1/(k+1), 2/(2k+1)) of width 1/((k+1)(2k+1)).  A generalized-AP DEFICIT family
(AP {1..k} + perturbations / lifts, or a 2D GAP) has M > 1/(k+1); does the rise
land in the window (gap member) or overshoot to >= 2/(2k+1)?
"""
import itertools, random
from fractions import Fraction as F
from math import gcd

def exact_M(W):
    dens = set()
    for v, w in itertools.combinations(W, 2):
        dens.add(v + w)
        if v != w:
            dens.add(abs(v - w))
    for v in W:
        dens.add(2 * v)
    best = F(0); seen = set()
    for s in dens:
        if s == 0:
            continue
        for j in range(1, s):
            t = F(j, s)
            if t in seen:
                continue
            seen.add(t)
            mv = min(abs(v * t - round(v * t)) for v in W)
            if mv > best:
                best = mv
    return best

def primitive(W):
    g = 0
    for v in W:
        g = gcd(g, v)
    return tuple(sorted(v // g for v in W)) if g > 1 else tuple(sorted(W))

def gen_ap_deficits(k, seed=0):
    """generalized-AP DEFICIT families with k elements: (a) single-lifts of the
    AP {1..k}; (b) 2-element block-lifts; (c) 2D GAPs (two-AP structure)."""
    rng = random.Random(seed)
    AP = list(range(1, k + 1))
    out = set()
    q1, q2 = k + 1, 2 * k + 1
    # (a) single lifts by +q1, +q2, +2q1, ...
    for i in range(k):
        for d in (q1, q2, 2 * q1, q1 + q2):
            W = list(AP); W[i] += d
            out.add(primitive(tuple(sorted(set(W)))))
    # (b) block-lifts: remove two, add two lifted (the 2/25-attainer species)
    for i, j in itertools.combinations(range(k), 2):
        for di, dj in [(q1, q1), (q2, q2), (q1, q2)]:
            W = list(AP); W[i] += di; W[j] += dj
            out.add(primitive(tuple(sorted(set(W)))))
    # (c) 2D generalized APs: {a + i*d1 + j*d2} (opus's n=7 shape: AP + boundary)
    for d1 in range(2, k):
        for d2 in range(1, 6):
            for L1 in range(2, k):
                base = [1 + i * d1 for i in range(L1)]
                extra = [1 + i * d1 + d2 for i in range(k - L1)]
                W = primitive(tuple(sorted(set(base + extra))))
                if len(W) == k:
                    out.add(W)
    # (d) random near-AP
    for _ in range(3000):
        W = list(AP)
        for _i in range(rng.randint(1, 3)):
            W[rng.randrange(k)] += rng.choice([q1, q2, 2*q1]) * rng.randint(1, 2)
        out.add(primitive(tuple(sorted(set(W)))))
    return [W for W in out if len(W) == k]

def main():
    print("=" * 82)
    print("STRUCTURE x WIDTH: generalized-AP M-rise vs window, across k")
    print("=" * 82)
    print(f"  {'k':>3} {'floor':>7} {'top 2/(2k+1)':>12} {'window w':>12} "
          f"{'min in-gap M':>14} {'min rise/w':>11} {'gap?':>6}")
    for k in range(6, 14):
        floor = F(1, k + 1); top = F(2, 2 * k + 1); w = top - floor
        fams = gen_ap_deficits(k, seed=k)
        in_gap = []
        min_rise = None
        for W in fams:
            M = exact_M(W)
            if M > floor:
                rise = M - floor
                if min_rise is None or rise < min_rise:
                    min_rise = rise
                if floor < M < top:
                    in_gap.append((M, W))
        ming = min(in_gap)[0] if in_gap else None
        rise_over_w = (min_rise / w) if min_rise else None
        print(f"  {k:>3} {'1/'+str(k+1):>7} {'2/'+str(2*k+1):>12} "
              f"{str(w):>12} {(str(ming) if ming else '-- (empty)'):>14} "
              f"{(f'{float(rise_over_w):.3f}' if rise_over_w else '--'):>11} "
              f"{('YES '+str(len(in_gap)) if in_gap else 'empty'):>6}")
    print()
    print("  READING: 'min rise/w' = (min M-rise)/(window width).  If < 1, the")
    print("  smallest deficit-rise fits INSIDE the window => gap NONEMPTY.  If the")
    print("  min-rise reaches the TOP (2/(2k+1)), rise/w >= 1 => gap EMPTY.  The")
    print("  crossover k is where the structured families stop fitting = opus's")
    print("  STRUCTURE x WIDTH closure; k=12 (window 1/325) is the target.")

if __name__ == "__main__":
    main()
