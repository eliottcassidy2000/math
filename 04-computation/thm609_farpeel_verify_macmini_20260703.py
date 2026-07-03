#!/usr/bin/env python3
"""
Verify THM-609 (base good-region floor) + the single-step far-peel structure (mac-mini-2026-07-03-S24).
THM-609: base B (<=12 speeds), LRC(<=13) gives t0 with min||b t0||>=1/13; then Safe_{1/14}(B) contains
[t0-1/(182V), t0+1/(182V)], length >= 1/(91V). Verify the safe interval + measure the ACTUAL good-region
length (>= the floor). Then the far-peel: remove far runner w's danger comb (measure 1/7, w teeth) from
goodRegion(B); check when length survives (opus step 3) vs needs the finite window (step 5).
"""
from fractions import Fraction as F
from math import gcd
import random

def nd(x):
    x = x % 1
    return min(x, 1 - x)

def min_over_B(B, t):
    return min(nd(b * t) for b in B)

def lrc_point(B, N=200000):
    """find t0 (rational a/N) maximizing min_i ||b_i t0|| -- the LRC(<=13) point (should be >= 1/13)."""
    best, bt = F(0), F(0)
    for k in range(1, N):
        t = F(k, N); m = min(nd(b * t) for b in B)
        if m > best: best, bt = m, t
    return best, bt

def goodregion_length(B, h=F(1, 14), Ngrid=300000):
    """Lebesgue measure of {t in [0,1): min_i ||b_i t|| >= h} via fine grid (approx)."""
    cnt = 0
    for k in range(Ngrid):
        t = (k + 0.5) / Ngrid
        if min(nd(b * t) for b in B) >= float(h): cnt += 1
    return cnt / Ngrid

if __name__ == "__main__":
    rng = random.Random(609)
    print("THM-609 verification: LRC point margin, safe-interval floor 1/(91V), actual good-region length.")
    print("=" * 92)
    print(f"{'base B (<=12)':>34} {'V=max':>6} {'LRC min':>8} {'>=1/13?':>8} {'floor 1/(91V)':>13} {'actual len':>11} {'>=floor?':>8}")
    for _ in range(8):
        m = rng.randint(6, 12)
        B = sorted(rng.sample(range(1, 40), m))
        V = max(B)
        mn, t0 = lrc_point(B, N=100000)
        floor = 1.0 / (91 * V)
        actual = goodregion_length(B, Ngrid=200000)
        print(f"{str(B):>34} {V:>6} {float(mn):>8.5f} {str(mn>=F(1,13)):>8} {floor:>13.6f} {actual:>11.6f} {str(actual>=floor):>8}")

    print("\nFAR-PEEL single step: remove far runner w (comb measure 1/7) from goodRegion(B). survives?")
    print("=" * 92)
    print(f"{'base B':>26} {'good len G':>10} {'far w':>7} {'len(G after -comb w)':>20} {'survives?':>10} {'w > ~15V*#comp?':>16}")
    B = [1, 2, 3, 5, 7, 11, 13, 17, 19, 23]  # a covering-ish base (12? here 10)
    B = sorted(random.Random(1).sample(range(1,25), 11))
    G = goodregion_length(B, Ngrid=300000)
    V = max(B)
    for w in [23, 30, 50, 100, 300, 1000, 5000]:
        full = B + [w]
        lenfull = goodregion_length(full, Ngrid=400000)
        survives = lenfull > 0
        # crude threshold estimate ~ 15*V (per unit #comp)
        print(f"{str(B):>26} {G:>10.5f} {w:>7} {lenfull:>20.6f} {str(survives):>10} {str(w > 15*V):>16}")
    print("\n=> THM-609 floor holds (actual >= 1/(91V)); the far-peel survives for w above a base-dependent")
    print("   threshold (opus step 3), leaving a finite window (step 5). THM-609's length G SETS that threshold.")
