#!/usr/bin/env python3
"""
lrc_tree_entropy_s543.py    oracle-2026-06-01-S543

Attack the ENTROPY of loneliness on the p-adic (base-p) tree.

The base-p tree = base-p expansions of t in [0,1): a depth-d node is an interval
[j/p^d, (j+1)/p^d). The descent (child = append a digit) is the x->p expanding
map; its entropy is log p. The CENTER = shift (S541, t-rotation / +1 odometer on
Z_p) has ZERO entropy (equicontinuous), so the LRC-relevant entropy lives TRANSVERSE
to the center, in this descent direction.

Two box-entropies of the SAFE (lonely) set S = { t : min_i ||v_i t|| >= 1/n }:
  N_touch(d) = #depth-d nodes that CONTAIN a safe point (intersect S)
  N_full(d)  = #depth-d nodes ENTIRELY safe (subset of S)
  h_touch = lim (1/d) log_p N_touch ,  h_full = lim (1/d) log_p N_full

Order-parameter picture to test:
  generic (|S|>0):      h_touch = h_full = 1   (lonely with positive measure)
  tight  (|S| measure 0, wall-only AP):  h_touch = 0, N_full = 0  (critical line)
  LRC-violating (no safe point): N_touch = 0    (would be 'h = -inf')
So LRC  <=>  N_touch(d) > 0 for some/all d  <=>  h_touch >= 0; the tight AP is the
h_touch: 1 -> 0 critical drop.
"""
from fractions import Fraction
from math import log
from functools import reduce
from math import gcd

def frac(x): return x - (x.numerator // x.denominator)
def d0(x):
    f = frac(x); return min(f, 1 - f)

def safe_intervals(speeds, n):
    """exact safe set as (list of [a,b] closed safe intervals, list of isolated safe points)."""
    thr = Fraction(1, n)
    W = set([Fraction(0), Fraction(1)])
    for s in speeds[1:]:
        if s == 0: continue
        for m in range(0, abs(s) + 1):
            for sg in (1, -1):
                t = Fraction(n * m + sg, n * abs(s))
                if 0 <= t <= 1: W.add(t)
    Wl = sorted(W)
    intervals = []
    for a, b in zip(Wl, Wl[1:]):
        mid = (a + b) / 2
        if min(d0(Fraction(s) * mid) for s in speeds[1:]) >= thr:
            intervals.append((a, b))
    # merge touching intervals
    merged = []
    for a, b in intervals:
        if merged and merged[-1][1] == a: merged[-1] = (merged[-1][0], b)
        else: merged.append((a, b))
    # isolated safe wall points (tight boundary loneliness) not covered by an interval
    pts = []
    covered = lambda t: any(a <= t <= b for a, b in merged)
    for t in Wl:
        if not covered(t) and min(d0(Fraction(s) * t) for s in speeds[1:]) >= thr:
            pts.append(t)
    return merged, pts

def box_counts(merged, pts, p, D):
    Ntouch, Nfull = [], []
    for d in range(1, D + 1):
        cell = Fraction(1, p ** d)
        nt = nf = 0
        for j in range(p ** d):
            lo, hi = j * cell, (j + 1) * cell
            full = any(a <= lo and hi <= b for a, b in merged)
            touch = full or any(a < hi and lo < b for a, b in merged) or any(lo <= t < hi for t in pts)
            if full: nf += 1
            if touch: nt += 1
        Ntouch.append(nt); Nfull.append(nf)
    return Ntouch, Nfull

def slope_log(counts, p):
    """h = (1/d) log_p N at the deepest d with N>0 (box-entropy estimate)."""
    for d in range(len(counts), 0, -1):
        if counts[d - 1] > 0:
            return log(counts[d - 1]) / (d * log(p))
    return float("-inf")

def primitive(s): return reduce(gcd, [x for x in s if x]) == 1

def main():
    print("LRC loneliness ENTROPY on the base-p tree (oracle-S543)\n")
    cases = [
        ("AP n=5 (tight)", (0,1,2,3,4), 5),
        ("generic n=5", (0,1,3,5,7), 5),
        ("AP n=4 (tight)", (0,1,2,3), 4),
        ("generic n=4", (0,1,3,5), 4),
        ("AP n=6 (tight)", (0,1,2,3,4,5), 6),
        ("generic n=6", (0,1,2,4,7,8), 6),
    ]
    for p in (2, 3):
        print(f"--- base p = {p} (descent entropy log p; depth D=7) ---")
        print(f"   {'system':<18} {'|S|':>10} {'h_full':>7} {'h_touch':>8}  N_touch(d=1..7)")
        for name, sp, n in cases:
            if not primitive(sp): continue
            merged, pts = safe_intervals(sp, n)
            meas = sum(b - a for a, b in merged)
            nt, nf = box_counts(merged, pts, p, 7)
            hf, ht = slope_log(nf, p), slope_log(nt, p)
            print(f"   {name:<18} {str(meas):>10} {hf:>7.3f} {ht:>8.3f}  {nt}")
        print()
    print("READING: generic systems have h_full = h_touch = 1 (the lonely set has full")
    print("box-dimension 1 in every base p -> loneliness with positive measure). The TIGHT")
    print("AP has |S|=0, N_full=0 (h_full=-inf) and N_touch BOUNDED (h_touch=0): the safe set")
    print("is a finite set of wall points -> the critical/zero-entropy line. So the loneliness")
    print("box-entropy h_touch is the order parameter: LRC <=> h_touch >= 0 (a safe point")
    print("exists); the extremal AP sits exactly at the h_touch: 1->0 drop. The center=shift")
    print("(odometer, zero entropy) is invisible here; the entropy is the transverse x->p")
    print("descent, independent of p (dim 1 generic / dim 0 tight in every base).")

if __name__ == "__main__":
    main()
