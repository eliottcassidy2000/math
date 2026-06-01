#!/usr/bin/env python3
"""
lrc_all_observers_s521.py   claudebox-2026-06-01-S521

LRC from all n angles at once (reflection:
07-reflections/lrc-all-runners-are-observers-s521.md).

n runners, distinct speeds w_1..w_n. Full danger graph: edge between runners within
1/n. Runner k lonely <=> k isolated. LRC <=> every vertex isolated at some time.
Findings: (1) loneliness ROTATES (min #lonely=0; max=n only at the regular polygon);
(2) the n observer-frames D_k={w_i-w_k} have different difficulty -- collapsed
(repeated) frames are easy, distinct/AP-like frames are tight (HIDDEN TIGHTNESS).
"""
from fractions import Fraction as F
from math import gcd

def dist(x):
    x = x % 1; return min(x, 1 - x)
def lonely_k(w, k, t, n):
    return all(dist(F(w[i]-w[k])*t) >= F(1, n) for i in range(n) if i != k)
def cells_w(w, n):
    W = set([F(0)])
    for i in range(n):
        for j in range(n):
            d = abs(w[i]-w[j])
            if d:
                for k in range(d+1):
                    W.add((F(k, d)+F(1, n*d)) % 1); W.add((F(k, d)-F(1, n*d)) % 1)
    W = sorted(x for x in W if 0 <= x < 1); W2 = W + [F(1)]
    return W + [(a+b)/2 for a, b in zip(W, W2[1:])]
def lonely_measure(w, k, n):
    W = set([F(0)])
    for i in range(n):
        if i == k: continue
        d = abs(w[i]-w[k])
        for kk in range(d+1):
            W.add((F(kk, d)+F(1, n*d)) % 1); W.add((F(kk, d)-F(1, n*d)) % 1)
    W = sorted(x for x in W if 0 <= x < 1); W2 = W + [F(1)]
    return sum((b-a) for a, b in zip(W, W2[1:]) if lonely_k(w, k, (a+b)/2, n))

def main():
    print("Loneliness rotates: min/max # simultaneously lonely (max=n only at regular polygon):")
    for w in [(0,1,2,3),(0,1,2,4),(0,1,3,4,9),(0,2,3,5,7)]:
        n = len(w); ts = cells_w(list(w), n)
        nums = [sum(1 for k in range(n) if lonely_k(list(w), k, t, n)) for t in ts]
        each = all(any(lonely_k(list(w), k, t, n) for t in ts) for k in range(n))
        print(f"  {w}: min#lonely={min(nums)}, max#lonely={max(nums)}, each lonely sometime={each}")
    print("\nHidden tightness: per-frame lonely-measure & reduced difference-set (0 = tight observer):")
    for w in [(0,1,2,3),(0,2,3,5,7)]:
        n = len(w)
        rows = sorted((float(lonely_measure(list(w), k, n)), w[k],
                       tuple(sorted(set(abs(w[i]-w[k]) for i in range(n) if i != k)))) for k in range(n))
        for meas, sp, Dk in rows:
            tag = " <- TIGHT frame" if meas == 0 else ""
            print(f"  {w} runner@{sp}: measure={meas:.4f} D={Dk} ({len(Dk)} distinct){tag}")
        print()

if __name__ == "__main__":
    main()
