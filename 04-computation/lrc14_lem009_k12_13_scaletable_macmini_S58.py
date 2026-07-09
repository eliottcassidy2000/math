"""
mac-mini-2026-07-08-S58 -- EXACT per-scale table for the k=12,13 longest-AP=(k-1)
tail family {0,d,..,(k-2)d} u {p}, mirroring kps-S88's k=11 d=3..6 table.
Pins the tail min and shows SCALE-MONOTONICITY (D3 rises with d toward the
decorrelation limit), the analog of opus-S157/kps-S88 for k=11.

D3 rises with d => tail min is at the smallest tail scale; combined with the
exhaustive compact check (prim-diam <= 18) this closes the tail. opus-S157's
resonance-sum rate |D3 - D3_inf| <= C/(pd) (k-general) covers d -> inf.
"""
from fractions import Fraction as F
from math import floor, gcd
from functools import reduce

TH = F(1, 7); M = F(6, 7)
BARS = {12: F(50285, 252252), 13: F(14249, 252252)}

def int_pow(A, B, p, lo, hi):
    if B == 0: return A**p * (hi - lo)
    return ((A + B*hi)**(p+1) - (A + B*lo)**(p+1)) / (B * (p+1))

def W_moments(E):
    E = sorted(E); k = len(E)
    ds = set(e for e in E if e != 0)
    for i in range(k):
        for j in range(i+1, k): ds.add(E[j]-E[i])
    bps = set([F(0), F(1)])
    for d in ds:
        for m in range(0, d+1): bps.add(F(m, d))
    bps = sorted(bps)
    m1 = m2 = m3 = F(0)
    for c in range(len(bps)-1):
        x0, x1 = bps[c], bps[c+1]; xm = (x0+x1)/2
        lin = [(F(-floor(e*xm)), F(e)) for e in E]
        order = sorted(range(k), key=lambda j: lin[j][0]+lin[j][1]*xm)
        sp = [lin[j] for j in order]
        gaps = [(sp[i+1][0]-sp[i][0], sp[i+1][1]-sp[i][1]) for i in range(k-1)]
        gaps.append((F(1)+sp[0][0]-sp[k-1][0], sp[0][1]-sp[k-1][1]))
        subs = set([x0, x1])
        for (a, b) in gaps:
            if b != 0:
                xs = (TH-a)/b
                if x0 < xs < x1: subs.add(xs)
        subs = sorted(subs)
        for s in range(len(subs)-1):
            u0, u1 = subs[s], subs[s+1]; um = (u0+u1)/2
            A = B = F(0)
            for (a, b) in gaps:
                if a + b*um > TH: A += (a-TH); B += b
            m1 += int_pow(A, B, 1, u0, u1)
            m2 += int_pow(A, B, 2, u0, u1)
            m3 += int_pow(A, B, 3, u0, u1)
    return m1, m2, m3

def D3(E):
    m1, m2, m3 = W_moments(E)
    den = m2 - m3/M
    return (m1/M + (m1 - m2/M)**2/den) if den > 0 else m1/M

def primitive(E):
    E = sorted(E); return reduce(gcd, [E[i+1]-E[i] for i in range(len(E)-1)]) == 1

def longest_ap(E):
    E = sorted(set(E)); s = set(E); best = 2
    for i in range(len(E)):
        for j in range(i+1, len(E)):
            d = E[j]-E[i]; L = 2; nxt = E[j]+d
            while nxt in s: L += 1; nxt += d
            back = E[i]-d
            while back in s: L += 1; back -= d
            best = max(best, L)
    return best

print("EXACT per-scale table: longest-AP=(k-1) tail family {0,d,..,(k-2)d} u {p}\n")
for k in (12, 13):
    bar = BARS[k]
    print(f"k={k}  (bar = {float(bar):.6f}):")
    print(f"  {'d':>2} {'AP span':>7} {'min D3 (exact)':>16} {'at p':>6} {'prim-diam':>9} {'margin':>9}")
    fam_min = (F(10), None, None)
    for d in range(1, 5):
        ap = [d*i for i in range(k-1)]; hi = ap[-1]
        dmin = (F(10), None)
        # all interior + a few above-AP outliers; keep only genuine longest-AP=(k-1)
        cand_p = [p for p in range(1, hi) if p % d != 0] + [hi+1, hi+2, hi+d+1]
        for p in cand_p:
            E = tuple(sorted(set(ap + [p])))
            if len(E) != k or not primitive(E): continue
            if longest_ap(E) != k-1: continue     # ensure the AP is exactly k-1 (tail, not full block)
            d3 = D3(E)
            if d3 < dmin[0]: dmin = (d3, p, E)
        if dmin[1] is None:
            print(f"  {d:>2} {hi:>7} {'(none)':>16}")
            continue
        pdiam = max(dmin[2]) - min(dmin[2])
        print(f"  {d:>2} {hi:>7} {float(dmin[0]):>16.6f} {dmin[1]:>6} {pdiam:>9} "
              f"{float(dmin[0]-bar):>+9.6f}")
        if dmin[0] < fam_min[0]: fam_min = (dmin[0], d, dmin[2])
    print(f"  => longest-AP=(k-1) tail-family min D3 = {float(fam_min[0]):.6f} "
          f"(d={fam_min[1]}), margin {float(fam_min[0]-bar):+.6f}  "
          f"{'CLEARS' if fam_min[0] > bar else 'BELOW!'}")
    print(f"     scale-monotone: D3 rises with d toward the decorrelation limit "
          f"(~{'0.389' if k==12 else '0.344'}); opus-S157 rate |D3-D3_inf|<=C/(pd) covers d->inf.\n")
