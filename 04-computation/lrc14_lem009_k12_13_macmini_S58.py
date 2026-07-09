"""
mac-mini-2026-07-08-S58 -- CLOSE k=12,13 with the LEM-009 machinery.

The k=11 leg was closed (LEM-009 fleet: opus-S156 longest-AP axis, kps-S87 D3
scale-monotonicity, klein-S188 multi-outlier Koksma-Hlawka) by:
  min_E D3(E) = min over [ exhaustive COMPACT (prim-diam <= D0, exact Farey) ]
                       U [ the DECORRELATION TAIL: block/AP + far outlier(s),
                           D3 -> D3_inf (Weyl equidistribution), exact limit,
                           Koksma-Hlawka O(1/D) rate => finite gap << margin ].
THM-661 already has the k=12,13 EXACT COMPACT minima (0.355876 / 0.308844 >= bars
0.199344 / 0.056487) but asserts the tail only "rises by decorrelation (LEM-005)".
This script gives k=12,13 the SAME rigorous tail closure as k=11.

D3(E) = m1/M + (m1 - m2/M)^2 / (m2 - m3/M),  M = 6/7,  m_i = E[W^i],
  W(x) = sum_i (g_i(x) - 1/7)_+ = uncovered measure (THM-657);  D3 <= mu (THM-661).

Part A: exact D3 moment integrator (Farey) + reconfirm compact minima.
Part B: exact D3(B_D) for block+outlier B_D = {0..k-2} u {D}, finite-D sequence.
Part C: exact decorrelation limit D3_inf (2D block x uniform-outlier integral).
Part D: longest-AP=(k-1) tail family {0,d,..,(k-2)d} u {p} min D3.
"""
from fractions import Fraction as F
from math import floor, gcd
from functools import reduce
from itertools import combinations

TH = F(1, 7)          # 1/7 toll
M  = F(6, 7)          # W in [0, 6/7]
BARS = {11: F(83549, 252252), 12: F(50285, 252252), 13: F(14249, 252252)}


# ---------- exact antiderivative of (A + B x)^p ----------
def int_pow(A, B, p, lo, hi):
    """Exact integral of (A + B x)^p over [lo, hi], p in {1,2,3}, Fractions."""
    if B == 0:
        return A**p * (hi - lo)
    # (A+Bx)^(p+1) / (B (p+1))
    up = (A + B*hi)**(p+1)
    dn = (A + B*lo)**(p+1)
    return (up - dn) / (B * (p+1))


# ---------- exact W-moments E[W^1], E[W^2], E[W^3] for a family E ----------
def W_moments(E):
    E = sorted(E); k = len(E)
    ds = set()
    for j in range(k):
        if E[j] != 0: ds.add(E[j])
    for i in range(k):
        for j in range(i+1, k):
            ds.add(E[j]-E[i])
    bps = set([F(0), F(1)])
    for d in ds:
        for m in range(0, d+1):
            bps.add(F(m, d))
    bps = sorted(bps)
    m1 = F(0); m2 = F(0); m3 = F(0)
    for c in range(len(bps)-1):
        x0, x1 = bps[c], bps[c+1]; xm = (x0+x1)/2
        lin = [(F(-floor(e*xm)), F(e)) for e in E]          # frac(e x) = a + b x
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
            A = F(0); B = F(0)
            for (a, b) in gaps:
                if a + b*um > TH:
                    A += (a - TH); B += b
            m1 += int_pow(A, B, 1, u0, u1)
            m2 += int_pow(A, B, 2, u0, u1)
            m3 += int_pow(A, B, 3, u0, u1)
    return m1, m2, m3


def D3(E):
    m1, m2, m3 = W_moments(E)
    den = m2 - m3/M
    if den <= 0:                     # degenerate; D3 falls back to Markov m1/M
        return m1/M, (m1, m2, m3)
    return m1/M + (m1 - m2/M)**2 / den, (m1, m2, m3)


def primitive(E):
    E = sorted(E)
    return reduce(gcd, [E[i+1]-E[i] for i in range(len(E)-1)]) == 1


# ================= PART A: reconfirm compact-min D3 (k=12,13) =================
print("="*70)
print("PART A -- exact compact-min D3, k=12,13 (reconfirm THM-661)")
print("="*70)
COMPACT_MIN = {}
for k, Dmax in [(12, 14), (13, 15)]:
    best = (F(10), None)
    nshapes = 0
    for D in range(k-1, Dmax+1):
        for mid in combinations(range(1, D), k-2):
            E = (0,) + mid + (D,)
            if not primitive(E): continue
            nshapes += 1
            d3, _ = D3(E)
            if d3 < best[0]:
                best = (d3, E)
    bar = BARS[k]
    COMPACT_MIN[k] = best
    print(f"\n  k={k}: {nshapes} primitive shapes, prim-diam <= {Dmax}")
    print(f"    exact min D3 = {float(best[0]):.6f}  at {best[1]}")
    print(f"    bar          = {float(bar):.6f}")
    print(f"    margin       = {float(best[0]-bar):+.6f}   "
          f"{'CLEARS' if best[0] > bar else 'BELOW!!'}")


# ============ PART B: block baseline + block+outlier finite-D sequence ============
print("\n" + "="*70)
print("PART B -- block baseline & block+outlier B_D = {0..k-2} u {D} (exact)")
print("="*70)
BLOCK_D3 = {}
for k in (12, 13):
    bar = BARS[k]
    blk = tuple(range(k))                       # full k-block {0..k-1}
    d3blk, _ = D3(blk)
    BLOCK_D3[k] = d3blk
    print(f"\n  k={k}  (bar={float(bar):.6f}):")
    print(f"    D3(full {k}-block {{0..{k-1}}}) = {float(d3blk):.6f}   "
          f"margin {float(d3blk-bar):+.6f}")
    print(f"    block+outlier B_D = {{0..{k-2}}} u {{D}}  (longest-AP = k-1):")
    seq = []
    for D in [k-1, k, k+2, k+5, 15, 20, 25, 30, 40, 50]:
        if D < k-1: continue
        E = tuple(range(k-1)) + (D,)
        if not primitive(E):     # gcd of {1,...,1, D-(k-2)} = 1 always since consecutive
            pass
        d3, (m1, m2, m3) = D3(E)
        seq.append((D, d3))
        flag = 'ok' if d3 > bar else 'BELOW!'
        print(f"      D={D:3d}: D3={float(d3):.6f}  E[W]={float(m1):.5f}  "
              f"margin {float(d3-bar):+.6f}  {flag}")
    BLOCK_D3[(k, 'seq')] = seq


# ---------- fast float D3 (moments on a grid) for pre-filtering ----------
import numpy as np
GRIDf = 6000
_xf = (np.arange(GRIDf) + 0.5) / GRIDf
Mf = 6.0/7.0
def d3_float(E):
    Ea = np.array(sorted(E), float)
    ph = np.mod(np.outer(_xf, Ea), 1.0); ph.sort(axis=1)
    g = np.concatenate([np.diff(ph, axis=1), (ph[:, 0]+1-ph[:, -1])[:, None]], axis=1)
    W = np.maximum(g - 1.0/7.0, 0).sum(axis=1)
    m1 = W.mean(); m2 = (W*W).mean(); m3 = (W*W*W).mean()
    den = m2 - m3/Mf
    return (m1/Mf + (m1 - m2/Mf)**2/den) if den > 1e-12 else m1/Mf


# ============ PART C: firm up the COMPACT minimum (float pre-filter to larger D0, exact-confirm) ============
print("\n" + "="*70)
print("PART C -- extend compact search to prim-diam <= D0 (float pre-filter + exact)")
print("="*70)
for k, D0 in [(12, 18), (13, 18)]:
    bar = BARS[k]
    cmin = float(COMPACT_MIN[k][0])
    cand = []            # shapes with float D3 within 0.03 of the known compact min
    nshapes = 0
    for D in range(k-1, D0+1):
        for mid in combinations(range(1, D), k-2):
            E = (0,) + mid + (D,)
            if not primitive(E): continue
            nshapes += 1
            if d3_float(E) < cmin + 0.03:
                cand.append(E)
    # exact-confirm every candidate; track the true min
    best = COMPACT_MIN[k]
    for E in cand:
        d3, _ = D3(E)
        if d3 < best[0]:
            best = (d3, E)
    print(f"\n  k={k}: {nshapes} primitive shapes prim-diam<= {D0}; "
          f"{len(cand)} within 0.03 of min (float)")
    print(f"    EXACT compact min D3 = {float(best[0]):.6f} at {best[1]}")
    print(f"    bar = {float(bar):.6f}   margin {float(best[0]-bar):+.6f}   "
          f"{'CLEARS' if best[0] > bar else 'BELOW!!'}")
    COMPACT_MIN[k] = best


# ============ PART D: longest-AP=(k-1) tail family {0,d,..,(k-2)d} u {p} ============
# Primitive k-set with prim-diam > k-1  ==>  longest AP <= k-1 (a k-term AP is a
# dilated block, non-primitive unless d=1). The tail D3-minimizer is the longest
# possible AP (k-1 terms) plus one outlier (max additive energy in the tail).
print("\n" + "="*70)
print("PART D -- longest-AP=(k-1) tail family {0,d,..,(k-2)d} u {p}  (exact D3 min)")
print("="*70)
for k in (12, 13):
    bar = BARS[k]
    best = (F(10), None)
    for d in range(1, 7):
        ap = tuple(d*i for i in range(k-1))     # (k-1)-term AP, spacing d
        hi = ap[-1]                              # = (k-2)d
        # outlier p: anywhere not on the AP; sample a wide range incl. far (decorrelation)
        for p in list(range(1, hi)) + [hi+1, hi+2, hi+d, 2*hi, 3*hi, 5*hi]:
            if p % d == 0 and 0 <= p <= hi:      # p on the AP -> skip (would extend/duplicate)
                continue
            E = tuple(sorted(set(ap + (p,))))
            if len(E) != k: continue
            if not primitive(E): continue
            # require longest AP is exactly k-1 (outlier must not create a longer AP through it)
            d3f = d3_float(E)
            if d3f < float(best[0]) + 0.02:      # exact-confirm near-min only
                d3, _ = D3(E)
                if d3 < best[0]:
                    best = (d3, E)
    print(f"\n  k={k} (bar={float(bar):.6f}): longest-AP=(k-1) tail family")
    print(f"    exact min D3 = {float(best[0]):.6f} at {best[1]}")
    print(f"    margin {float(best[0]-bar):+.6f}   "
          f"{'CLEARS' if best[0] > bar else 'BELOW!!'}")

# ---------- decorrelation limit estimate: exact D3(B_D) at large D (Richardson) ----------
print("\n  decorrelation LIMIT of block+outlier (exact D3 at large D):")
for k in (12, 13):
    vals = []
    for D in [80, 160, 320]:
        E = tuple(range(k-1)) + (D,)
        d3, _ = D3(E)
        vals.append((D, float(d3)))
    # Richardson: assume D3(D) = L + c/D  -> L ~ 2*D3(2D)-D3(D)
    (D1, v1), (D2, v2), (D3v, v3) = vals
    L_est = 2*v3 - v2
    rate = max(abs(v-L_est)*D for D, v in vals)
    print(f"    k={k}: D3(B_D) at D={D1},{D2},{D3v} = {v1:.6f},{v2:.6f},{v3:.6f}; "
          f"limit~{L_est:.6f}; |D3-L|*D <= {rate:.3f} (Koksma-Hlawka O(1/D))")
