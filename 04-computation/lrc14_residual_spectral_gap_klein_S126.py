#!/usr/bin/env python3
"""
klein-2026-07-04-S126 (HYP-4080) - THE RESIDUAL SLIVER: is there a SPECTRAL GAP above 1/12?

The last open sliver of the m=2,f=2 confinement (covering-min rigidity) is the small-tightener x
near-AP corner: M(2U u {w1,w2}) >= 1/12 for LOOSE 11-runner even parts U (M(U)>1/12) with small
tighteners.  opus's Lemma 3 (klein S125 Lean) closes LARGE tighteners; the residual is where
M(U)-1/12 is tiny (Lemma 3's threshold u_max/(6(M(U)-1/12)) blows up).

KEY QUESTION (the crux "bound v_max(U)"): is
    delta0 := inf { M(U) - 1/12 : U an 11-runner set, M(U) > 1/12 }  >  0  ?
IF delta0 > 0 (a Lagrange-style GAP above the tight value 1/12), then loose U have M(U)-1/12 >=
delta0, so Lemma 3 fires for all tighteners > u_max/(6 delta0), bounding the residual's tighteners
relative to u_max -- a real closure lever.  IF M(U) accumulates at 1/12 (NO gap), the residual is
genuinely LRC-hard (no spectral-gap shortcut).

This script tests it exactly (Fractions): the M(U) spectrum just above 1/12 for 11-runner sets
(perturbations of the tight AP {1..11}), AND large-dilation perturbations (does M(U)-1/12 -> 0 as
the dilation grows?).  Honest either way.
"""
from fractions import Fraction as F
from math import gcd
import itertools

def cdist_q(a, q):
    r = a % q
    return min(r, q - r)

def Mval(S, Qcap):
    """exact M(S) = max_t min_{v in S} ||v t||, scanned over Farey denominators Q<=Qcap."""
    best = F(0)
    for Q in range(2, Qcap + 1):
        for a in range(1, Q // 2 + 1):
            if gcd(a, Q) != 1: continue
            m = min(F(cdist_q(v * a, Q), Q) for v in S)
            if m > best: best = m
    return best

AP = list(range(1, 12))   # {1,...,11}, the tight 11-runner set
print("="*78)
print("TIGHT ANCHOR: M({1..11}) should = 1/12 (LRC(12) extremal)")
print("="*78)
mAP = Mval(AP, 60)
print(f"  M({{1..11}}) = {mAP} (~{float(mAP):.6f});  1/12 = {float(F(1,12)):.6f};  tight? {mAP==F(1,12)}")

print()
print("="*78)
print("TEST 1: M(U) spectrum just above 1/12 -- perturb {1..11} (replace one element by m<=24).")
print("  Collect M(U) values in (1/12, 1/12 + 0.02]; is there a GAP just above 1/12?")
print("="*78)
vals = {}
for drop in range(1, 12):
    for m in range(1, 25):
        U = sorted(set(AP) - {drop} | {m})
        if len(U) != 11: continue
        g = 0
        for v in U: g = gcd(g, v)
        if g != 1: continue
        M = Mval(U, 60)
        if F(1,12) < M <= F(1,12) + F(1,50):
            vals[M] = vals.get(M, 0) + 1
spec = sorted(vals.keys())
print(f"  distinct M(U) values in (1/12, 1/12+0.02]: {len(spec)}")
for M in spec[:12]:
    print(f"    M={M} (~{float(M):.6f})  gap above 1/12 = {float(M-F(1,12)):.6f}  (#sets {vals[M]})")
if spec:
    print(f"  SMALLEST gap above 1/12 in this family: {float(spec[0]-F(1,12)):.6f} at M={spec[0]}")

print()
print("="*78)
print("TEST 2 (the decisive one): does M(U)-1/12 -> 0 for LARGE-DILATION perturbations?")
print("  U_c = c*{1..11} with ONE speed shifted: {c,2c,...,10c, 11c + d}. Tight AP dilated is")
print("  exactly 1/12; a small shift d makes it loose. Track M(U_c)-1/12 as c grows.")
print("  If -> 0, NO spectral gap => residual is LRC-hard (accumulates at 1/12).")
print("="*78)
print(f"{'c':>4} {'d':>3} {'U_c max':>8} {'M(U_c)':>16} {'M-1/12':>12}")
for d in [1, -1, 2]:
    for c in [1, 2, 3, 5, 8]:
        U = [c*k for k in range(1, 11)] + [11*c + d]
        U = sorted(set(U))
        if len(U) != 11: continue
        g = 0
        for v in U: g = gcd(g, v)
        # (allow imprimitive here; we only track M which is scale-sensitive via the shift)
        Qcap = min(3*(max(U)+1), 700)
        M = Mval(U, Qcap)
        print(f"{c:>4} {d:>3} {max(U):>8} {str(M):>16} {float(M-F(1,12)):>12.6f}")
print("  (If M-1/12 shrinks as c grows, the tight value 1/12 is an ACCUMULATION point => no gap.)")

print()
print("="*78)
print("READING: a positive smallest-gap that does NOT shrink with dilation => spectral gap =>")
print("residual bounded/finite (closure lever). A shrinking gap => no shortcut, LRC-hard sliver.")
print("="*78)
print("DONE")
