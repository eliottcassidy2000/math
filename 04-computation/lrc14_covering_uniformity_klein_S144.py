#!/usr/bin/env python3
"""
klein-2026-07-06-S144 - HEIGHT-UNIFORMITY OF THE (C) COVERING SYSTEM (the last open node).

Proof map: (C) = "AP {1..12} is the unique primitive 12-family (up to dilation) with M<2/25".
Reduced (kps S43/44, opus S124-126) to a FINITE COVERING SYSTEM: every NON-AP primitive 12-family
is "cleared" at some modulus q<=39, i.e. M_q(W) := max_{c coprime q} min_i cdist(v_i c, q)/q >= 2/25
(a rational_point_margin certificate => M(W) >= M_q(W) >= 2/25, LOOSE). The OPEN critical-path node
(proof map): the covering is UNIFORM OVER ALL HEIGHTS -- min-clearing-q does NOT grow with entry size.

This script is the decisive VALIDITY CHECK + mechanism probe:
  (1) For broad + HIGH-HEIGHT + ADVERSARIAL non-AP families, compute min-clearing-q (smallest q<=60
      with M_q >= 2/25). Find any needing q>39, or uncoverable to q<=60 (a covering GAP => kps refuted).
  (2) The critical adversarial class: CRT-LIFTS W == (dilated AP) mod L for L=lcm(subset of small q)
      -- these agree with the AP at many moduli, the natural way to try to obstruct many q at once.
  (3) mod-25-saturated blockers (case 2, the hard residual): confirm they clear at some q<=39.
Report the histogram of min-clearing-q and the WORST families. If max clearing-q stays bounded
(<= ~14 per kps S44) and no gaps, the covering is height-uniform (validated).
Exact (Fractions / integer cdist).
"""
from fractions import Fraction as F
from math import gcd, lcm
from random import seed, randint, sample

TWO25 = F(2, 25)

def cd(a, q):
    r = a % q
    return min(r, q - r)

def Mq(W, q):
    """q-modular covering-min = max_{c in [1,q/2], gcd(c,q)=1} min_i cdist(v_i c, q) / q."""
    best = 0
    for c in range(1, q // 2 + 1):
        if gcd(c, q) != 1:
            continue
        m = min(cd(v * c, q) for v in W)
        if m > best:
            best = m
    return F(best, q)

def min_clearing_q(W, qmax=60):
    """smallest q in [2,qmax] with M_q(W) >= 2/25, or None."""
    for q in range(2, qmax + 1):
        if Mq(W, q) >= TWO25:
            return q
    return None

def tgcd(W):
    g = 0
    for v in W: g = gcd(g, v)
    return g

def is_prim(W):
    return tgcd(W) == 1 and len(set(W)) == 12 and all(v > 0 for v in W)

AP = list(range(1, 13))
seed(20260706)

print(f"2/25 = {float(TWO25):.5f}. Clearing-q = smallest q<=60 with M_q(W) >= 2/25.")
print(f"sanity: AP {{1..12}} min-clearing-q = {min_clearing_q(AP)} (should be None: AP is NOT cleared -- it IS the exception)")
print("=" * 84)

def survey(label, families):
    hist = {}
    gaps = []           # families with no clearing q <= 39
    big = []            # families needing q in (14, 39]
    n = 0
    for W in families:
        W = sorted(W)
        if not is_prim(W) or W == AP:
            continue
        n += 1
        q = min_clearing_q(W, 60)
        if q is None or q > 39:
            gaps.append((q, W))
        else:
            hist[q] = hist.get(q, 0) + 1
            if q > 14:
                big.append((q, W))
    print(f"[{label}] {n} non-AP primitive families")
    print(f"   min-clearing-q histogram: {dict(sorted(hist.items()))}")
    print(f"   max clearing-q (<=39): {max(hist) if hist else 'n/a'}")
    if big:
        print(f"   needing q in (14,39]: {len(big)}  e.g. {big[:3]}")
    if gaps:
        print(f"   *** {len(gaps)} COVERING GAPS (no clear q<=39): {gaps[:5]}")
    else:
        print(f"   NO covering gaps: every non-AP family cleared at some q<=39. [uniformity holds here]")
    return gaps

# 1. broad random, moderate + HIGH height
survey("random height<=60", [sorted(sample(range(1, 61), 12)) for _ in range(4000)])
survey("random height<=5000", [sorted(sample(range(1, 5001), 12)) for _ in range(3000)])
survey("random height<=1e6", [sorted(sample(range(1, 1000001), 12)) for _ in range(2000)])

# 2. CRT-LIFTS: W == c*AP mod L, but non-AP as integers (the adversarial obstruct-many-q class)
def crt_lifts(Lmods, ntry, hmax):
    fams = []
    L = 1
    for q in Lmods: L = lcm(L, q)
    for _ in range(ntry):
        c = 1
        # W_i = c*i + L*k_i, random small k_i, keep primitive & distinct
        ks = [randint(0, hmax) for _ in range(12)]
        W = sorted(set((c * i + L * k) for i, k in zip(range(1, 13), ks)))
        if len(W) == 12 and tgcd(W) == 1 and W != AP:
            fams.append(W)
    return fams
survey("CRT-lift ==AP mod lcm(2..12)", crt_lifts(range(2,13), 3000, 30))
survey("CRT-lift ==AP mod 25", crt_lifts([25], 3000, 200))
survey("CRT-lift ==AP mod lcm(13,25)", crt_lifts([13,25], 3000, 50))
survey("CRT-lift ==AP mod lcm(2..12,25)", crt_lifts(list(range(2,13))+[25], 3000, 20))

# 3. mod-25-saturated blockers (case 2 -- the hard residual): families with M_25 < 2/25, non-AP
def sat_blockers(ntry, hmax):
    fams = []
    for _ in range(ntry):
        W = sorted(sample(range(1, hmax), 12))
        if is_prim(W) and W != AP and Mq(W, 25) < TWO25:
            fams.append(W)
    return fams
survey("mod-25-saturated blockers (h<=300)", sat_blockers(60000, 300))
survey("mod-25-saturated blockers (h<=3000)", sat_blockers(60000, 3000))

print("\nREADING: if EVERY non-AP family (incl. high-height + CRT-lifts + mod-25 blockers) clears at")
print("some bounded q (<=39, ideally <=14), the (C) covering is height-uniform -- the open node holds.")
print("A single covering GAP would refute kps S44's finite covering. The max clearing-q is the mechanism.")
print("DONE")
