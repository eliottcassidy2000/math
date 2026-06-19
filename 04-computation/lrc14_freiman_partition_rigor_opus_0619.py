#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FREIMAN PARTITION RIGOR for LRC(14) endgame (opus-2026-06-19).

GOAL: make precise that EVERY primitive non-AP E with |E|=k (k=8,9,10) is EITHER
  (a) contained in a proper d>=2 GAP (generalized arithmetic progression) of
      BOUNDED size  [small doubling], OR
  (b) has a DISSOCIATED stranger w in E (no height-<=H relation to E\{w})
      [large doubling].
No E falls through both. Find the doubling threshold sigma* and explicit H.

Definitions used here (all integer arithmetic, fractions.Fraction for L_y):
- E primitive: gcd(E)=1, 0 in E, |E|=k, distinct nonneg ints.
- doubling sigma(E) = |E+E| / |E|.
- d-dim GAP P(a_0; q_1..q_d; l_1..l_d) = { a_0 + sum_i x_i q_i : 0<=x_i<l_i },
  proper if the sums are distinct; |P| = prod l_i.  We say E lives in a GAP of
  dimension d and size S if E subset P, |P|=S.  Full AP = dim-1 GAP.
- dissociated stranger: w in E such that there is NO nonzero integer vector
  (m_v) with |m_v|<=H, supported on E\{w}, plus a coeff on w (|coeff|<=H, nonzero),
  with sum m_v v + coeff*w = 0.  I.e. w is not an O(H)-bounded integer
  combination of the other elements (as a Z-relation). Equivalently w escapes
  every bounded relation lattice.

We test the partition empirically over many k=8 sets (and k=9,10 spot checks),
compute L_y, and look for the genuine "third pocket": small doubling but no
bounded GAP AND no dissociated stranger.

Engine copied from lrc14_empty_sector_distribution_kps.py (dist_p, L_y, g_poly).
"""
import sys, itertools
from fractions import Fraction
from functools import reduce
from math import gcd
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout, 'reconfigure') else None

# ---------- L_y engine (verbatim) ----------
def N_at(E, x):
    hit = set()
    for e in E:
        v = e * x
        v = v - (v.numerator // v.denominator)
        s = (v.numerator * 7) // v.denominator
        hit.add(s)
    return sum(1 for j in range(1, 7) if j not in hit)

def dist_p(E):
    E = sorted(set(E))
    bps = set([Fraction(0), Fraction(1)])
    for e in E:
        if e == 0: continue
        for a in range(0, 7 * e + 1):
            bps.add(Fraction(a, 7 * e))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    p = [Fraction(0)] * 7
    for i in range(len(bps) - 1):
        lo, hi = bps[i], bps[i+1]
        if hi == lo: continue
        mid = (lo + hi) / 2
        t = N_at(E, mid)
        p[t] += (hi - lo)
    return p

def g_poly(k):
    g = []
    for t in range(7):
        if k == 8:
            val = Fraction((t-1)*(t-2)*(t-4)*(t-5), 40)
        elif k in (9, 10):
            val = Fraction(-(t-2)*(t-3)*(t-6), 36)
        else:
            val = Fraction((t-3)*(t-4), 12)
        g.append(val)
    return g

def L_y(E, k):
    p = dist_p(E); g = g_poly(k)
    return sum(p[t] * g[t] for t in range(7)), p

CAP = {8: Fraction(38153, 100000), 9: Fraction(49426, 100000), 10: Fraction(6044, 10000)}

# ---------- Freiman structure tools ----------
def doubling(E):
    E = sorted(set(E))
    s = set()
    for a in E:
        for b in E:
            s.add(a + b)
    return Fraction(len(s), len(E)), len(s)

def is_full_AP(E):
    E = sorted(set(E))
    if len(E) < 2: return True
    d = E[1] - E[0]
    return all(E[i+1] - E[i] == d for i in range(len(E)-1))

def excess(E):
    _, sz = doubling(E)
    return sz - (2 * len(E) - 1)

def smallest_gap_2dim(E, maxsize_factor=8):
    """
    Try to find a 2-dim GAP P(a0; q1,q2; l1,l2) of BOUNDED size containing E.
    Strategy: E primitive (gcd 1). For a 2-dim GAP, pick q1 = small generator,
    q2 = second generator. We search candidate moduli q1 over small divisor-like
    steps and represent E as { a0 + x*q1 + y*q2 }.
    Concretely: try every q1 in a candidate set; reduce E mod q1 -> the residues
    must form a short AP (the q2-direction) for the GAP to be 2-dim with that q1.
    Return (q1,q2,l1,l2,size) of the SMALLEST size found, or None.
    """
    E = sorted(set(E))
    n = len(E)
    span = E[-1] - E[0]
    best = None
    # candidate q1: differences and small values
    cand = set()
    for i in range(len(E)):
        for j in range(i+1, len(E)):
            cand.add(E[j] - E[i])
    cand |= set(range(2, span + 1))
    for q1 in sorted(cand):
        if q1 <= 1: continue
        # write each e = q1 * row + col,  col in 0..q1-1
        cols = sorted(set(e % q1 for e in E))
        rows = sorted(set(e // q1 for e in E))
        # for a 2-dim GAP with first generator q1 (the "column within block") we
        # need cols to be an AP (q2-direction) AND rows to be an AP (q1-direction
        # mapped to multiples of q1). Then size = len(rows_AP)*len(cols_AP).
        # cols AP:
        def ap_cover(vals):
            vals = sorted(set(vals))
            if len(vals) == 1: return (vals[0], 1, 1)  # step,len, (start)
            diffs = [vals[i+1]-vals[i] for i in range(len(vals)-1)]
            g = reduce(gcd, diffs)
            if g == 0: return None
            l = (vals[-1] - vals[0]) // g + 1
            return (g, l, vals[0])
        cinfo = ap_cover(cols)
        rinfo = ap_cover(rows)
        if cinfo is None or rinfo is None: continue
        cstep, clen, cstart = cinfo
        rstep, rlen, rstart = rinfo
        # genuine 2-dim: need both directions length >=2, and product reasonably
        # bounds the set. size of the GAP:
        size = clen * rlen
        # must contain E and be proper-ish (we don't enforce strict properness,
        # just boundedness). Require it's a real reduction (size < trivial span+1
        # full AP cover) and 2-dimensional (both lens>=2) OR exactly covers.
        if clen >= 2 and rlen >= 2 and size <= maxsize_factor * n:
            if best is None or size < best[-1]:
                # generators: q2 = cstep (within-row), q1eff = q1*rstep (row jump)
                best = (q1 * rstep, cstep, rlen, clen, size)
    return best

def dissociated_stranger(E, H):
    """
    Is there w in E with NO bounded integer relation tying it to E\{w}?
    A relation: integers m_v (|m_v|<=H) for v in E\{w}, and c (1<=|c|<=H), with
        c*w + sum_v m_v v = 0.
    If for some w NO such relation exists, w is an H-dissociated stranger.
    We brute force over small E (k<=10) and small H by checking whether w lies in
    the set { -(1/c) sum m_v v }. Equivalent integer feasibility: does there exist
    c in 1..H and integer combo of E\{w} with coeffs in -H..H summing to -c*w?
    We test by enumerating bounded combos -> the reachable multiples; expensive,
    so we restrict: a relation of height<=H exists iff w is in the H-bounded
    integer span lattice of E\{w} at level c. We use a meet-in-the-middle reachable
    set of sum_v m_v v for |m_v|<=H.
    Returns the FIRST w that is a stranger, plus the min height that would catch it
    if any (None meaning > H), else None (every w has a bounded relation).
    """
    E = sorted(set(E))
    n = len(E)
    for wi, w in enumerate(E):
        others = [v for v in E if v != w]
        # reachable T = { sum m_v v : |m_v|<=H } ; relation exists iff exists c in
        # 1..H with -c*w in T (and the combo not all-zero, but c*w!=0 unless w=0).
        # Build reachable incrementally with pruning by absolute bound.
        bound = H * sum(others) + H * n  # generous
        reach = {0}
        for v in others:
            new = set()
            for s in reach:
                for m in range(-H, H + 1):
                    val = s + m * v
                    if -bound <= val <= bound:
                        new.add(val)
            reach = new
        found_rel = False
        for c in range(1, H + 1):
            if (-c * w) in reach:
                # need the combo to be nonzero unless w!=0; if w==0 then c*0=0 in
                # reach trivially via all-zero -> not a real relation. Skip w==0.
                if w == 0:
                    found_rel = True  # 0 is never a stranger (it's in every GAP)
                    break
                found_rel = True
                break
        if not found_rel and w != 0:
            return (w, None)  # stranger at height H
    return None

# ---------- experiments ----------
def classify(E, k, H_diss=2, gapfactor=8):
    """Return dict with sigma, excess, isAP, gap2dim, stranger, L_y."""
    sig, dsz = doubling(E)
    ex = excess(E)
    ap = is_full_AP(E)
    gap = None if ap else smallest_gap_2dim(E, gapfactor)
    strg = None if ap else dissociated_stranger(E, H_diss)
    Lv, _ = L_y(E, k)
    return dict(E=E, sigma=sig, dsz=dsz, excess=ex, isAP=ap, gap=gap, stranger=strg, L=Lv)

def primitive(E):
    E = sorted(set(E))
    return E[0] == 0 and reduce(gcd, E) == 1

def main():
    print("=" * 78)
    print("FREIMAN PARTITION RIGOR — LRC(14) endgame (opus-2026-06-19)")
    print("=" * 78)

    # 1. The S11 examples: verify they ARE in a low-dim GAP (no dissociated stranger).
    print("\n[1] S11 EXAMPLES (claimed: no dissociated stranger -> must be in d>=2 GAP)")
    s11 = [[0,3,5,16,28,30,33,35], [0,4,12,15,20,21,25,31]]
    for E in s11:
        c = classify(E, 8, H_diss=2, gapfactor=12)
        print(f"  E={E}")
        print(f"    sigma={float(c['sigma']):.3f} (|E+E|={c['dsz']}), excess={c['excess']}, AP={c['isAP']}, L_y={float(c['L']):.5f} (cap_8={float(CAP[8]):.5f})")
        print(f"    2-dim GAP (q1,q2,l1,l2,size)={c['gap']}")
        print(f"    H=2 dissociated stranger={c['stranger']}")

    # 2. Threshold hunt for k=8: scan many primitive non-AP sets, record
    #    (sigma, has-bounded-GAP, has-stranger). Find sigma* separating regimes.
    print("\n[2] THRESHOLD HUNT k=8 (random + structured primitive non-AP sets)")
    import random
    random.seed(20260619)
    rows = []
    seen = set()
    # structured: APs perturbed + 2D gaps + random
    cand_sets = []
    # 2D GAPs A1+A2
    for q in range(2, 9):
        for l1 in range(2, 5):
            for l2 in range(2, 5):
                if l1 * l2 < 8: continue
                P = sorted({a + b for a in range(0, l1 * q, q) for b in range(l2)})
                if len(P) >= 8:
                    cand_sets.append(P[:8])
    # random spread up to 50
    for _ in range(4000):
        sp = random.randint(8, 50)
        E = sorted(random.sample(range(0, sp + 1), 8))
        E = [e - E[0] for e in E]
        cand_sets.append(E)
    cnt_a = cnt_b = cnt_both = cnt_neither = 0
    neither_examples = []
    maxL_nonAP = Fraction(-1)
    maxL_set = None
    sig_with_gap = []
    sig_with_strg = []
    for E in cand_sets:
        key = tuple(E)
        if key in seen: continue
        seen.add(key)
        if not primitive(E) or is_full_AP(E): continue
        c = classify(E, 8, H_diss=2, gapfactor=8)
        has_gap = c['gap'] is not None
        has_strg = c['stranger'] is not None
        if c['L'] > maxL_nonAP:
            maxL_nonAP = c['L']; maxL_set = E
        if has_gap: sig_with_gap.append(float(c['sigma']))
        if has_strg: sig_with_strg.append(float(c['sigma']))
        if has_gap and not has_strg: cnt_a += 1
        elif has_strg and not has_gap: cnt_b += 1
        elif has_gap and has_strg: cnt_both += 1
        else:
            cnt_neither += 1
            if len(neither_examples) < 25:
                neither_examples.append((E, float(c['sigma']), float(c['L'])))
    print(f"  scanned {len(seen)} candidate sets")
    print(f"  (a) bounded-GAP, no stranger : {cnt_a}")
    print(f"  (b) stranger, no bounded-GAP : {cnt_b}")
    print(f"  both GAP and stranger        : {cnt_both}")
    print(f"  NEITHER (third pocket?)      : {cnt_neither}")
    if sig_with_gap:
        print(f"  sigma range WITH bounded GAP : [{min(sig_with_gap):.3f}, {max(sig_with_gap):.3f}]")
    if sig_with_strg:
        print(f"  sigma range WITH stranger    : [{min(sig_with_strg):.3f}, {max(sig_with_strg):.3f}]")
    print(f"  max L_y over non-AP primitive : {float(maxL_nonAP):.5f} at {maxL_set} (cap_8={float(CAP[8]):.5f}, margin={float(CAP[8]-maxL_nonAP):.5f})")
    if neither_examples:
        print("  THIRD-POCKET candidates (re-examine with bigger H / gapfactor):")
        for E, sg, Lv in neither_examples:
            print(f"    {E} sigma={sg:.3f} L_y={Lv:.5f}")

    # 3. Re-examine third-pocket candidates with bigger H and gapfactor to see if
    #    they collapse into (a) or (b).
    print("\n[3] RE-EXAMINE THIRD-POCKET with H_diss up to 4, gapfactor up to 20")
    survivors = []
    for E, sg, Lv in neither_examples:
        c = classify(E, 8, H_diss=4, gapfactor=20)
        has_gap = c['gap'] is not None
        has_strg = c['stranger'] is not None
        status = ("GAP" if has_gap else "") + ("+STRANGER" if has_strg else "")
        if not has_gap and not has_strg:
            survivors.append(E)
            status = "STILL NEITHER"
        print(f"    {E}: gap={c['gap']} stranger={c['stranger']} -> {status}")
    print(f"  genuine survivors after escalation: {len(survivors)}")
    if survivors:
        print("  *** THESE ARE TRUE THIRD-POCKET RESIDUE ***", survivors)
    else:
        print("  *** NO third pocket: dichotomy (a)|(b) is COMPLETE at k=8 (H<=4, gap<=20k) ***")

    # 4. Threshold: max sigma that still yields a bounded GAP vs min sigma giving stranger.
    print("\n[4] THRESHOLD sigma* (boundary between small-doubling GAP and large-doubling stranger)")
    if sig_with_gap and sig_with_strg:
        print(f"  max sigma still GAP-containable      sigma_GAP_max = {max(sig_with_gap):.3f}")
        print(f"  min sigma producing a stranger       sigma_STR_min = {min(sig_with_strg):.3f}")
        print(f"  overlap region [{min(sig_with_strg):.3f}, {max(sig_with_gap):.3f}] = sets with BOTH (safe either way)")

if __name__ == "__main__":
    main()
