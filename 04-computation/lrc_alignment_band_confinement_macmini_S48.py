#!/usr/bin/env python3
"""
mac-mini-2026-07-04-S48 -- HYP-4094: the ALIGNMENT-BAND GEOMETRY of the
loose-base confinement (hcomp's sole remaining case).

Setup: primitive compressed covering 13-set V = B ∪ {w}, base B loose
(M(B) >= 2/25 by the S47 spectrum), killer w <= 13 * max(B) (compressed).
M(V) < 1/13 requires the killer's radius-1/13 comb to FULLY COVER the base's
lonely set L_B(1/13).

GEOMETRY: the killer's teeth are DISJOINT (width 2/(13w), gap 11/(13w));
an interval is fully covered iff it fits inside ONE tooth:
  component J = [a,b] covered  <=>  exists k: w*a >= k - 1/13 and w*b <= k + 1/13
  <=>  ||w * c_J||_dist <= 1/13 - w*|J|/2    (c_J = midpoint)
=> for EACH component: a residue BAND on w mod den(c_J); ALL components
simultaneously: a CRT band system; plus covering pins (q | w for q missed by
B); inside the compressed window. Confinement = the intersection is empty or
a small finite candidate list, each checked exactly.

This script runs the whole pipeline EXACTLY on live loose bases.
"""
from fractions import Fraction as F
from math import gcd
import sys
sys.path.insert(0, '04-computation')
from lonely_profile import profile

R13 = F(1, 13)

def dist(x):
    f = x - int(x)
    if f < 0: f += 1
    return min(f, 1 - f)

def components_at(S, r):
    """Exact components of L_S(r) as [a,b] closed intervals (via danger union)."""
    iv = []
    for v in S:
        half = F(r) / v
        for k in range(v):
            c = F(k, v); a, b = c - half, c + half
            if a < 0: iv.append((a + 1, F(1))); iv.append((F(0), b))
            elif b > 1: iv.append((a, F(1))); iv.append((F(0), b - 1))
            else: iv.append((a, b))
    iv.sort()
    merged = []
    for a, b in iv:
        if merged and a <= merged[-1][1]:
            if b > merged[-1][1]: merged[-1] = (merged[-1][0], b)
        else: merged.append((a, b))
    comps = []
    for i in range(len(merged)):
        b_cur = merged[i][1]
        a_next = merged[(i + 1) % len(merged)][0] + (1 if i + 1 == len(merged) else 0)
        if a_next > b_cur:
            comps.append((b_cur, a_next))
    return comps

def band_residues(cJ, halfwidth_fn, wmax):
    """Admissible w in [1, wmax] with dist(w*cJ) <= 1/13 - w*|J|/2 (w-dependent)."""
    ok = []
    for w in range(1, wmax + 1):
        beta = R13 - halfwidth_fn * w
        if beta < 0: continue
        if dist(w * cJ) <= beta:
            ok.append(w)
    return ok

def confinement_pipeline(B, label):
    B = sorted(B)
    pB = profile(B, F(1, 3))
    MB = pB.M()
    comps = components_at(B, R13)
    b = max(B)
    wmax = 13 * b   # compressed window (v_max <= 13 * v_second; v_second <= max(B))
    print(f"--- base {label}: {B}")
    print(f"    M(B) = {MB} = {float(MB):.6f}; components of L_B(1/13): {len(comps)}")
    # covering pins: q in 2..14 with no multiple in B must divide w
    pins = [q for q in range(2, 15) if not any(v % q == 0 for v in B)]
    print(f"    covering pins on w: multiples of {pins}")
    # per-component admissible sets
    sets = []
    for (a, bb) in comps:
        cJ = (a + bb) / 2
        halfw = (bb - a) / 2
        adm = set(band_residues(cJ, halfw, wmax))
        sets.append(adm)
    inter = set(range(1, wmax + 1))
    for s in sets: inter &= s
    for q in pins:
        inter = {w for w in inter if w % q == 0}
    # compressed + primitive + w > max(B) (killer is the max) + V covering
    cands = sorted(w for w in inter if w > b and gcd(w, __import__('math').gcd(*B) if len(B)>1 else B[0]) >= 1)
    cands = [w for w in cands if all(any(v % q == 0 for v in B + [w]) for q in range(2, 15))]
    from functools import reduce
    cands = [w for w in cands if reduce(gcd, B + [w]) == 1]
    print(f"    band-intersection candidates in window (w <= {wmax}): {len(cands)}"
          f"{' -- EMPTY: loose case CLOSED for this base' if not cands else ''}")
    bad = []
    for w in cands:
        V = B + [w]
        pV = profile(V, F(1, 3))
        MV = pV.M()
        if MV is not None and MV < R13:
            bad.append((w, MV))
    if cands:
        print(f"    exact check of candidates: {len(bad)} with M(V) < 1/13 "
              f"{[(w, str(m)) for w, m in bad] if bad else '-- ALL SAFE: confinement holds'}")
    return len(cands), bad

print("=" * 78)
print("ALIGNMENT-BAND CONFINEMENT: loose bases (M(B) >= 2/25), compressed window")
print("=" * 78)
total_bad = 0
for B, lab in [
    ([1,2,3,4,5,6,7,8,9,10,11,24], "the 2/25 base {1..11,24}"),
    ([1,2,3,4,5,6,7,8,9,10,11,13], "{1..11,13}"),
    ([1,2,3,4,5,6,7,8,9,10,12,13], "{1..10,12,13}"),
    ([2,3,4,5,6,7,8,9,10,11,12,13], "{2..13}"),
]:
    n, bad = confinement_pipeline(B, lab)
    total_bad += len(bad)
print(f"\nTOTAL violations across tested loose bases: {total_bad}")
print("(criterion: M(V) < 1/13 forces w into the CRT band intersection; the")
print(" pipeline enumerates it exactly and checks each candidate.)")
