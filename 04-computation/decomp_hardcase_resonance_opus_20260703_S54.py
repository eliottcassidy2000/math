#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
IS THE HARD DECOMPOSITION CASE FORCED TO BE RESONANT?  (structural decomposition of hlarge, OPEN-Q toward finishing)

opus-2026-07-03-S54. The renormalization architecture (HYP-4041): (1) ratio<=13 => spread13 [DONE];
(2) else near-equal cluster => THM-608/scale_separation [DONE]; (3) recurse; (4) bounded => census.
The OPEN residual (mac-mini): the "slow-runner (e.g. 1) + WIDE-far spread" that violates THM-608's
near-equal condition (ii). QUESTION: is that wide-far part FORCED by covering to be a 13-spaced (resonant,
Eisenstein Phi_6) comb -- so my scale_separation_phase (S53) closes it? If the hard covering families are
resonant, the decomposition closes: ratio<=13 -> spread13; near-equal cluster -> THM-608; resonant comb ->
scale_separation_phase; bounded -> census.

METHOD: search covering 13-families with a SLOW runner (contains 1) + LARGE far runners (magnitude > M),
that are NOT spread13 (ratio > 13). For each, find its lonely witness t = a/q (max_t min_v), record the
denominator q and whether q is Eisenstein (q = Phi_6(n) = n^2-n+1 for some n, or a divisor/multiple), and
whether the far runners are 13-spaced / share a residue structure (the resonance signature 14*13=-1 mod 183
generalizes to: far runners w_i with w_i * t = tight AP).
"""
import sys, itertools
from fractions import Fraction as Fr
from math import gcd
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
def norm(x):
    f = x - (x.numerator // x.denominator); return min(f, 1 - f)
def covering(S):
    return all(any(v % q == 0 for v in S) for q in range(2, 15))
def M_and_witness(S, Dmax):
    dens = set()
    for v in S: dens.add(2 * v)
    for v in S:
        for w in S:
            dens.add(v + w)
            if v != w: dens.add(abs(v - w))
    best = Fr(0); bt = None
    for b in sorted(d for d in dens if 0 < d <= Dmax):
        for a in range(1, b):
            if gcd(a, b) != 1: continue
            md = min(norm(Fr(a * v, b)) for v in S)
            if md > best: best = md; bt = Fr(a, b)
    return best, bt
def is_eisenstein(q):
    # q = Phi_6(n) = n^2-n+1 for some n, or q | such (the covering-min denominators are Phi_6-flavored)
    hits = []
    for n in range(2, 60):
        p6 = n * n - n + 1
        if p6 == q: hits.append(("=Phi6(%d)" % n, p6))
        elif p6 % q == 0 and p6 != q: hits.append(("|Phi6(%d)" % n, p6))
        elif q % p6 == 0 and q != p6: hits.append(("Phi6(%d)|" % n, p6))
    return hits

print("=" * 98)
print(" HARD DECOMPOSITION CASE: covering 13-families with a slow runner {1} + wide-far, NOT spread13")
print("=" * 98)
# Construct: base {1..k} (contains the slow runner 1) + far runners that make it covering with a big gap.
# Far runners chosen as multiples covering the high moduli, with magnitude >> base (ratio > 13).
# Use the deep-well-family generator: {1..12} + far runner = 13*m (covers 13) tuned to also hit 14.
results = []
# family A: the deep well {1..12, 182} (r=1 far) -- the canonical resonant case
famsA = [ tuple(sorted(set(range(1,13)) | {182})) ]
# family B: slow base {1..6} + far block covering 7..14 (>=7 far, the hard >=7-far case)
#   need multiples of 7,8,9,10,11,12,13,14 among 7 far runners > 100. Use q*round(K/q) alignment (mac-mini band-blocker).
K = 2*3*5*7*11*13  # 30030
for scale in [1]:
    far = []
    for q in [7,8,9,10,11,12,13,14]:
        far.append(q * round(K*scale / q))   # aligned far runner divisible by q, near K*scale
    far = sorted(set(far))
    base = [1,2,3,4,5,6]
    S = tuple(sorted(set(base) | set(far)))
    if len(S) == 13 and covering(S):
        famsB = [S]
    else:
        # trim/pad to 13 covering
        famsB = []
        for drop in itertools.combinations(base, max(0,len(S)-13)):
            S2 = tuple(x for x in S if x not in drop)
            if len(S2)==13 and covering(S2): famsB.append(S2); break
        if not famsB and len(S) < 13:
            famsB = [tuple(sorted(set(S) | set(range(1, 1+ (13-len(S))))))] if covering(tuple(sorted(set(S)|set(range(1,1+(13-len(S))))))) else []

for tag, fams in [("A deep-well r=1", famsA), ("B >=7-far aligned", famsB)]:
    for S in fams:
        if len(S) != 13 or not covering(S):
            print(f"  [{tag}] skip (|S|={len(S)}, covering={covering(S) if len(S)<=13 else '?'})"); continue
        mx, mn = max(S), min(S)
        ratio = mx / mn
        M, wt = M_and_witness(S, 2*mx + 2)
        q = wt.denominator if wt else None
        eis = is_eisenstein(q) if q else []
        far = [v for v in S if v > 13]
        # far residues mod 13 and the pairwise differences (13-spacing signature)
        far_mod13 = sorted(set(v % 13 for v in far))
        far_diffs = sorted(set(far[j]-far[i] for i in range(len(far)) for j in range(i+1,len(far))))
        spaced13 = all(d % 13 == 0 for d in far_diffs) if len(far) >= 2 else (len(far)==1)
        print(f"\n  [{tag}] S = {S}")
        print(f"     ratio {mx}/{mn} = {ratio:.1f} (spread13? {ratio<=13}); #far(>13) = {len(far)}: {far}")
        print(f"     LONELY: M = {M} = {float(M):.5f} at t = {wt} (q = {q}); >=1/14? {M>=Fr(1,14)}")
        print(f"     q Eisenstein-Phi_6? {eis if eis else 'no'}")
        print(f"     far mod 13 = {far_mod13} (single residue = 13-spaced-collapse? {len(far_mod13)==1}); 13-spaced diffs? {spaced13}")
        results.append((tag, S, ratio, len(far), q, bool(eis), spaced13))

print("\n" + "=" * 98)
print(" READING")
print("=" * 98)
print("  If the hard (>=7-far, ratio>13) covering families are LONELY at an EISENSTEIN Phi_6 denominator with")
print("  a 13-spaced (single-residue-mod-13) far part, then the decomposition CLOSES: every hard case is")
print("  resonant => scale_separation_phase (S53). If instead they are lonely at generic q, the far part is")
print("  handled by the near-equal renormalization (THM-608) or splits further. Either way each far cluster")
print("  is one of the two peelable types (near-equal / resonant) -- the structural-decomposition claim.")
for r in results:
    print(f"   {r[0]}: ratio {r[2]:.1f}, #far {r[3]}, q={r[4]}, Eisenstein={r[5]}, 13-spaced-far={r[6]}")
print("DONE.")
