#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_q108_diverse_angles_kpswf3.py   (kind-pasteur 2026-06-21, THREAD D)

DIVERSE new-angle generator for the LRC(14) cover bound.  Six genuinely distinct
attack angles, each with a 5-minute exact-arithmetic probe and a promise rating.

Core object (canon):  for a finite set E of nonneg integers,
   measS7(E) = Lebesgue measure of { x in [0,1) : {floor(7*frac(e*x)) : e in E} == Z/7 }
            = p0(E)  = the "cover probability" the LRC(14) cover bound is about.
We need:  for "wide" E (in particular |E|=k with a balanced/multi-cluster shape),
   p0(E) <= cap_k    (cap_8=2243/5880, cap_9=1979/4004, cap_10=55/91).

This file probes 6 NEW lenses (NOT Delsarte-LP / torus-discrepancy / death-chain /
Freiman, which are the existing main threads):

 (1) MATROID / TUTTE of the relation lattice Lambda(E)   [cover encoded as a matroid]
 (2) TOURNAMENT-SIDE bound: Lambda(E) is a cycle space (HYP-2719) -> do H/OCF/chromatic
     invariants bound corr?
 (3) TRANSFER-OPERATOR / death-chain spectral gap 1/7 -> exponential mixing bound
 (4) KOKSMA-HLAWKA with an explicit bounded-variation test function for the cover event
 (5) LP / Delsarte directly on the CORRELATED cover (not the decorrelated plateau)
 (6) FOURIER / GAUSS-SUM bound exploiting QR(7)

EXACT rational arithmetic where measures are involved.
"""
import itertools
from fractions import Fraction as Fr
from math import gcd, cos, pi
import cmath

P = 7

# ============================================================================
# CORE: exact measS7(E) via breakpoint cells in x.
# ============================================================================
def sector(yf):  # yf a Fraction in [0,1)
    return int(P * yf)

def breakpoints(E):
    """All x in [0,1) where some floor(7 e x) jumps: x = t/(7 e)."""
    bp = {Fr(0), Fr(1)}
    for e in E:
        if e == 0:
            continue
        for t in range(0, P * e):
            bp.add(Fr(t, P * e))
    return sorted(bp)

def cell_sectorset(E, mid):
    return frozenset(sector((e * mid) % 1) for e in E)

def measS7(E):
    """Exact rational measure of x with full 7-sector cover by E."""
    E = [int(e) for e in E if int(e) != 0]
    xs = breakpoints(E)
    total = Fr(0)
    for a, b in zip(xs, xs[1:]):
        mid = (a + b) / 2
        if len(cell_sectorset(E, mid)) == P:
            total += (b - a)
    return total

CAP = {8: Fr(2243, 5880), 9: Fr(1979, 4004), 10: Fr(55, 91)}

# ============================================================================
# Helper: the per-x "sector multiset profile" -- which residues each e lands in.
# ============================================================================
def profile_cells(E):
    """Return list of (length, tuple-of-sectors) over x-cells. tuple aligned to E order."""
    E = [int(e) for e in E if int(e) != 0]
    xs = breakpoints(E)
    out = []
    for a, b in zip(xs, xs[1:]):
        mid = (a + b) / 2
        secs = tuple(sector((e * mid) % 1) for e in E)
        out.append((b - a, secs))
    return out


# ============================================================================
# ANGLE 1 -- MATROID / TUTTE of the cover incidence.
# ----------------------------------------------------------------------------
# Idea: For a fixed x, e "covers sector s" iff floor(7 e x)=s.  The cover event is
# "the 7 sectors are all hit".  Think of the bipartite incidence  E x Z/7  with an
# edge (e,s) present on the x-cell where e lands in s.  Covering all 7 sectors is a
# TRANSVERSAL / set-cover.  The relevant matroid is the TRANSVERSAL MATROID of the
# bipartite "lands-in" relation; the cover event = the ground set Z/7 is spanned.
# Probe: compute, per x-cell, the bipartite incidence and ask whether the cover
# event is governed by a matroid rank (rank of the transversal matroid on Z/7).
# If cover <=> transversal-matroid-rank = 7, then p0 = measure{rank=7}, and the
# Tutte/rank generating function might give a clean cap.  We test: is "rank=7"
# equivalent to "cover", and does the rank distribution have low-order structure?
# ============================================================================
def angle1_matroid(E):
    cells = profile_cells(E)
    # For each cell: the set of sectors hit = image of E. cover <=> |image|=7.
    # transversal matroid on Z/7 with sets {sector(e)} (singletons) has rank = |image|.
    # so rank=7 <=> cover  TRIVIALLY (singletons).  The interesting matroid is on E:
    # the "partition matroid" where e's are grouped by which sector they hit; a
    # transversal picking one e per sector.  cover <=> a system of distinct reps exists
    # for the 7 sectors among the e->sector assignment <=> all 7 sectors are images.
    # => still trivial. Tutte of a partition matroid factorizes: T = prod_s (stuff).
    # Test the claimed factorization of p0 over independent "sector-blocks":
    rankdist = {}
    for length, secs in cells:
        r = len(set(secs))
        rankdist[r] = rankdist.get(r, Fr(0)) + length
    p0 = rankdist.get(7, Fr(0))
    return p0, rankdist


# ============================================================================
# ANGLE 2 -- TOURNAMENT-SIDE bound.
# ----------------------------------------------------------------------------
# HYP-2719: Lambda(E) IS a cycle space.  A natural tournament-side handle: build,
# for each x-cell, the "collision graph" G_x on E where e~e' iff sector(e)=sector(e')
# (they land in the same residue).  cover (all 7 hit) <=> the collision graph has
# exactly 7 color classes used out of |E| vertices, i.e. the PARTITION of E into
# sector-classes uses all 7 colors.  The number of sector-classes = chromatic-like
# invariant.  Probe: is p0 controlled by E[ #distinct-sector-classes ]?  Compute the
# distribution of (#sector classes) over x, and the "collision count" alpha (pairs in
# same sector) -- the tournament-side OCF analog.  If a Turan-type bound
#   #classes <= 7  with equality measure = p0  AND  p0 <= f(alpha) for an explicit f,
# we have a tournament-flavored bound.
# ============================================================================
def angle2_tournament(E):
    cells = profile_cells(E)
    k = len(E)
    # E[#classes], E[#collision-pairs], and p0
    Eclasses = Fr(0)
    Ecoll = Fr(0)
    p0 = Fr(0)
    classdist = {}
    for length, secs in cells:
        c = len(set(secs))
        classdist[c] = classdist.get(c, Fr(0)) + length
        Eclasses += length * c
        # collisions = sum_s C(count_s,2)
        cnt = {}
        for s in secs:
            cnt[s] = cnt.get(s, 0) + 1
        coll = sum(v * (v - 1) // 2 for v in cnt.values())
        Ecoll += length * coll
        if c == P:
            p0 += length
    return p0, Eclasses, Ecoll, classdist


# ============================================================================
# ANGLE 3 -- TRANSFER-OPERATOR / DEATH-CHAIN spectral gap.
# ----------------------------------------------------------------------------
# Process the elements e_1<...<e_k in increasing order; state = set of sectors hit so
# far (subset of Z/7).  Adding e_i moves the "hit-set" by adding sector(e_i x).  Over
# random x the increments are NOT independent, but in the DECORRELATED (well-separated)
# regime each new far element adds a ~uniform sector independent of the past.  Then the
# hit-set is a coupon-collector / "death chain" on subsets of Z/7 with absorbing state
# "all 7".  The complement count (#missing) is a pure-death chain m -> m with kill
# prob m/7 each step.  cover-by-k = P(absorbed by step k).  This gives an EXPLICIT
# upper bound p0 <= P(coupon collector collects 7 in <=k draws) when increments are
# (sub)uniform.  Probe: compute the iid coupon bound C(k) and compare to the true
# decorrelated plateau P2-like and to cap_k.  The "spectral gap 1/7" = smallest death
# rate = the slowest mode.
# ============================================================================
def coupon_cover_prob(k, P=7):
    """Exact P(7 uniform-coupon types all seen in k iid draws) via inclusion-exclusion."""
    # = sum_{j=0}^{7} (-1)^j C(7,j) ((7-j)/7)^k
    from math import comb
    tot = Fr(0)
    for j in range(P + 1):
        tot += (-1) ** j * comb(P, j) * Fr((P - j) ** k, P ** k)
    return tot

def angle3_deathchain(kmax=12):
    rows = []
    for k in range(7, kmax + 1):
        c = coupon_cover_prob(k)
        rows.append((k, c))
    return rows


# ============================================================================
# ANGLE 4 -- KOKSMA-HLAWKA with explicit bounded-variation test function.
# ----------------------------------------------------------------------------
# The cover indicator 1[cover] is the integrand.  KH says |integral - average over a
# point set| <= V(f) * D*(points).  But the cleaner KH move: write 1[all 7 sectors hit]
# as a product/inclusion-exclusion of 1D indicators and bound EACH missing-sector
# probability.  p0 = P(no sector missed) = 1 - P(union of "sector s missed").
# By Bonferroni (truncated inclusion-exclusion), an UPPER bound on p0 uses the EVEN
# partial sums:  p0 <= 1 - S1 + S2  where S1 = sum_s P(s missed), S2 = sum_{s<s'}
# P(s,s' both missed).  Each term is an exact measure (a union of arcs).  This is a
# rigorous, cheap UPPER bound that needs only pairwise miss-probabilities.  Probe:
# compute the Bonferroni-2 upper bound and compare to cap_k; if it already beats cap,
# we have an elementary proof route via 2nd-order Bonferroni.
# ============================================================================
def miss_prob(E, missing):
    """exact measure of x such that NONE of E lands in any sector in `missing` set."""
    E = [int(e) for e in E if int(e) != 0]
    xs = breakpoints(E)
    miss = set(missing)
    total = Fr(0)
    for a, b in zip(xs, xs[1:]):
        mid = (a + b) / 2
        secs = set(sector((e * mid) % 1) for e in E)
        if not (secs & miss):   # none of E hit a missing sector => those sectors are empty
            total += (b - a)
    return total

def angle4_bonferroni(E):
    # S1 = sum_s P(sector s empty); S2 = sum_{s<s'} P(s,s' both empty)
    S1 = Fr(0)
    for s in range(P):
        S1 += miss_prob(E, {s})
    S2 = Fr(0)
    for s in range(P):
        for t in range(s + 1, P):
            S2 += miss_prob(E, {s, t})
    S3 = Fr(0)
    for s in range(P):
        for t in range(s + 1, P):
            for u in range(t + 1, P):
                S3 += miss_prob(E, {s, t, u})
    p0_true = measS7(E)
    lower_b2 = 1 - S1 + S2 - S3          # Bonferroni order-3 (lower bound on p0)
    upper_b2 = 1 - S1 + S2               # Bonferroni order-2 (upper bound on p0)
    upper_b0 = 1 - S1                     # order-1 (upper bound, weakest)
    return p0_true, upper_b0, upper_b2, lower_b2, S1, S2, S3


# ============================================================================
# ANGLE 5 -- LP / DELSARTE on the CORRELATED cover (sketch / diagnostic).
# ----------------------------------------------------------------------------
# The decorrelated Delsarte LP (mac-mini) bounds the plateau P2.  Here we test whether
# a *correlated* refinement is even needed: compute, for the worst balanced bases, the
# gap between true p0 and the decorrelated plateau, i.e. how much correlation actually
# contributes.  If correlation only ever DECREASES p0 (corr <= 0 for wide E), then the
# decorrelated Delsarte bound is already valid for the true cover (no correlated LP
# needed).  Probe: sign of corr = measS7(E) - decorrelated-plateau across many wide E.
# ============================================================================
def decorrelated_plateau_general(E):
    """P(cover) if each e's sector were independent uniform given the x-marginal.
       We approximate by the 'each element lands uniform & independent' coupon value,
       but element sectors are NOT uniform for small e.  Instead use the exact
       per-element sector law and treat elements as independent: tensor the marginals.
       This is the natural 'decorrelated' surrogate."""
    E = [int(e) for e in E if int(e) != 0]
    # per-element marginal law of sector(e x): for e, sector = floor(7 e x); as x~U,
    # 7 e x mod 7 ... actually sector(e x)=floor(7 frac(e x)); frac(e x)~U[0,1) for e>=1,
    # so sector(e x) ~ Uniform(Z/7) EXACTLY for every e>=1.  So the decorrelated model
    # is: k independent Uniform(Z/7) -> coupon cover prob with k draws.
    k = len(E)
    return coupon_cover_prob(k)

def angle5_corr_sign(bases):
    rows = []
    for name, E in bases.items():
        p0 = measS7(E)
        dec = decorrelated_plateau_general(E)
        rows.append((name, len(E), p0, dec, p0 - dec))
    return rows


# ============================================================================
# ANGLE 6 -- FOURIER / GAUSS-SUM bound exploiting QR(7).
# ----------------------------------------------------------------------------
# The cover event "all 7 sectors hit" has a Fourier expansion over (Z/7)^? .  Write the
# miss-set indicator via characters of Z/7.  P(sector s empty) involves the product over
# e of (prob e avoids s).  The QR(7) structure (residues {1,2,4}) is special: a Gauss
# sum bound on sum_e omega^{...} controls how concentrated the sectors are.  Probe:
# compute the character-sum / Weyl-sum  W_h(E) = (1/|E|) sum_e e(h * e * x-grid) and
# relate its magnitude to the cover deficit.  Specifically: the cover probability is
# governed by the discrepancy  max_h |sum_e ...|; a small Gauss-sum bound => sectors
# near-equidistributed => cover near coupon value.  We test the QR-structured E
# (E = quadratic residues / a Sidon-like set) for anomalously LOW or HIGH p0.
# ============================================================================
def angle6_fourier(E, nx=4900):
    """Diagnostic Weyl sums on a fine x-grid; report max_h |W_h| as an equidistribution
       proxy and the cover prob.  (nx multiple of 7 for clean sampling.)"""
    Efl = [int(e) for e in E if int(e) != 0]
    # sample x on grid, compute per-h Weyl sum of e(h * sector/7) averaged over e and x
    maxW = 0.0
    for h in range(1, P):
        acc = 0j
        for ix in range(nx):
            x = (ix + Fr(1, 2)) / nx
            for e in Efl:
                s = sector((e * x) % 1)
                acc += cmath.exp(2j * pi * h * s / P)
        maxW = max(maxW, abs(acc) / (nx * len(Efl)))
    p0 = measS7(E)
    return p0, maxW


# ============================================================================
def main():
    print("#" * 84)
    print("# THREAD D -- DIVERSE NEW-ANGLE PROBES for the LRC(14) cover bound")
    print("#" * 84)

    # Worst-case / representative WIDE bases (balanced shapes that approach the cap).
    bases = {
        "k=8 even AP [2,4,6,8,10,12,14,16]": [2,4,6,8,10,12,14,16],
        "k=9 even AP [2..18]":               [2,4,6,8,10,12,14,16,18],
        "k=8 consec [1..8]":                 [1,2,3,4,5,6,7,8],
        "k=9 consec [1..9]":                 [1,2,3,4,5,6,7,8,9],
        "k=10 consec [1..10]":               list(range(1,11)),
        "k=9 balanced 2-cluster":            [1,2,3,4, 50,51,52,53,54],  # two far clusters
        "k=8 QR(7)-flavored [1,2,3,4,5,6,7,9]": [1,2,3,4,5,6,7,9],
    }

    # ---- ANGLE 1 ----
    print("\n" + "="*84)
    print("ANGLE 1 -- MATROID/TUTTE of cover incidence (rank distribution per x-cell)")
    print("="*84)
    for name, E in bases.items():
        p0, rd = angle1_matroid(E)
        rdf = {r: round(float(v),4) for r,v in sorted(rd.items())}
        cap = CAP.get(len(E))
        print(f"  {name:42s} p0={float(p0):.5f}"
              f"{' cap='+str(round(float(cap),4)) if cap else ''}")
        print(f"      rank(=#sectors hit) dist: {rdf}")
    print("  NOTE: 'rank=7 <=> cover' is TRIVIAL for singleton transversal matroid; the")
    print("        rank distribution = #sectors-hit distribution.  Tutte adds no new")
    print("        structure beyond the coupon/cover count.  -> see rating in summary.")

    # ---- ANGLE 2 ----
    print("\n" + "="*84)
    print("ANGLE 2 -- TOURNAMENT-SIDE: #sector-classes (chromatic-like) & collisions")
    print("="*84)
    for name, E in bases.items():
        p0, Ecl, Eco, cd = angle2_tournament(E)
        print(f"  {name:42s} p0={float(p0):.5f}  E[#classes]={float(Ecl):.4f}"
              f"  E[#collisions]={float(Eco):.4f}")

    # ---- ANGLE 3 ----
    print("\n" + "="*84)
    print("ANGLE 3 -- DEATH-CHAIN / coupon cover bound (iid uniform sectors)")
    print("="*84)
    print("  iid coupon-collector cover prob C(k) (decorrelated, exact):")
    for k, c in angle3_deathchain(12):
        cap = CAP.get(k)
        tag = f"  cap_{k}={float(cap):.5f}  C(k){'<=' if c<=cap else '>'}cap" if cap else ""
        print(f"    k={k:2d}: C(k)={float(c):.6f}{tag}")

    # ---- ANGLE 4 ----
    print("\n" + "="*84)
    print("ANGLE 4 -- KOKSMA-HLAWKA / BONFERRONI upper bound on p0 (exact arcs)")
    print("="*84)
    print(f"  {'base':42s}{'p0_true':>9}{'B1(up)':>9}{'B2(up)':>9}{'B3(lo)':>9}{'cap':>8}")
    for name, E in bases.items():
        p0, u0, u2, l2, S1, S2, S3 = angle4_bonferroni(E)
        cap = CAP.get(len(E))
        capf = f"{float(cap):.4f}" if cap else "  -  "
        beat = "<=cap!" if (cap and u2 <= cap) else ""
        print(f"  {name:42s}{float(p0):>9.5f}{float(u0):>9.4f}{float(u2):>9.4f}"
              f"{float(l2):>9.4f}{capf:>8}  {beat}")

    # ---- ANGLE 5 ----
    print("\n" + "="*84)
    print("ANGLE 5 -- correlation SIGN: p0_true - decorrelated coupon plateau")
    print("="*84)
    for name, k, p0, dec, corr in angle5_corr_sign(bases):
        cap = CAP.get(k)
        print(f"  {name:42s} k={k} p0={float(p0):.5f} coupon={float(dec):.5f}"
              f"  corr={float(corr):+.5f}  {'(corr<0)' if corr<0 else '(corr>=0!)'}")

    # ---- ANGLE 6 ----
    print("\n" + "="*84)
    print("ANGLE 6 -- FOURIER/Gauss-sum equidistribution proxy (max_h |Weyl|) + QR(7)")
    print("="*84)
    for name, E in bases.items():
        p0, mw = angle6_fourier(E, nx=2100)
        print(f"  {name:42s} p0={float(p0):.5f}  max_h|Weyl|={mw:.5f}")

    print("\nDONE.")

if __name__ == "__main__":
    main()
