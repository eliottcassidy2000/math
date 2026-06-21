#!/usr/bin/env python3
"""
ROUTE 5, part 5 -- DOES the tournament OCF level-sum mechanism transfer to
the LRC measS7?

TOURNAMENT (mechanism PROVED at p=7,11, structural p=19,23):
  H = sum_j 2^j alpha_j  (independence poly of odd-cycle conflict graph at
  fugacity x=2). Level 1 (cycle counts) favors Paley/flat; levels j>=2
  (disjoint packings) favor the AP/interval; the 2^j growth tips it to the
  AP for large p. The "c_k flip" = the level where the AP overtakes.

LRC:  measS7 = P(all 6 inner sectors hit) over the parameter y in [-1/14,1/14]
  pooled over resonances a=1..6. By INCLUSION-EXCLUSION on the 6 hit-events:
     measS7 = sum_{S subset of [6]} (-1)^|S| * MISS(S),
  where MISS(S) = (pooled) measure that ALL sectors in S are simultaneously
  empty (missed). This is the EXACT LRC analogue of the OCF: a signed sum
  over subsets, weighted (here by (-1)^|S| rather than 2^|S|).

This script computes, EXACTLY (Fractions), per resonance a and pooled:
  - the joint occupancy law of the 6 inner sectors (sectors 1..6; sector 0
    is the e=0 perfect short, always occupied),
  - MISS_k(E) = sum_{|S|=k} MISS(S) = the level-k miss mass,
  - the inclusion-exclusion reconstruction measS7 = sum_k (-1)^k MISS_k,
  - the Walsh coefficients qhat_S and verifies measS7 = (1/64) sum_S qhat_S.

Then it tests the TRANSFER QUESTION across full-residue k=8 shapes:
  - Is consec the AP-extremal corner of the level masses MISS_k?
  - Which level k carries the consec advantage (the analogue of the
    tournament's "first flip" level)?
  - Is there a c_k * e_k(spectral-moment) form where consec is the corner?

Author: mac-mini-2026-06-21-S20 (ROUTE 5)
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
from functools import reduce
from collections import defaultdict
sys.path.insert(0, '04-computation')
from lrc14_route1_conductance_minimax_opus_0621 import (
    sec, breakpoints, is_full_residue, primitive, cond_vec, consec, measS7,
)
sys.stdout.reconfigure(line_buffering=True)

HALF = F(1, 14)


def occ_set(E, a, y):
    """Set of occupied sectors (0..6) at (a,y)."""
    return set(sec(e, a, y) for e in E)


def miss_level_masses(E):
    """Compute, pooled over a=1..6, the measure of {exactly the sectors in S
    among 1..6 are MISSED (empty), the rest hit}, aggregated to level masses.

    We want MISS(S) = measure that ALL sectors in S are empty (regardless of
    others). Then inclusion-exclusion:
      measS7 = P(all 1..6 hit) = sum_S (-1)^|S| MISS(S).
    Sector 0 is the e=0 short (always occupied) -- but careful: full cover
    needs ALL 7 sectors incl 0. With e=0 in E, sector 0 always occupied, so
    'all 7 hit' == 'all of 1..6 hit'. We verify this matches measS7.
    """
    # MISS(S) for S subset of {1..6}: measure where sectors in S all empty.
    miss = defaultdict(F)  # frozenset(S) -> measure
    for a in range(1, 7):
        bps = breakpoints(E, a, HALF)
        for lo, hi in zip(bps, bps[1:]):
            if hi <= lo:
                continue
            mid = (lo + hi) / 2
            occ = occ_set(E, a, mid)
            empty = frozenset(s for s in range(1, 7) if s not in occ)
            length = hi - lo
            # this interval contributes to MISS(S) for every S subset of empty
            # but we accumulate by the EXACT empty set; convert to MISS(S) after
            miss[("EXACT", empty)] += length
    # Now MISS(S) = sum over exact-empty E' with S subset of E'
    exact = {k[1]: v for k, v in miss.items()}
    # level masses by exact empty size
    exact_level = defaultdict(F)
    for E_set, v in exact.items():
        exact_level[len(E_set)] += v
    # MISS(S) aggregated to level: MISS_k = sum_{|S|=k} MISS(S)
    # MISS(S) = sum_{E' superset S} exact(E')
    miss_S = defaultdict(F)
    allsec = list(range(1, 7))
    for k in range(0, 7):
        for S in itertools.combinations(allsec, k):
            Sset = frozenset(S)
            tot = F(0)
            for E_set, v in exact.items():
                if Sset <= E_set:
                    tot += v
            miss_S[Sset] = tot
    miss_level = defaultdict(F)
    for Sset, v in miss_S.items():
        miss_level[len(Sset)] += v
    return exact, exact_level, miss_S, miss_level


def reconstruct_measS7(miss_level):
    """measS7 = sum_k (-1)^k MISS_k. (MISS_0 = total pooled measure = 6*(1/7)
    since y ranges over [-1/14,1/14] length 1/7 for each of 6 resonances.)"""
    return sum((-1) ** k * miss_level[k] for k in range(0, 7))


def walsh_coeffs(exact):
    """qhat_S = E[prod_{j in S}(2 h_j - 1)] over the pooled measure, where
    h_j = 1 if sector j hit. Normalize by total mass M0 = 6/7.
    Here 'E' is the pooled-measure average. Returns dict S->qhat, and the
    sum check (1/64) sum_S qhat_S * M0  vs measS7 (with the convention used
    in HYP-2758, measS7 = (1/64) sum_S qhat_S where qhat over the NORMALIZED
    occupancy law)."""
    M0 = F(6, 7)  # total pooled measure (6 resonances x interval length 1/7)
    # build occupancy law: for each exact-empty set, h_j = 1 iff j not in empty
    qhat = {}
    allsec = list(range(1, 7))
    for k in range(0, 7):
        for S in itertools.combinations(allsec, k):
            Sset = frozenset(S)
            acc = F(0)
            for E_set, v in exact.items():
                pr = 1
                for j in Sset:
                    hj = 0 if j in E_set else 1
                    pr *= (2 * hj - 1)
                acc += pr * v
            qhat[Sset] = acc / M0  # normalized expectation
    return qhat, M0


def main():
    print("#" * 78)
    print("# ROUTE 5 part 5: LRC inclusion-exclusion level-sum vs OCF transfer")
    print("#" * 78)

    C = consec(8)
    print(f"\nconsec_8 = {C}")
    ms = measS7(C)
    print(f"measS7(consec_8) = {ms} = {float(ms):.6f}")

    exact, exact_level, miss_S, miss_level = miss_level_masses(C)
    print("\nEXACT-empty level masses (measure with EXACTLY k sectors empty):")
    for k in sorted(exact_level):
        print(f"  k={k}: {exact_level[k]} = {float(exact_level[k]):.6f}")
    print("  (k=0 = measure fully covered = measS7 contribution)")

    print("\nMISS_k = sum_{|S|=k} P(all of S empty):")
    for k in sorted(miss_level):
        print(f"  MISS_{k} = {miss_level[k]} = {float(miss_level[k]):.6f}")

    recon = reconstruct_measS7(miss_level)
    print(f"\ninclusion-exclusion measS7 = sum_k (-1)^k MISS_k = {recon} = {float(recon):.6f}")
    print(f"matches measS7? {recon == ms}")

    qhat, M0 = walsh_coeffs(exact)
    s_qhat = sum(qhat.values())
    print(f"\nWalsh check: (M0/64) * sum_S qhat_S = {(M0/64)*s_qhat} "
          f"vs measS7={ms}  match={ (M0/64)*s_qhat == ms }")
    nonneg = all(q >= 0 for q in qhat.values())
    print(f"all qhat_S >= 0 (Fourier-positive)? {nonneg}")

    # ---- compare consec vs a few other full-residue k=8 shapes ----
    print("\n--- level-mass comparison: consec vs other full-residue k=8 shapes ---")
    shapes = [
        ("consec", list(range(8))),
        ("0..6,8", [0,1,2,3,4,5,6,8]),
        ("0..6,9", [0,1,2,3,4,5,6,9]),
        ("0..6,13", [0,1,2,3,4,5,6,13]),
        ("0,2,4,6,8,10,12,7", [0,2,4,6,8,10,12,7]),  # AP-ish even+filler
        ("0..5,7,8", [0,1,2,3,4,5,7,8]),
    ]
    print(f"  {'shape':>22} {'measS7':>12} | MISS_1 MISS_2 MISS_3 MISS_4 MISS_5 MISS_6")
    rows = []
    for name, E in shapes:
        if not is_full_residue(E):
            print(f"  {name:>22}  NOT full-residue, skip")
            continue
        _, _, _, mlev = miss_level_masses(E)
        ms_e = measS7(E)
        rows.append((name, ms_e, mlev))
        mvals = " ".join(f"{float(mlev[k]):.4f}" for k in range(1, 7))
        print(f"  {name:>22} {float(ms_e):>12.6f} | {mvals}")

    # Which level discriminates? consec should MINIMIZE the alternating-sum
    # of misses appropriately. Report per-level consec-minus-other.
    print("\n  consec advantage per MISS level (consec - other), "
          "weighted by (-1)^k:")
    base = rows[0][2]
    for name, ms_e, mlev in rows[1:]:
        contribs = [((-1)**k)*(base[k]-mlev[k]) for k in range(1,7)]
        net = sum(contribs)
        print(f"  vs {name:>20}: per-level (-1)^k*(consec-other) = "
              f"{[f'{float(c):+.4f}' for c in contribs]}  net={float(net):+.6f} "
              f"(measS7 diff={float(rows[0][1]-ms_e):+.6f})")


if __name__ == "__main__":
    main()
