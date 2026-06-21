#!/usr/bin/env python3
"""
lrc14_angleE_cumulant_telescope_macmini_0620s8.py  (mac-mini-2026-06-20-S8)  ANGLE E part 2

Linear regroupings of the 64 signed J-terms all FAIL (part 1).  The failures are LARGE
(+0.19..+0.48), so no coarsening of the sign-weighted sum is consec-extremal.

NEW IDEA.  The signed sum measS7 = Sum_M (-1)^|M| J(M) is the inclusion-exclusion COVER
probability.  Equivalently it is a value of the COLOR-COUNT generating function.  Define the
"occupancy generating polynomial".  For a clock vector, each x in [0,1) realises a SUBSET
C(x) subset Z/7 of HIT residues (always 0 in C(x)).  Then
   J(M) = meas{ C(x) cap M = empty } = meas{ C(x) subset Z/7 \ M },
   measS7 = meas{ C(x) = Z/7 }.
Let  p_S := meas{ C(x) = S }  for S subset Z/7 with 0 in S.  Then
   J(M) = Sum_{S: S cap M = empty} p_S,     measS7 = p_{Z/7}.
The p_S are the "color-pattern occupation measures" (the Potts microstate weights).
measS7 = p_full is a SINGLE microstate weight -- the all-colors-present pattern.

So the cleanest object is the VECTOR (p_S)_{0 in S}.  We test:
  (T1) Is consec's p_full the max?  (definition -- yes by crux, recheck)
  (T2) MAJORIZATION: does consec's pattern vector majorize / dominate adversaries in a way
       that forces p_full max?  Test: among patterns S with |S|>=6, does consec dominate?
  (T3) the "deficiency" telescoping: p_full = 1 - Sum_{S != full} p_S, and group the missing
       patterns by which/how many residues are absent.  Missing-1 patterns (exactly one color
       absent) are the dominant deficiency.  measS7 = 1 - Sum_j q_j + (higher), where
       q_j = meas{exactly residue j absent}?  No -- use the cleaner CHAIN:
       1 - measS7 = meas{ some residue absent } = meas{ C(x) != full }.
       Define A_j = {residue j absent}.  1 - measS7 = meas(Union A_j).  We want to MINIMIZE
       meas(Union A_j) over E (consec should minimize the uncovered measure).
       Bonferroni:  meas(Union) >= Sum meas(A_j) - Sum meas(A_i cap A_j)  (lower bd, wrong dir)
                    meas(Union) <= Sum meas(A_j)                            (upper bd, right dir for min? no)
       We want consec to MINIMIZE meas(Union A_j).  Test whether consec minimizes Sum_j meas(A_j)
       (the level-1 Bonferroni) -- and whether the CORRECTION terms are controlled.
  (T4) the absent-residue COUNT distribution: let N(x)=7-|C(x)| = #absent colors.  measS7=P(N=0).
       The full distribution of N is a 8-vector (N=0..7).  Compare consec vs adversaries.
"""
import sys, itertools
from fractions import Fraction as F
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

def hit_pattern_measures(E):
    """return dict S(frozenset, 0 in S) -> p_S = meas{C(x)=S}."""
    E = sorted(set(E)); bps = {F(0), F(1)}
    for e in E:
        if e == 0: continue
        for m in range(0, 7 * e + 1): bps.add(F(m, 7 * e))
    bps = sorted(bps); pat = defaultdict(lambda: F(0))
    for i in range(len(bps) - 1):
        x0, x1 = bps[i], bps[i + 1]
        if x1 <= x0: continue
        xm = (x0 + x1) / 2
        S = frozenset(int(((e * xm) % 1) * 7) for e in E)
        pat[S] += x1 - x0
    return pat

def absent_count_dist(pat):
    """N(x)=7-|C(x)|; return list d[0..7]=meas{N=c}."""
    d = [F(0)] * 8
    for S, p in pat.items():
        d[7 - len(S)] += p
    return d

def absent_residue_meas(pat):
    """A_j = meas{residue j absent} for j=0..6 (A_0=0 since e=0 pins residue0)."""
    A = [F(0)] * 7
    for S, p in pat.items():
        for j in range(7):
            if j not in S: A[j] += p
    return A

def pair_absent_meas(pat):
    """A_{ij}=meas{residues i,j both absent}."""
    A = {}
    for i, j in itertools.combinations(range(7), 2):
        A[(i, j)] = F(0)
    for S, p in pat.items():
        for i, j in itertools.combinations(range(7), 2):
            if i not in S and j not in S: A[(i, j)] += p
    return A

def measS7(pat):
    return pat[frozenset(range(7))]

def main():
    k = 8
    consec = list(range(k))
    pc = hit_pattern_measures(consec)
    s7c = measS7(pc)
    print("=" * 96)
    print(f"ANGLE E part 2: occupation-measure / absent-count telescoping, k={k}")
    print("=" * 96)
    print(f"\nmeasS7(consec)={s7c}={float(s7c):.6f}")

    dc = absent_count_dist(pc)
    print("\n(T4) absent-color count distribution N=0..7 for consec:")
    print("   " + ", ".join(f"N={c}:{float(dc[c]):.4f}" for c in range(8)))

    Ac = absent_residue_meas(pc)
    print("\n(T3) per-residue absent measure A_j (consec), want consec to MINIMIZE these:")
    print("   " + ", ".join(f"A_{j}={float(Ac[j]):.4f}" for j in range(7)))
    print(f"   Sum_j A_j (level-1 Bonferroni for uncovered) = {float(sum(Ac)):.5f}")
    print(f"   1-measS7 (true uncovered) = {float(1-s7c):.5f}")

    # adversary pool
    N = 11
    pool = [(0,) + c for c in itertools.combinations(range(1, N + 1), k - 1)]
    print(f"\nadversary pool (0..{N}): {len(pool)} sets")

    # T3a: does consec MINIMIZE Sum_j A_j ?
    sumA_c = sum(Ac)
    beat_sumA = 0; worst = None
    # T3b: does consec MINIMIZE the true uncovered meas (= 1 - measS7)? (this IS the crux, restated)
    beat_true = 0
    # T4a: does consec STOCHASTICALLY DOMINATE adversaries on N (P(N<=c) >= adv for all c)?
    # i.e. consec's absent-count is "smaller" -> more covered. Test prefix sums.
    cdf_c = [sum(dc[:c+1]) for c in range(8)]  # P(N<=c)
    dominates_cdf = 0; cdf_fail_detail = None
    for E in pool:
        pe = hit_pattern_measures(E)
        Ae = absent_residue_meas(pe)
        if sum(Ae) < sumA_c - F(0):
            beat_sumA += 1
            if worst is None or sum(Ae) < worst[1]:
                worst = (E, sum(Ae))
        if measS7(pe) > s7c:
            beat_true += 1
        de = absent_count_dist(pe)
        cdf_e = [sum(de[:c+1]) for c in range(8)]
        # consec stochastically dominates (more covered) iff cdf_c[c] >= cdf_e[c] for all c
        if all(cdf_c[c] >= cdf_e[c] for c in range(8)):
            pass
        else:
            dominates_cdf += 1
            if cdf_fail_detail is None:
                # find first c where adv has more mass at low N
                for c in range(8):
                    if cdf_e[c] > cdf_c[c]:
                        cdf_fail_detail = (E, c, cdf_e[c] - cdf_c[c]); break

    print(f"\n(T3a) adversaries with SMALLER Sum_j A_j than consec: {beat_sumA}/{len(pool)}")
    if worst: print(f"      most-beating: E={worst[0]} SumA={float(worst[1]):.5f} (consec {float(sumA_c):.5f})")
    print(f"(crux) adversaries with LARGER measS7 than consec (should be 0): {beat_true}/{len(pool)}")
    print(f"\n(T4) adversaries NOT stochastically dominated by consec on absent-count N: {dominates_cdf}/{len(pool)}")
    if cdf_fail_detail:
        E, c, d = cdf_fail_detail
        print(f"      e.g. E={E} has more mass at N<={c} by {float(d):+.5f}  (so consec does NOT stoch-dominate)")
    else:
        print(f"      *** consec STOCHASTICALLY DOMINATES every adversary on N! (P(N<=c) larger for all c) ***")

if __name__ == "__main__":
    main()
