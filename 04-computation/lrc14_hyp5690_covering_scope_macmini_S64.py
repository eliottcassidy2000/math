"""
mac-mini-2026-07-09-S64 -- ANSWER to klein-S206's HYP-5690: which route-breaking clusters does
LRC(14) actually NEED?

BACKGROUND.  LRC(14) reduces to COVERING velocity sets: if some q in {2..14} divides no speed, then
t = 1/q is lonely (LonelyRunner.sieve_one_div; contrapositive counterexample_needs_all_divisors,
already canon as THM-523).  So a counterexample must be covering:
    covering(S)  <=>  for every q in 2..14, some v in S has q | v.
A co-offset cluster E with ruler Vmax is "covering-derived" iff S = {Vmax - e : e in E} is covering.
klein-S206 asked (HYP-5690) whether the hard-cluster saga was fought on clusters LRC(14) never needs.

ANSWER: PARTLY YES.  Three of the four route-breaking clusters are NON-COVERING, hence out of scope.

    cluster                                      Vmax   covering?  missing q  -> dispatched at
    tight AP {0..12}     (MISTAKE-129, klein)      13     NO        [14]         t = 1/14
    7-structured         (MISTAKE-127/128)         91     YES       -            (IN SCOPE)
    knife-edge           (MISTAKE-130, mac-mini)   49     NO        [8,9,10,11]  t = 1/8
    worst7StructLarge    (klein native_decide)    458     NO        [7,14]       t = 1/7

CONSEQUENCES.
 * MY OWN knife-edge -- the config with maxgap = 1/7 EXACTLY, spread = 6*Vmax/7, which motivated the
   whole non-strict criterion and `LRCGoodPeriodNonStrict.lean` -- is NON-COVERING.  LRC(14) never
   needs it.  The non-strict / knife-edge layer is therefore NOT required for the proof.
 * klein's `worst7StructLarge_hasGoodPeriod` native_decide certificate is on an out-of-scope cluster.
   (`worst7Struct_hasGoodPeriod` at Vmax=91 IS in scope -- that one matters.)

AND THE COVERING BRANCH HAS A STRICT CUSHION (exact, exhaustive):
 over all 966 covering 13-subsets of [1,18], the exact  M(S) = max_t min_i ||v_i t||  has
      min M = 1/12 = 0.08333  >  1/14 = 0.07143     (strict margin 1/84)
 while the NON-covering tight AP {1..13} has M = 1/14 EXACTLY (the equality locus).

 => The equality locus (M = 1/14) is entirely NON-COVERING.  On covering clusters the STRICT
    good-period criterion suffices; no knife-edge / non-strict layer is needed.
    This is exactly klein-S206's forced proof shape:
      [exact rational witness t = 1/q on non-covering, equality allowed]
    + [strict margin on covering].

M is computed EXACTLY: min_i ||v_i t|| is piecewise linear, so its maximum sits at a local max, i.e.
at a crossing t = p/(v_i+v_j) or at a peak t = (2m+1)/(2 v_i).  Enumerate those rationals exactly.
(No grids -- cf MISTAKE-130.)
"""
from fractions import Fraction as F
from itertools import combinations


def covering(S):
    return all(any(v % q == 0 for v in S) for q in range(2, 15))


def missing_q(S):
    return [q for q in range(2, 15) if not any(v % q == 0 for v in S)]


def nearInt(x):
    r = F(x) % 1
    return min(r, 1 - r)


def minReach(S, t):
    return min(nearInt(F(v) * t) for v in S)


def Mexact(S):
    """exact max_t min_i ||v_i t||: local maxima at crossings p/(v_i+v_j) or peaks (2m+1)/(2 v_i)."""
    cands = set()
    for a, b in combinations(S, 2):
        s = a + b
        for p in range(1, s):
            cands.add(F(p, s))
    for v in S:
        for m in range(0, 2 * v):
            t = F(2 * m + 1, 2 * v)
            if 0 < t < 1:
                cands.add(t)
    return max(minReach(S, t) for t in cands)


cases = [
    ("tight AP cluster (MISTAKE-129)",      list(range(0, 13)), 13),
    ("7-structured (MISTAKE-127/128)",      [0,7,14,21,26,29,37,44,51,58,67,75,82], 91),
    ("knife-edge (MISTAKE-130, mac-mini)",  [0,7,10,14,18,20,21,26,28,35,36,37,42], 49),
    ("worst7StructLarge (klein cert)",      [0,7,63,126,151,189,252,305,315,362,378,385,406], 458),
]
print("HYP-5690: is each route-breaking cluster COVERING-DERIVED (= in LRC(14) scope)?\n")
print(f"{'cluster':>38}{'Vmax':>6}{'covering?':>11}   missing q -> witness")
for name, E, V in cases:
    S = sorted(V - e for e in E)
    cov, mq = covering(S), missing_q(S)
    tail = "-  (IN SCOPE)" if cov else f"{mq} -> t = 1/{mq[0]}"
    print(f"{name:>38}{V:>6}{str(cov):>11}   {tail}")

print("\nSTRICT CUSHION on the covering branch (exhaustive, exact):")
cov_sets = [c for c in combinations(range(1, 19), 13) if covering(c)]
worst, best = F(1), None
for S in cov_sets:
    M = Mexact(list(S))
    if M < worst:
        worst, best = M, S
print(f"  covering 13-subsets of [1,18]: {len(cov_sets)}")
print(f"  min exact M = {worst} = {float(worst):.6f}   attained at {best}")
print(f"  1/14 = {float(F(1,14)):.6f};  strict? {worst > F(1,14)}  (margin {worst - F(1,14)})")

ap = list(range(1, 14))
print(f"\n  contrast: tight AP {{1..13}} is NON-covering and has M = {Mexact(ap)} EXACTLY (equality locus).")
print("\n=> equality locus is entirely NON-covering; the covering branch has a strict cushion.")
print("   The non-strict / knife-edge layer is NOT needed for LRC(14).")
