#!/usr/bin/env python3
"""
LRC(14) S3 residual -- ANGLE = DIRECT COVERING-SYSTEM BOUND (kind-pasteur-2026-06-18-S5-wf).

GOAL: attempt a SINGLE uniform lower bound  mu(E) >= c0 > 0  with NO spread split,
via the 0-anchored gap / strip-avoidance measure.

SETUP (recap, all upstream-PROVED THM-527/528):
  E = co-offset set, 0 in E, |E| = k <= 13, gcd(E)=1 WLOG (L1).
  mu(E) = meas{ x in [0,1): circular max-gap of {frac(e x): e in E} > 2/7 }.
  Uniform floor inf_E mu(E) >= c0 > 0  ==>  (intersect G_P) ==> LRC(14).

THE 0-ANCHORED GAP IDEA.
  Since 0 in E, the point frac(0*x)=0 is ALWAYS present. Consider the fixed open arc
  A = (1/7, 3/7) of length 2/7. The endpoint 0 of A's complement... let's be precise.
  Define
     mu_0(E) := meas{ x : frac(e x) NOTIN (1/7, 3/7)  for ALL e in E }.
  If every frac(e x) avoids the OPEN arc (1/7,3/7), then that whole open arc of length
  2/7 lies inside ONE gap between consecutive points (the arc contains no point), so the
  circular max-gap is >= 2/7. For the STRICT inequality > 2/7 we instead use the CLOSED
  arc [1/7, 3/7]: if all points avoid the closed arc [1/7,3/7] (length 2/7), the open
  super-arc strictly containing it of some length > 2/7 is point-free, giving max-gap
  > 2/7 ... actually avoiding closed [1/7,3/7] only guarantees max-gap >= 2/7 still,
  because points could sit AT 1/7 and 3/7. We handle the strictness carefully below.

  KEY CONTAINMENT (the lemma we use):
     { x : frac(e x) in [0,1/7) U (3/7, 1)  for all e } subset { maxgap > 2/7 }.
  Reason: all points lie in the closed complementary arc B = [3/7,1] U [0,1/7] (mod 1)
  which is the arc from 3/7 forward to 1/7, of length 1/7 + (1 - 3/7) = 1/7+4/7 = 5/7.
  Wait recompute: complement of (1/7,3/7) is [3/7,1) U [0,1/7] going the long way; its
  length is 1 - 2/7 = 5/7. If all k points lie in a CLOSED arc of length 5/7, the
  remaining open arc of length 2/7 is point-free => some gap >= 2/7. For STRICT > 2/7,
  require all points in the OPEN arc (3/7, 1/7+1) i.e. strictly avoid CLOSED [1/7,3/7].
  Then the point-free open arc has length strictly > 2/7. Good. So define

     mu_0(E) := meas{ x : frac(e x) in (3/7, 1) U {0} U (0, 1/7)  ... }  -- messy at 0.

  Cleaner: since 0 in E forces the point 0 to be present, and 0 is NOT in [1/7,3/7],
  use the OPEN-avoidance:
     S0(E) := { x : for all e in E, frac(e x) NOTIN [1/7, 3/7] }   (closed danger arc)
  On S0(E), all points avoid closed [1/7,3/7], so the open arc (1/7,3/7) (length 2/7,
  but actually a bit more is point free) ... The point-free region is at least the open
  (1/7,3/7); whether maxgap is STRICTLY > 2/7 needs the points nearest 1/7 and 3/7 to be
  strictly outside. On S0 they are (closed danger excluded => points are < 1/7 or > 3/7
  strictly OR equal to boundary? closed danger = [1/7,3/7] EXCLUDED means frac < 1/7 or
  frac > 3/7 strictly). So nearest points are strictly outside => point-free open arc
  strictly contains [1/7,3/7]?? No: nearest point below could be 1/7 - eps, nearest above
  3/7 + eps, point-free open arc = (1/7-eps, 3/7+eps) length 2/7 + 2eps > 2/7. STRICT.

  THEREFORE:  S0(E) subset {maxgap > 2/7},  so  mu(E) >= meas(S0(E)) =: mu_0(E).
  And mu_0(E) = meas{ x : ||e x|| ... }  no -- it's  frac(e x) in [0,1/7) U (3/7,1).
  This is a k-fold STRIP-AVOIDANCE measure. We study its uniform floor.

  NOTE the danger arc [1/7,3/7] is NOT symmetric about 0; this is the "0-anchored" arc.
  By reflection L0 (x<->1-x) we could also use [4/7,6/7]; we may take the MAX or AVERAGE.

This script:
  (A) EXACT mu(E) tool, RIGOROUS (adds gap=2/7 breakpoints).
  (B) EXACT mu_0(E) tool for the fixed danger arc [1/7,3/7].
  (C) Verify mu_0(E) <= mu(E) (containment) on many shapes.
  (D) Search for the INFIMUM of mu_0(E): is it > 0 uniformly? Find the worst shapes.
  (E) The density/covering count: sum of strip densities = k * (5/7) "good" vs cover.
  (F) Honest verdict on whether a single uniform mu_0 >= c0 closes B(k).
"""
from fractions import Fraction as F
from math import gcd, comb
from itertools import combinations, product
import random

G0 = F(2, 7)          # max-gap threshold
DLO, DHI = F(1, 7), F(3, 7)   # 0-anchored closed danger arc [1/7,3/7]


# ===========================================================================
# (A) RIGOROUS EXACT mu(E)
#   x -> frac(e x) piecewise linear; order changes at x = m/(e_i-e_j).
#   Within an order cell, the sorted points are fixed-order linear functions, each gap
#   is linear in x, so {maxgap > 2/7} on the cell is a union of sub-intervals bounded by
#   the order breakpoints AND the points where some gap = 2/7 (also rational, linear solve).
#   We compute exactly by collecting ALL breakpoints (collisions + gap=2/7 crossings),
#   then on each atomic interval evaluate maxgap at the midpoint (constant order, and
#   maxgap-vs-2/7 sign constant since we've added all =2/7 crossings).
# ===========================================================================
def _frac(q):
    return q - (q.__floor__())

def _all_breakpoints(E):
    E = sorted(set(E))
    k = len(E)
    bp = set([F(0), F(1)])
    # collisions: (e_i - e_j) x in Z
    for i in range(k):
        for j in range(i + 1, k):
            d = E[j] - E[i]
            if d == 0:
                continue
            for m in range(0, d + 1):
                bp.add(F(m, d))
    # gap = 2/7 crossings: for an ordered pair (e_a, e_b), the (signed) circular gap
    #   frac(e_b x) - frac(e_a x) (mod 1) = 2/7. Since within a cell frac(e x)=e x - floor,
    #   the difference is (e_b-e_a) x - integer. So (e_b-e_a) x = integer + 2/7.
    #   => x = (n + 2/7)/(e_b - e_a) for integer n. Also the WRAP gap uses 1 - sum, but
    #   a single gap equal 2/7 is captured by some pair difference = 2/7 mod 1.
    #   We add x = (n + 2/7)/D and (n - 2/7)/D = (n+5/7)/D for all pair diffs D and n in range.
    diffs = set()
    for i in range(k):
        for j in range(k):
            if i != j:
                diffs.add(E[j] - E[i])
    for D in diffs:
        if D == 0:
            continue
        aD = abs(D)
        # x in [0,1): (n + 2/7)/aD in [0,1) => n in [-2/7, aD-2/7]
        n = 0
        while True:
            x1 = (F(n) + F(2, 7)) / aD
            x2 = (F(n) + F(5, 7)) / aD
            added = False
            if F(0) <= x1 < F(1):
                bp.add(x1); added = True
            if F(0) <= x2 < F(1):
                bp.add(x2); added = True
            if x1 >= F(1) and x2 >= F(1):
                break
            n += 1
            if n > aD + 2:
                break
    return sorted(b for b in bp if F(0) <= b < F(1))

def maxgap_at(E, x):
    pts = sorted(set(_frac(e * x) for e in E))
    if len(pts) == 1:
        return F(1)
    gaps = [pts[t + 1] - pts[t] for t in range(len(pts) - 1)]
    gaps.append(pts[0] + 1 - pts[-1])
    return max(gaps)

def mu_exact(E):
    bp = _all_breakpoints(E)
    pts_bp = bp + [F(1)]
    tot = F(0)
    for a, b in zip(bp, pts_bp[1:]):
        if b <= a:
            continue
        mid = (a + b) / 2
        if maxgap_at(E, mid) > G0:
            tot += (b - a)
    return tot


# ===========================================================================
# (B) EXACT mu_0(E): meas{ x : frac(e x) NOTIN [1/7,3/7] for all e in E }.
#   = meas{ x : for all e, frac(e x) in [0,1/7) U (3/7,1) }.
#   For each e, the BAD set B_e = {x: frac(e x) in [1/7,3/7]} is a union of e intervals
#   each of length (2/7)/e: x in [ (n+1/7)/e, (n+3/7)/e ] for n=0..e-1.
#   GOOD_e = complement. mu_0 = meas( intersection_e GOOD_e ).
#   We compute exactly via breakpoints = all (n+1/7)/e and (n+3/7)/e.
# ===========================================================================
def mu0_exact(E, dlo=DLO, dhi=DHI):
    E = sorted(set(E))
    bp = set([F(0), F(1)])
    for e in E:
        if e == 0:
            continue
        for n in range(e + 1):
            for off in (dlo, dhi):
                x = (F(n) + off) / e
                if F(0) <= x < F(1):
                    bp.add(x)
    bp = sorted(bp)
    tot = F(0)
    for a, b in zip(bp, bp[1:]):
        if b <= a:
            continue
        mid = (a + b) / 2
        ok = True
        for e in E:
            if e == 0:
                continue
            fr = _frac(e * mid)
            if dlo <= fr <= dhi:
                ok = False
                break
        if ok:
            tot += (b - a)
    return tot


# Reflected danger arc [4/7,6/7] (0-anchored on the other side), and we may take the
# union-good = avoid BOTH? No: avoiding a single arc suffices for one gap. Using the arc
# that is "most avoidable". We also compute mu_0 with arc centered to MAXIMIZE measure.
def mu0_best(E):
    """max over the two natural 0-anchored danger arcs [1/7,3/7] and [4/7,6/7]."""
    a = mu0_exact(E, F(1, 7), F(3, 7))
    # [4/7,6/7]:
    b = mu0_exact(E, F(4, 7), F(6, 7))
    return max(a, b), a, b


# ===========================================================================
# Helpers
# ===========================================================================
def normalize(E):
    E = sorted(set(E))
    g = 0
    for e in E:
        g = gcd(g, e)
    if g == 0:
        return E
    return [e // g for e in E]

def header(t):
    print("\n" + "=" * 74)
    print(t)
    print("=" * 74)


# ===========================================================================
# MAIN
# ===========================================================================
if __name__ == "__main__":
    header("(C) CONTAINMENT CHECK: mu_0(E) <= mu(E) on many shapes")
    random.seed(1)
    bad_contain = 0
    tested = 0
    examples = []
    # consecutive, perforated, AP, random bounded
    shape_gens = []
    for k in range(3, 10):
        shape_gens.append(list(range(k)))                      # consecutive 0..k-1
    for k in range(3, 9):
        shape_gens.append([0] + [2 * i + 1 for i in range(k - 1)])  # 0 + odds
    for _ in range(200):
        k = random.randint(3, 11)
        hi = random.randint(k, 3 * k + 5)
        E = normalize([0] + random.sample(range(1, hi + 1), k - 1))
        shape_gens.append(E)
    worst0 = None
    for E in shape_gens:
        E = normalize(E)
        if len(E) < 3:
            continue
        m0 = mu0_exact(E)
        m = mu_exact(E)
        tested += 1
        if m0 > m:
            bad_contain += 1
            if len(examples) < 5:
                examples.append((E, m0, m))
        if worst0 is None or m0 < worst0[1]:
            worst0 = (E, m0, m)
    print(f"tested {tested} shapes; containment violations (mu_0 > mu): {bad_contain}")
    for ex in examples:
        print("  VIOLATION", ex)
    print(f"worst (smallest) mu_0 seen: E={worst0[0]} mu_0={worst0[1]}={float(worst0[1]):.5f} mu={worst0[2]}={float(worst0[2]):.5f}")

    header("(D) INFIMUM SEARCH for mu_0(E): is the 0-anchored floor positive & uniform?")
    # systematic: consecutive 0..k-1 for k up to 13, and the known min-mu shapes
    print("Consecutive E = {0,1,...,k-1}:")
    for k in range(3, 14):
        E = list(range(k))
        m0, ma, mb = mu0_best(E)
        m = mu_exact(E)
        print(f"  k={k:2d}: mu_0(arc1)={float(mu0_exact(E)):.5f}={mu0_exact(E)}  "
              f"mu_0(best)={float(m0):.5f}  mu={float(m):.5f}={m}")

    header("(D2) WORST-CASE mu_0 over broad random + structured scan (find inf)")
    random.seed(7)
    overall_worst = None
    overall_worst_best = None
    n_scan = 4000
    for _ in range(n_scan):
        k = random.randint(3, 13)
        hi = random.randint(k, 4 * k + 10)
        pool = random.sample(range(1, hi + 1), k - 1)
        E = normalize([0] + pool)
        if len(E) < 3:
            continue
        m0 = mu0_exact(E)
        mbest, ma, mb = mu0_best(E)
        if overall_worst is None or m0 < overall_worst[1]:
            overall_worst = (E, m0)
        if overall_worst_best is None or mbest < overall_worst_best[1]:
            overall_worst_best = (E, mbest)
    print(f"scanned {n_scan} random bounded shapes (k<=13)")
    print(f"  worst mu_0(arc1) : E={overall_worst[0]}")
    print(f"                     mu_0={overall_worst[1]}={float(overall_worst[1]):.6f}")
    print(f"  worst mu_0(best) : E={overall_worst_best[0]}")
    print(f"                     mu_0={overall_worst_best[1]}={float(overall_worst_best[1]):.6f}")

    header("(E) THE COVERING / DENSITY ACCOUNTING")
    # Each strip B_e = {frac(e x) in [1/7,3/7]} has density 2/7. Good_e density 5/7.
    # If the e were 'independent', mu_0 ~ (5/7)^k -> 0 as k grows. So a NAIVE
    # independence heuristic predicts mu_0 -> 0. The question: does ARITHMETIC STRUCTURE
    # (the e_i being integer multiples sharing x) keep mu_0 bounded below?
    for k in range(3, 14):
        print(f"  k={k:2d}: (5/7)^k = {float((F(5,7))**k):.6f}  [naive independence prediction]")
    print("  -> If mu_0 tracked (5/7)^k it would VANISH. Check actual vs naive:")
    for k in range(3, 14):
        E = list(range(k))
        print(f"     k={k:2d}: actual mu_0(consec)={float(mu0_exact(E)):.6f}  naive={float((F(5,7))**k):.6f}")
