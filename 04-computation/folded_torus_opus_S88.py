#!/usr/bin/env python3
"""
opus-2026-07-05-S88 -- HYP-4166: THE FOLDED TORUS picture of cluster deserts.

Fold a desert of the cluster {v_1 < ... < v_c} by the slowest comb: x = v_1*t mod 1,
period index k = floor(v_1*t).  Then:
  * comb 1 covers the fixed band x in [-rho, rho] every period;
  * comb j covers, within period k, teeth centered at x = v_1*m/v_j mod 1 for the m with
    m/v_j in period k -- equivalently a LINE structure: successive periods shift the
    tooth pattern by alpha_j = d_j/v_1 mod 1 (the normalized difference = the SLOPE),
    teeth x-width 2*rho*(v_1/v_j), count v_j/v_1 (1 or 2 for ratio < 26/11).
  * PER-PERIOD MASS IDENTITY: each comb supplies exactly 2*rho of x-measure per period.
    Demand |G| = 1-2*rho = 11/13 per period; c-1 suppliers of 2/13 each => need
    2(c-1) >= 11: c >= 6.5: the cluster ceiling, geometrically.
  * A desert of length l = N/v_1 <=> the (k,x)-strip G x {0..N-1} fully covered
    => a 2(c-1)/11-efficient covering system => (target) robust-DMNR forces slope
    rationality at denominators <= N.

VERIFY on the known desert (consecutive 7-cluster at 1/7): the fold, the mass identity,
the slack ledger per period, the slope picture (alpha_j = j/W: near-zero slopes!), and
HOW the covering pattern uses its slack.  Then measure: N (desert length in periods)
vs the predicted K = 13.
"""
from fractions import Fraction as F

RHO = F(1, 13)

def fold_check(cluster, t_lo, t_hi):
    v1 = cluster[0]
    N = int((t_hi - t_lo) * v1)
    print(f"cluster {cluster}: desert [{float(t_lo):.5f}, {float(t_hi):.5f}] = {N} periods of 1/{v1}")
    # per-period ledger: for each period, the x-mass supplied by each comb inside G
    G_lo, G_hi = RHO, 1 - RHO
    total_overlap = F(0)
    uncovered_any = 0
    for k in range(N):
        a = t_lo + F(k, v1)
        b = a + F(1, v1)
        # collect folded teeth of all combs on [a, b), as x-intervals
        ivs = []
        for j, w in enumerate(cluster):
            m0, m1 = int(w * a) - 1, int(w * b) + 1
            for m in range(m0, m1 + 1):
                ta, tb = (F(m) - RHO) / w, (F(m) + RHO) / w
                lo, hi = max(ta, a), min(tb, b)
                if lo < hi:
                    ivs.append(((lo - a) * v1, (hi - a) * v1, j))
        # mass per comb
        mass = {}
        for lo, hi, j in ivs:
            mass[j] = mass.get(j, F(0)) + (hi - lo)
        # union within G
        gs = sorted((max(lo, G_lo), min(hi, G_hi)) for lo, hi, j in ivs if hi > G_lo and lo < G_hi)
        merged = []
        for lo, hi in gs:
            if merged and lo <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], hi))
            else:
                merged.append((lo, hi))
        cov = sum(hi - lo for lo, hi in merged)
        supply_G = sum(min(hi, G_hi) - max(lo, G_lo) for lo, hi, j in ivs if hi > G_lo and lo < G_hi)
        overlap = supply_G - cov
        total_overlap += overlap
        if cov < G_hi - G_lo - F(1, 10**9):
            uncovered_any += 1
        if k < 3 or k == N - 1:
            print(f"  period {k}: comb masses {[str(mass.get(j, 0)) for j in range(len(cluster))]}"
                  f" | G-supply {supply_G} G-covered {cov} overlap {overlap}")
    print(f"  periods with uncovered G: {uncovered_any} (must be 0 inside a desert)")
    print(f"  total overlap (slack used): {total_overlap} = {float(total_overlap):.4f}; "
          f"budget = N*(2(c-1)-11)/13 = {float(F(N * (2*(len(cluster)-1) - 11), 13)):.4f}")
    # slopes
    print(f"  slopes alpha_j = d_j/v_1: {[str(F(w - cluster[0], cluster[0])) for w in cluster[1:]]}")

# the known desert: consecutive 7-cluster at 1/7 (from S87: eps in ~[-1/546, 1/91])
W = 400
cluster = list(range(W, W + 7))
t0 = F(1, 7)
lo, hi = t0 - F(1, 546), t0 + F(1, 91)
fold_check(cluster, lo, hi)
print()
# a NON-desert stretch for contrast (generic location)
print("contrast: generic location (should show uncovered periods immediately):")
fold_check(cluster, F(1, 5), F(1, 5) + F(6, W))
