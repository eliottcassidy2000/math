"""
lrc14_cap_is_pair_avoidance_kps.py  (kind-pasteur-2026-06-27-S31ag)

HYP-3090 refinement: cap_k = C(k+1,2)/C(14,2). With j = |P| = 13-k (small-part size),
this is min_{|P|=j} meas(lonely(P)) = C(14-j,2)/C(14,2) = (14-j)(13-j)/(14*13).

That is a PAIR-AVOIDANCE probability (2 of 14 apex points both miss the j danger combs).
Test directly: compute meas(lonely(P)) = meas{x: ||p x|| >= 1/14 for all p in P} exactly,
minimize over j-element integer sets P, and check min == C(14-j,2)/C(14,2). Identify the
minimizer (is it the AP {1..j}?). Confirms the cap is the order-2 (pairwise) prediction
exactly for j<=3 and reveals where (j>=4) it breaks.
"""
import sys, itertools
from fractions import Fraction as F
from math import comb

def lonely_measure(P, thr=F(1, 14)):
    """meas{x in [0,1): ||p x|| >= thr for all p in P}, exact."""
    P = [p for p in P if p != 0]
    if not P:
        return F(1)
    bps = set([F(0), F(1)])
    for p in P:
        for k in range(p + 1):
            lo = F(k) / p - thr / p
            hi = F(k) / p + thr / p
            for v in (lo, hi):
                if 0 <= v <= 1:
                    bps.add(v)
    bps = sorted(bps)
    total = F(0)
    for a, b in zip(bps, bps[1:]):
        if b <= a:
            continue
        mid = (a + b) / 2
        if all(min((p*mid) % 1, 1 - (p*mid) % 1) >= thr for p in P):
            total += b - a
    return total

if __name__ == "__main__":
    sys.stdout.reconfigure(line_buffering=True)
    C142 = comb(14, 2)
    print("min_{|P|=j} meas(lonely(P)) vs C(14-j,2)/C(14,2):")
    print(f"{'j':>2} {'predicted C(14-j,2)/91':>22} {'computed min':>16} {'minimizer P':>22} match")
    for j in range(1, 6):
        pred = F(comb(14 - j, 2), C142)
        best = None; bestP = None
        # search P subset of {1..16}, size j (the apex denominator is 14, search a bit beyond)
        for P in itertools.combinations(range(1, 17), j):
            m = lonely_measure(P)
            if best is None or m < best:
                best = m; bestP = P
        match = "EXACT" if best == pred else f"dev={float(best-pred):+.5f}"
        print(f"{j:>2} {str(pred):>22} {str(best):>16} {str(bestP):>22} {match}")
    print()
    print("Interpretation: cap is the PAIRWISE (order-2) avoidance prob, exact for j<=3;")
    print("j>=4 (k<=9) is where 3+-point correlation dips it below = the LRC(14) difficulty.")
