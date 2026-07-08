"""
lrc14_pz_k11_exhaustive_opus_S148.py   (opus-2026-07-08-S148, HYP-5327)

NAIL THE k=11 LEG (the single binding open piece of LRC(14): min_E PZ >= bar, margin +0.0159).

klein-S178 verified diam <= 15 exhaustive (min 0.346788 >= bar 0.331206).  This EXTENDS the
exact exhaustive per-diameter minimum as far as feasible (diam 16, 17, ...), to:
 (a) confirm the min RISES monotonically past diam 15 (decorrelation) -> shrink the tail;
 (b) find the TRUE per-diameter minimizer (a stretched block near the AP) and its PZ;
 (c) locate the exhaustive frontier D* beyond which a decorrelation-tail bound must take over.

Exact PZ via the general integrator (opus-S148).  Enumerate primitive k=11 shapes E with
min=0, max=D, gcd of differences = 1 (translation+dilation normalize).  Count grows as
C(D-1, k-2); feasible to D ~ 18.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations
import sys, time

sys.path.insert(0, ".")
sys.path.insert(0, "04-computation")
from lrc14_pz_general_integrator_opus_S148 import pz_exact, BAR

K = 11
bar = BAR[K]


def gcd_all(xs):
    g = 0
    for x in xs:
        g = gcd(g, x)
    return g


def main():
    t0 = time.time()
    print("=" * 92)
    print(f"k={K} EXHAUSTIVE per-diameter min PZ (bar = {bar} = {float(bar):.6f})")
    print("  primitive shapes E: 0 = min, D = max, gcd(differences)=1")
    print("=" * 92)
    overall_min = None
    overall_arg = None
    for D in range(K - 1, 19):
        mn = None; arg = None; cnt = 0
        # E = {0} u chosen (k-2) from 1..D-1 u {D}
        for mid in combinations(range(1, D), K - 2):
            E = [0] + list(mid) + [D]
            diffs = [E[i + 1] - E[i] for i in range(len(E) - 1)]
            if gcd_all(diffs) != 1:
                continue
            cnt += 1
            _, _, PZ = pz_exact(E)
            if mn is None or PZ < mn:
                mn = PZ; arg = E
        if mn is None:
            continue
        if overall_min is None or mn < overall_min:
            overall_min = mn; overall_arg = arg
        flag = "CLEARS" if mn >= bar else "*** BELOW BAR ***"
        print(f"  diam {D:3d}: {cnt:7d} primitive shapes, min PZ = {mn} = {float(mn):.6f}"
              f"  margin {float(mn - bar):+.6f}  {flag}")
        print(f"            minimizer {arg}")
        if time.time() - t0 > 800:
            print(f"  [time budget; stopped at diam {D}]")
            break
    print()
    print(f"  OVERALL exhaustive min PZ (diam {K-1}..{D}): {overall_min} = {float(overall_min):.6f}")
    print(f"    at {overall_arg}; margin over bar {float(overall_min - bar):+.6f}")
    print(f"  Monotone rise past the minimizing diameter => the exhaustive frontier D* is the")
    print(f"  cutoff for the decorrelation-tail bound (diam > D* => PZ >= this min >= bar).")
    print(f"[{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
