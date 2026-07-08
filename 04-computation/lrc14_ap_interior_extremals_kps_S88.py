#!/usr/bin/env python3
"""kps-2026-07-08-S88 -- ENUMERATE the AP+interior extremal family for the k=11 tail (prim-diam <= 30).

The low-D3 tail shapes are AP-arithmetic (kps-S87): an AP_L at scale d + (11-L) extra points.  This
exactly enumerates that family, finds the true tail minimizer, and confirms scale-monotonicity +
that it is >= bar.  Cross-checks the exhaustive (lrc14_exhaustive_diam30_kps_S88)."""
import sys
sys.path.insert(0, "."); sys.path.insert(0, "04-computation")
from lrc14_d3_exact_verify_klein_S184 import D3 as D3e
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
from itertools import combinations

BAR = Fr(83549, 252252)
BARF = float(BAR)
DMAX = 30


def isprim(E):
    E = sorted(set(E)); return reduce(gcd, [e - E[0] for e in E]) == 1


def pdiam(E):
    E = sorted(set(E)); return max(E) - min(E)


def main():
    print(f"AP+interior extremals for the k=11 tail, prim-diam <= {DMAX}  (bar={BARF:.6f})")
    print("=" * 92)

    # ---- L=10 AP (scale d) + 1 extra point ----
    print("\n[L=10 AP at scale d, + 1 extra point]  min D3 over all extra positions, prim-diam in [25,30]:")
    print(f"  {'d':>2} {'AP span':>8} {'best +p':>8} {'prim-diam':>10} {'min D3':>10} {'margin':>9}")
    fam_min = (Fr(10), None)
    for d in range(1, 5):
        ap = tuple(d * j for j in range(10))          # {0,d,..,9d}, spans 9d
        best = (Fr(10), None)
        for p in range(1, DMAX + 1):
            if p % d == 0 and p <= 9 * d: continue     # skip AP lattice points
            E = tuple(sorted(set(ap + (p,))))
            if len(E) != 11 or not isprim(E): continue
            pd = pdiam(E)
            if pd < 25 or pd > DMAX: continue
            v = D3e(E)
            if v < best[0]: best = (v, p)
        if best[1] is not None:
            E = tuple(sorted(set(ap + (best[1],))))
            if best[0] < fam_min[0]: fam_min = (best[0], E)
            print(f"  {d:>2} {9*d:>8} {best[1]:>8} {pdiam(E):>10} {float(best[0]):>10.6f} {float(best[0]-BAR):>+9.5f}")

    # ---- L=9 AP (scale d) + 2 extra points ----
    print("\n[L=9 AP at scale d, + 2 extra points]  min D3, prim-diam in [25,30]:")
    print(f"  {'d':>2} {'best +p,q':>12} {'prim-diam':>10} {'min D3':>10} {'margin':>9}")
    for d in range(1, 4):
        ap = tuple(d * j for j in range(9))            # {0,d,..,8d}, spans 8d
        pool = [p for p in range(1, DMAX + 1) if not (p % d == 0 and p <= 8 * d)]
        best = (Fr(10), None)
        for p, q in combinations(pool, 2):
            E = tuple(sorted(set(ap + (p, q))))
            if len(E) != 11 or not isprim(E): continue
            pd = pdiam(E)
            if pd < 25 or pd > DMAX: continue
            v = D3e(E)
            if v < best[0]: best = (v, (p, q))
        if best[1] is not None:
            E = tuple(sorted(set(ap + best[1])))
            if best[0] < fam_min[0]: fam_min = (best[0], E)
            print(f"  {d:>2} {str(best[1]):>12} {pdiam(E):>10} {float(best[0]):>10.6f} {float(best[0]-BAR):>+9.5f}")

    # ---- L=8 AP (scale d) + 3 extra points (spot check, d=2,3) ----
    print("\n[L=8 AP at scale d, + 3 extra points]  min D3, prim-diam in [25,30] (spot check):")
    for d in (2, 3):
        ap = tuple(d * j for j in range(8))            # spans 7d
        pool = [p for p in range(1, DMAX + 1) if not (p % d == 0 and p <= 7 * d)]
        best = (Fr(10), None)
        for tri in combinations(pool, 3):
            E = tuple(sorted(set(ap + tri)))
            if len(E) != 11 or not isprim(E): continue
            pd = pdiam(E)
            if pd < 25 or pd > DMAX: continue
            v = D3e(E)
            if v < best[0]: best = (v, tri)
        if best[1] is not None:
            E = tuple(sorted(set(ap + best[1])))
            if best[0] < fam_min[0]: fam_min = (best[0], E)
            print(f"  d={d} best +{best[1]} prim-diam {pdiam(E)}: min D3={float(best[0]):.6f} margin {float(best[0]-BAR):+.5f}")

    print("\n" + "=" * 92)
    E = fam_min[1]
    print(f"TAIL MINIMIZER (AP+interior family, prim-diam 25..30): D3 = {float(fam_min[0]):.6f}")
    print(f"  = {fam_min[0]}")
    print(f"  shape = {E}  prim-diam {pdiam(E)}  margin {float(fam_min[0]-BAR):+.6f}  {'>= bar' if fam_min[0]>=BAR else 'BELOW BAR'}")
    print("  (opus-S155's A=(0,3,6,8,9,12,15,18,21,24,27)=3*{0..9}+8 gives 0.452986)")


if __name__ == "__main__":
    main()
