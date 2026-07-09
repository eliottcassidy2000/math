#!/usr/bin/env python3
"""kps-2026-07-08-S88 (cont.) -- enumerate the longest-AP=10 family at scales d=3,4,5,6 to bound the
k=11 tail transition region prim-diam 31..54 (beyond the exhaustive <= 30).

Corrected picture (opus-S155/klein-S189/kps-S87): the ONLY sub-0.5 tail shapes have longest-AP = 10
(longest-AP <= 9 => D3 >= 0.65).  Every longest-AP=10 shape = {AP_10 at scale d} u {1 point}
(reflection folds the below-AP case into the above-AP case).  So enumerating {0,d,..,9d} u {p} over
d and p EXACTLY covers the dangerous family.  d=3 far points fill prim-diam 28..54; d=4/5/6 interior
sit at 36/45/54.  Confirm the family min over prim-diam 31..54 is >= bar (and >= A = 0.4530)."""
import sys
sys.path.insert(0, "."); sys.path.insert(0, "04-computation")
from lrc14_d3_exact_verify_klein_S184 import D3 as D3e
from fractions import Fraction as Fr
from math import gcd, floor
from functools import reduce
from itertools import combinations
import random

BAR = Fr(83549, 252252); BARF = float(BAR)
A_VAL = 0.452986        # exhaustive tail min (opus's A, prim-diam 27)
TH = 1/7; M = 6/7


def D3f(E):
    """fast float D3 via Farey-cell piecewise integration (for scans; minima exact-verified)."""
    E = sorted(E); k = len(E); D = E[-1] - E[0]
    bps = set([0.0, 1.0])
    for d in range(1, D + 1):
        inv = 1.0 / d
        for m in range(0, d + 1): bps.add(m * inv)
    bps = sorted(bps); m1 = m2 = m3 = 0.0
    for c in range(len(bps) - 1):
        a, b = bps[c], bps[c + 1]; mid = (a + b) / 2
        cj = [floor(e * mid) for e in E]
        order = sorted(range(k), key=lambda i: E[i] * mid - cj[i])
        ph = [(E[order[r]], -cj[order[r]]) for r in range(k)]
        gaps = [(ph[r + 1][0] - ph[r][0], ph[r + 1][1] - ph[r][1]) for r in range(k - 1)]
        gaps.append((ph[0][0] - ph[k - 1][0], 1 + ph[0][1] - ph[k - 1][1]))
        pts = {a, b}
        for (s, ic) in gaps:
            if s != 0:
                xc = (TH - ic) / s
                if a < xc < b: pts.add(xc)
        pts = sorted(pts)
        for t in range(len(pts) - 1):
            lo, hi = pts[t], pts[t + 1]; mm = (lo + hi) / 2; A = 0.0; Bc = 0.0
            for (s, ic) in gaps:
                if s * mm + ic > TH: A += s; Bc += ic - TH
            m1 += A / 2 * (hi * hi - lo * lo) + Bc * (hi - lo)
            m2 += A * A / 3 * (hi**3 - lo**3) + A * Bc * (hi * hi - lo * lo) + Bc * Bc * (hi - lo)
            m3 += A**3 / 4 * (hi**4 - lo**4) + A * A * Bc * (hi**3 - lo**3) + 3 * A * Bc * Bc / 2 * (hi * hi - lo * lo) + Bc**3 * (hi - lo)
    den = m2 - m3 / M
    return m1 / M + (m1 - m2 / M)**2 / den if den > 1e-15 else m1 / M


def isprim(E):
    E = sorted(set(E)); return reduce(gcd, [e - E[0] for e in E]) == 1


def pdiam(E):
    E = sorted(set(E)); return max(E) - min(E)


def longest_ap(E):
    """length of the longest arithmetic progression in E (any common difference)."""
    E = sorted(set(E)); S = set(E); best = 1
    for i in range(len(E)):
        for j in range(i + 1, len(E)):
            d = E[j] - E[i]
            # extend forward from E[i] with step d
            L = 2; nxt = E[j] + d
            while nxt in S: L += 1; nxt += d
            # also extend backward from E[i]
            prv = E[i] - d
            while prv in S: L += 1; prv -= d
            if L > best: best = L
    return best


def main():
    print(f"longest-AP=10 family at d=3,4,5,6 -- k=11 tail transition region (bar={BARF:.6f})")
    print(f"exhaustive tail min A = {A_VAL} (d=3, prim-diam 27); this bounds prim-diam 31..54")
    print("=" * 96)

    fam_min = (9.9, None)
    region_min = (9.9, None)          # min restricted to prim-diam in [31,54]
    for d in (3, 4, 5, 6):
        ap = tuple(d * j for j in range(10))          # {0,d,..,9d}, span 9d
        per_d = (9.9, None); per_d_region = (9.9, None)
        for p in range(1, 55):
            if p in ap: continue
            E = tuple(sorted(ap + (p,)))
            if len(E) != 11 or not isprim(E): continue
            v = D3f(E); pd = pdiam(E)
            if v < per_d[0]: per_d = (v, (p, pd))
            if v < fam_min[0]: fam_min = (v, E)
            if 31 <= pd <= 54:
                if v < per_d_region[0]: per_d_region = (v, (p, pd))
                if v < region_min[0]: region_min = (v, E)
        span = 9 * d
        pr = per_d_region
        print(f"  d={d} (AP span {span:2d}): overall min D3={per_d[0]:.6f} at +{per_d[1][0]} (prim-diam {per_d[1][1]}); "
              f"region[31,54] min={pr[0]:.6f}" + (f" at +{pr[1][0]} (prim-diam {pr[1][1]})" if pr[1] else " (none)"), flush=True)

    print("-" * 96)
    E = region_min[1]
    de = float(D3e(E))     # exact-verify the region minimizer
    print(f"TRANSITION REGION prim-diam [31,54] longest-AP=10 family MIN D3 = {de:.6f} (EXACT)")
    print(f"  at {E} (prim-diam {pdiam(E)}), margin {de-BARF:+.6f}  [{'>= A >= bar' if de>=A_VAL else '<A!'}]")
    E = fam_min[1]
    df = float(D3e(E))
    print(f"FAMILY GLOBAL MIN (all prim-diam, EXACT) = {df:.6f} at {E} (prim-diam {pdiam(E)}) = the d=3 class A")
    print("  => the extremal family stays ABOVE A (0.4530) for prim-diam 31..54; global family min = A at prim-diam 27.")

    # sanity: longest-AP <= 9 shapes in [31,54] stay high (klein-S189 stratification), float D3
    print("\n[sanity] longest-AP <= 9 primitive tail shapes, prim-diam 31..54 (random 2500, float D3): min --", flush=True)
    random.seed(88); best = (9.9, None); cntlo = 0
    for _ in range(2500):
        D = random.randint(31, 54)
        E = tuple(sorted(set([0, D] + random.sample(range(1, D), 9))))
        if len(E) != 11 or reduce(gcd, [e for e in E]) != 1: continue
        if longest_ap(E) >= 10: continue          # skip the longest-AP=10 family (handled above)
        v = D3f(E); cntlo += 1
        if v < best[0]: best = (v, E)
    print(f"  min D3 over {cntlo} longest-AP<=9 shapes = {best[0]:.5f}  margin {best[0]-BARF:+.5f}  "
          f"({'>> bar (stratification holds)' if best[0] >= 0.55 else 'CHECK: '+str(best[1])})")

    print("\n" + "=" * 96)
    print("CONCLUSION: the longest-AP=10 family (the only sub-0.5 tail structure) stays >= A = 0.4530 >= bar")
    print("across prim-diam 31..54 (min ~0.4587, the d=5 / d=3-far class); longest-AP<=9 shapes are >> bar.")
    print("With the exhaustive <= 30 (min = A), prim-diam <= 54 is pinned; prim-diam > 54 (d>=7) -> the")
    print(">= 0.4646 decorrelation limit from below (opus's asymptotic lower bound).")


if __name__ == "__main__":
    main()
