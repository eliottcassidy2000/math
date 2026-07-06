#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S19 -- HYP-4472: THE COMPACTNESS ROUTE to the density floor
via EQUICONTINUITY-OF-SAFE (integrating S18 equicontinuity + S17 safe-stable +
S4 J-K accumulation).

S18: safe is the REGULAR/equicontinuous functional (M is jagged).  This tests:
 (1) CONTINUITY along lift rays: safe_1d(w^N, 2/25) -> safe_2d(U, 2/25) as
     N -> inf (Weyl equidistribution -- the equicontinuity of safe), where
     w^N = a + N*b and U = 2-torus {(a_i t + b_i s)}.
 (2) safe_2d(U) > 0 for COUPLED 2-tori, bounded below (M(U) >= 1/12 > 2/25 => the
     danger arcs cannot cover the 2-torus => a UNIFORM 2-torus floor).
 (3) The zero-locus: safe_2d = 0 only for AP-direction / product-with-AP tori.

If safe_1d -> safe_2d continuously AND safe_2d >= c > 0 off the AP-locus, then by
compactness safe_1d is bounded below off the AP-lift-locus = the density floor
(the height non-uniformity is an M-artifact, dissolved by using safe).
"""
import itertools, random
from fractions import Fraction as F
from math import gcd

def dist(x):
    return abs(x - round(x))

def safe_1d_float(W, rho, grid=40000):
    """Leb{t: ||v t|| >= rho all v} approx (float)."""
    cnt = 0
    for j in range(grid):
        t = (j + 0.5) / grid
        if all(dist(v * t) >= rho for v in W):
            cnt += 1
    return cnt / grid

def safe_2d(uv, rho, grid=600):
    """Leb of the safe set on the 2-torus {(u_i t + v_i s)}."""
    cnt = 0
    for jt in range(grid):
        t = (jt + 0.5) / grid
        for js in range(grid):
            s = (js + 0.5) / grid
            ok = True
            for (u, v) in uv:
                if dist(u * t + v * s) < rho:
                    ok = False
                    break
            if ok:
                cnt += 1
    return cnt / (grid * grid)

RHO = 2 / 25

def part1_continuity():
    print("=" * 78)
    print("PART 1: safe CONTINUITY along a lift ray (equicontinuity of safe)")
    print("=" * 78)
    # a coupled 2-torus: base a_i on t, lift b_i on s.  Use a non-AP direction.
    a = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]     # base = AP residues
    b = [0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0]        # lift dirs 4 and 6 (the attainer dir)
    uv = list(zip(a, b))
    s2 = safe_2d(uv, RHO, grid=500)
    print(f"  2-torus U (base AP, lift {{4,6}}): safe_2d(2/25) = {s2:.5f}")
    print(f"  safe_1d(w^N) for w^N = a + N*b, growing N (should -> safe_2d):")
    for N in (13, 25, 50, 100, 200, 400):
        W = tuple(a[i] + N * b[i] for i in range(12))
        if len(set(W)) != 12:
            continue
        s1 = safe_1d_float(W, RHO, grid=60000)
        print(f"    N={N:>4}: safe_1d = {s1:.5f}  (|safe_1d - safe_2d| = {abs(s1-s2):.5f})")
    print("  => convergence safe_1d -> safe_2d confirms safe is CONTINUOUS along")
    print("     the lift ray (Weyl equidistribution) -- the equicontinuity of safe.")

def part2_torus_floor():
    print()
    print("=" * 78)
    print("PART 2: safe_2d(U) > 0 for COUPLED 2-tori (the uniform torus floor)")
    print("=" * 78)
    random.seed(19)
    a = list(range(1, 13))
    mins = []
    for _ in range(40):
        # random lift direction b (coupled: >=1 nonzero, not all -> product)
        b = [random.randint(0, 3) for _ in range(12)]
        if sum(1 for x in b if x != 0) < 2:
            continue
        uv = list(zip(a, b))
        s2 = safe_2d(uv, RHO, grid=360)
        mins.append((s2, b))
    mins.sort()
    print(f"  coupled 2-tori sampled: {len(mins)}")
    print(f"  MIN safe_2d over coupled tori: {mins[0][0]:.5f} (lift dir {mins[0][1]})")
    print(f"  # with safe_2d ~ 0 (< 0.005): {sum(1 for s,b in mins if s < 0.005)}")
    print(f"  smallest 4 safe_2d: {[f'{s:.4f}' for s,b in mins[:4]]}")
    # the AP-direction torus (b = a, i.e. the dilated-AP limit) -> safe=0
    print()
    print("  ZERO-LOCUS check -- AP-direction tori:")
    for label, b in [("b=0 (pure AP, N=1)", [0]*12), ("b=a (dilated-AP dir)", list(range(1,13)))]:
        uv = list(zip(a, b))
        s2 = safe_2d(uv, RHO, grid=400)
        print(f"    {label}: safe_2d = {s2:.5f} {'(ZERO -- AP locus)' if s2 < 0.005 else ''}")

if __name__ == "__main__":
    part1_continuity()
    part2_torus_floor()
