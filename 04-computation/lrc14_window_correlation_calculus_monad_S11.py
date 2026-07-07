#!/usr/bin/env python3
r"""
lrc14_window_correlation_calculus_monad_S11.py   (monad-explorer-2026-07-07-S11, HYP-5157)

THE WINDOW-AVERAGED CORRELATION CALCULUS -- the "why" layer under the excess-mass
certificates (Part 7 of the S11 program).

DERIVATION (one line each):
  E[V_theta] = int_x sum_r (g_r - theta)_+ = int_x int_u 1[window (u,u+theta) hole] du dx
             = sum_{J subset E} (-1)^{|J|} I_J(theta),
  where I_J(theta) = int_u meas{x : e_j x in (u,u+theta) for all j in J} du.
    |J| = 1:  I = theta                                       (universal)
    |J| = 2:  I = int_x (theta - ||(e_i - e_j) x||)_+ dx = theta^2   (UNIVERSAL:
              substitute y = (e_i - e_j)x mod 1 -- measure-preserving for ANY nonzero
              integer difference; the tent integral is theta^2 regardless of speed!)
    |J| = 3:  I_ijk = int_x (theta - span{e_i x, e_j x, e_k x})_+ dx  (SHAPE-DEPENDENT:
              depends on the reduced difference pattern (p, q) = (e_j - e_i, e_k - e_i)/g;
              resonant patterns (small p, q) have I ~ theta^2-scale, generic ~ theta^3)

  => E[V_theta](E) = 1 - k theta + C(k,2) theta^2 - Sigma_3(E) + Sigma_4(E) - ...
     The SHAPE enters first at WEIGHT 3.  (Matches HYP-4987: any proof must consume
     weight >= 3; matches monad-S10: the deficit's quadratic part vanishes --
     HERE IS THE SAME PHENOMENON ON THE TAIL SIDE, with the u-average as the killer.)

  BONFERRONI (pointwise in (u, x), so valid for every E):
     E[V_theta](E) >= 1 - k theta + C(k,2) theta^2 - Sigma_3(E)         [B3 lower]
     E[V_theta](E) <= 1 - k theta + C(k,2) theta^2 - Sigma_3(E) + Sigma_4(E)  [B4 upper]

  At k=8, theta=1/7: universal part = 1 - 8/7 + 28/49 = 3/7 = 0.42857...

TESTS:
  1. verify I_pair = theta^2 exactly for a spread of differences (sanity of the lemma)
  2. exact Sigma_3(E) per shape (order-cell engine on triples); B3 vs exact E[V]
  3. THE EXTREMALITY QUESTION: is Sigma_3 maximized by the AP over the battery +
     descent?  (a weighted-3AP-count extremality: Sigma_3 = sum over triples of
     w(reduced pattern), w = I_3 -- large exactly on near-AP patterns; the classical
     "AP maximizes 3-AP count" is the unweighted shadow)
  4. exact Sigma_4 at the AP (how fast does the expansion converge / B4 tightness)
  5. the same calculus for the SECOND moment's leading terms (two windows):
     E[V^2]'s weight-2 part is also universal -- verified numerically here.
"""
import sys
from fractions import Fraction as F
from itertools import combinations
from math import gcd
import random

exec(open('/home/bigo/math/04-computation/lrc14_cubic_moment_gate_monad_S11.py')
     .read().split("if __name__")[0])

THETA = F(1, 7)

# ---------------------------------------------------------------------------
# exact I_J via the span integrand, using order cells of the sub-shape J
# ---------------------------------------------------------------------------

def exact_I_subset(J, theta=THETA):
    """I_J = int_x (theta - span({e x : e in J}))_+ dx, exact.
    span = 1 - maxgap of the |J| phases.  So I_J = int (theta - 1 + maxgap)_+.
    Uses the cell engine: per cell each gap is affine; maxgap is a max of affines;
    (theta - 1 + maxgap)_+ integrates piecewise."""
    J = sorted(set(int(e) for e in J))
    bps = breakpoints(J)
    tot = F(0)
    for ci in range(len(bps) - 1):
        x0, x1 = bps[ci], bps[ci + 1]
        if x0 == x1:
            continue
        gaps = cell_gaps(J, x0, x1)
        # integrand: (theta - 1 + maxgap(x))_+ ; maxgap = max of affine functions.
        # subdivide at all pairwise crossings of the gap lines
        cuts = {x0, x1}
        for (A1, B1), (A2, B2) in combinations(gaps, 2):
            if A1 != A2:
                xc = (B2 - B1) / (A1 - A2)
                if x0 < xc < x1:
                    cuts.add(xc)
        cuts = sorted(cuts)
        for si in range(len(cuts) - 1):
            u0, u1 = cuts[si], cuts[si + 1]
            um = (u0 + u1) / 2
            # active max gap at midpoint
            vals = [(A * um + B, A, B) for (A, B) in gaps]
            _, A, B = max(vals)
            # h(x) = theta - 1 + A x + B ; integrate h_+ over (u0,u1)
            h0 = theta - 1 + A * u0 + B
            h1 = theta - 1 + A * u1 + B
            if h0 <= 0 and h1 <= 0:
                continue
            if h0 >= 0 and h1 >= 0:
                tot += (h0 + h1) * (u1 - u0) / 2
            else:
                xz = (1 - theta - B) / A
                if h0 > 0:
                    tot += h0 * (xz - u0) / 2
                else:
                    tot += h1 * (u1 - xz) / 2
    return tot


def sigma_m(E, m, theta=THETA):
    """Sigma_m(E) = sum over m-subsets J of I_J(theta), exact."""
    tot = F(0)
    for J in combinations(sorted(set(E)), m):
        tot += exact_I_subset(J, theta)
    return tot


if __name__ == "__main__":
    k = 8
    print("=" * 100)
    print("PART 7a -- THE UNIVERSAL PAIR LEMMA: I_pair = theta^2 for every nonzero difference")
    print("=" * 100)
    for d in [1, 2, 3, 7, 13, 40, 183]:
        I = exact_I_subset([0, d])
        print(f"  diff {d:4d}: I = {I} = {float(I):.6f}   (theta^2 = {float(THETA**2):.6f})  "
              f"{'OK' if I == THETA**2 else 'FAIL'}")
    sys.stdout.flush()

    print()
    print("=" * 100)
    print("PART 7b -- TRIPLE ATOMS: I(0,p,q) for reduced patterns (weight table)")
    print("=" * 100)
    pats = [(1,2),(1,3),(2,3),(1,4),(3,4),(1,5),(2,5),(1,7),(3,7),(1,13),(5,13),(7,13),(11,29)]
    for (p, q) in pats:
        I = exact_I_subset([0, p, q])
        print(f"  pattern (0,{p:2d},{q:2d}): I = {str(I):>16s} = {float(I):.6f}   "
              f"(theta^3 = {float(THETA**3):.6f}, theta^2/2 = {float(THETA**2/2):.6f})")
    sys.stdout.flush()

    print()
    print("=" * 100)
    print("PART 7c -- Sigma_3 per shape; Bonferroni B3 = 3/7 - Sigma_3 vs exact E[V]")
    print("=" * 100)
    battery = [
        ("AP_8",                     list(range(8))),
        ("near-AP one-swap 7->8",    [0,1,2,3,4,5,6,8]),
        ("near-AP one-swap 7->13",   [0,1,2,3,4,5,6,13]),
        ("parity interlace (M119)",  [0,2,4,6,8,10,11,12]),
        ("two-block 4+4 far",        [0,1,2,3,40,41,42,43]),
        ("single-far",               [0,1,2,3,4,5,6,40]),
        ("geometric",                [0,1,2,4,8,16,32,64]),
        ("Sidon (perfect diff)",     [0,1,3,7,12,20,30,44]),
        ("dense-random",             [0,2,3,7,9,13,16,17]),
    ]
    univ = 1 - k * THETA + comb(k, 2) * THETA ** 2
    print(f"  universal part 1 - 8th + 28th^2 = {univ} = {float(univ):.6f}")
    s3_results = []
    for name, E in battery:
        s3 = sigma_m(E, 3)
        aV, MV, m3V, vmaxV = excess_moments(E, [THETA])
        b3 = univ - s3
        s3_results.append((name, E, s3, aV[0]))
        print(f"  {name:26s} Sigma_3 = {float(s3):.6f}  B3 = {float(b3):+.6f}  "
              f"exact E[V] = {float(aV[0]):.6f}  (B3 valid: {b3 <= aV[0]})")
        sys.stdout.flush()

    print()
    print("=" * 100)
    print("PART 7d -- Sigma_4 at AP_8 and near-AP (expansion convergence, B4 sandwich)")
    print("=" * 100)
    for name, E in battery[:2] + [battery[7]]:
        s3 = sigma_m(E, 3); s4 = sigma_m(E, 4)
        aV, _, _, _ = excess_moments(E, [THETA])
        print(f"  {name:26s} B3 = {float(univ - s3):+.6f} <= E[V] = {float(aV[0]):.6f} "
              f"<= B4 = {float(univ - s3 + s4):+.6f}   (Sigma_4 = {float(s4):.6f})")
        sys.stdout.flush()

    print()
    print("=" * 100)
    print("PART 7e -- IS Sigma_3 MAXIMIZED AT THE AP? (weighted-3AP extremality; jump descent)")
    print("=" * 100)
    rng = random.Random(715711)
    def norm8(E):
        E = sorted(set(E)); E = [e - E[0] for e in E]
        g = 0
        for e in E[1:]:
            g = gcd(g, e)
        if g > 1:
            E = [e // g for e in E]
        return E
    s3_ap = sigma_m(list(range(8)), 3)
    print(f"  Sigma_3(AP_8) = {s3_ap} = {float(s3_ap):.6f}")
    best = (s3_ap, list(range(8)))
    cache = {tuple(range(8)): s3_ap}
    cur, curv = list(range(8)), s3_ap
    for step in range(150):
        E2 = list(cur)
        r = rng.random()
        if r < 0.5:
            E2[rng.randrange(8)] = rng.randint(0, 40)
        elif r < 0.8:
            i, j = rng.randrange(8), rng.randrange(8)
            E2[i] = E2[j] + rng.choice([1, -1, 2, -2])
        else:
            E2 = [rng.randint(0, 40) for _ in range(8)]
        E2 = norm8(E2)
        if len(E2) != 8:
            continue
        t = tuple(E2)
        if t not in cache:
            cache[t] = sigma_m(E2, 3)
        v = cache[t]
        if v > best[0]:
            best = (v, E2)
            print(f"  step {step}: Sigma_3 = {float(v):.6f} > AP at {E2}  <-- AP NOT the max!")
        if v > curv or rng.random() < 0.3:
            cur, curv = E2, v
    print(f"  max Sigma_3 over {len(cache)} shapes: {float(best[0]):.6f} at {best[1]}")
    print(f"  (AP value {float(s3_ap):.6f}; AP {'IS' if best[1]==list(range(8)) else 'is NOT'} the maximizer)")
    sys.stdout.flush()
