#!/usr/bin/env python3
"""procgen_mlr_20260930_crowd.py -- Part B (crowded times) of the multiplicative lonely runner lane.

B1  Triangle Lemma (two generators, threshold R <= D/6): exhaustive verification on Z^2 windows,
    plus the sharpness control R ~ D/4 and the one-generator run lemma.
B2  Uniqueness Lemma (one small numerator per Diophantine window).
B3  Coset crowding bound: near fraction of every coset of <2,3> in (Z/D)^* vs the PROVED bound.
B4  Crowd census in boxes: theta-crowded times vs their best numerator; the single-triangle law;
    test of the uniform-U dichotomy (box version) on primes and Collatz gates.
Prints to stdout; raises on any failed check.
"""
import math
import random
import sys
import time
from fractions import Fraction as F

import numpy as np

import procgen_mlr_20260930_core as core
from procgen_mlr_20260930_core import check

T0 = time.time()


def hdr(s):
    print("\n" + "=" * 100 + "\n" + s + "\n" + "=" * 100)
    sys.stdout.flush()


def b1_triangle():
    hdr("B1. Triangle Lemma (R <= D/4): near set {|r(j,k)| < R} in Z^2 = disjoint union of triangles T(P0) with exact residues")
    rng = random.Random(20260930)
    total = 0
    ninst = 0
    for D in (101, 343, 1001, 1009, 3125, 4001, 7553):
        check(D % 2 and D % 3, "D must be coprime to 6")
        for den in (4, 6, 16, 40):
            R = D // den
            if R < 2:
                continue
            hs = range(1, D) if D <= 343 else [rng.randrange(1, D) for _ in range(40)]
            W = int(math.log2(R)) + 4
            for h in hs:
                total += core.verify_triangle_lemma(D, h, R, -W, 2 * W, -W, 2 * W)
                ninst += 1
        sys.stdout.flush()
    print(f"  verified (b) component = full triangle of its unique minimum with exact residues r = 2^x 3^y u,")
    print(f"  (c) apex numerator u coprime to 6, (d) shadows exact, pairwise disjoint, disjoint from the near set:")
    print(f"  {total} interior components in {ninst} (D, h, R) instances, R/D in {{1/4, 1/6, 1/16, 1/40}}, "
          f"D in {{101, 7^3, 7*11*13, 1009, 5^5, 4001, 7553}}")
    # sharpness control: the engine of the lemma is "adjacent near points are related by an exact x2 / x3";
    # for R <= D/4 this always holds, above D/4 the x3 step can wrap onto a near point
    for frac in (F(1, 4), F(3, 10), F(1, 3)):
        nviol = 0
        nadj = 0
        for D in (1009, 4001):
            R = int(D * frac)
            for h in range(1, 200):
                res = {(j, k): core.centred(h * core.mult(j, k, D), D) for j in range(-3, 9) for k in range(-3, 7)}
                for (j, k), r in res.items():
                    if abs(r) >= R:
                        continue
                    for (Q, g) in (((j + 1, k), 2), ((j, k + 1), 3)):
                        if Q in res and abs(res[Q]) < R:
                            nadj += 1
                            if res[Q] != g * r:
                                nviol += 1
        print(f"  engine check R = {frac} D: {nviol} of {nadj} adjacent near pairs violate the exact x2/x3 relation")
        if frac <= F(1, 4):
            check(nviol == 0, "exact relation violated at R <= D/4")
        else:
            check(nviol > 0, "expected violations above D/4")
    # one-generator run lemma: if |r(j)| < R <= D/4 for j = 0..t then r(j) = 2^j r(0) and |r(0)| < R 2^-t
    nrun = 0
    for D in (1009, 4001, 30011):
        R = D // 4
        for h in range(1, 3000):
            r = [core.centred(h * pow(2, j, D), D) for j in range(40)]
            j = 0
            while j < 40:
                if abs(r[j]) < R:
                    t = j
                    while t + 1 < 40 and abs(r[t + 1]) < R:
                        t += 1
                    for i in range(j, t + 1):
                        check(r[i] == 2 ** (i - j) * r[j], "run lemma: inexact doubling")
                    check(abs(r[j]) * 2 ** (t - j) < R, "run lemma bound")
                    nrun += 1
                    j = t + 1
                else:
                    j += 1
    print(f"  one-generator run lemma (R <= D/4): every near run [j0,j1] has r(j) = 2^(j-j0) r(j0) and "
          f"|r(j0)| < R 2^-(j1-j0): {nrun} runs checked")


def b2_uniqueness():
    hdr("B2. Uniqueness Lemma: two numerators |u0|,|u1| <= U at P0, P1 with U 2^|dj| 3^|dk| < D/2 are the same orbit point")
    rng = random.Random(7)
    npairs = 0
    nsame_apex = 0
    for D in (1009, 4001, 30011, 100003):
        hs = range(1, D) if D <= 4001 else [rng.randrange(1, D) for _ in range(3000)]
        for h in hs:
            for U in (3, 12):
                pts = []
                for j in range(-7, 8):
                    for k in range(-5, 6):
                        r = core.centred(h * core.mult(j, k, D), D)
                        if abs(r) <= U:
                            pts.append((j, k, r))
                for i in range(len(pts)):
                    for l in range(i + 1, len(pts)):
                        j0, k0, u0 = pts[i]
                        j1, k1, u1 = pts[l]
                        if U * 2 ** abs(j1 - j0) * 3 ** abs(k1 - k0) < D / 2:
                            check(F(u0) * F(2) ** (j1 - j0) * F(3) ** (k1 - k0) == u1, "uniqueness lemma violated")
                            npairs += 1
                            if u0 % 2 and u0 % 3 and u1 % 2 and u1 % 3:
                                nsame_apex += 1
    check(nsame_apex == 0, "two distinct 6-free numerators inside one window")
    print(f"  {npairs} pairs inside the window checked (D in {{1009, 4001, 30011, 100003}}, U in {{3, 12}}): all satisfy")
    print(f"  u1 = 2^(j1-j0) 3^(k1-k0) u0 exactly; no window contains two 6-free numerators (two apexes): {nsame_apex}")


def b3_coset_bound():
    hdr("B3. Coset crowding bound: near fraction of a coset hH (H = <2,3>) vs max_{u >= m, 6-free} sigma(R/u)/sigma(D/2u)")
    nco = 0
    worst = 0.0
    most = []
    for D in range(5, 1500):
        if D % 2 == 0 or D % 3 == 0:
            continue
        H, cos = core.cosets23(D)
        for den in (4, 8, 16, 32):
            R = D // den
            if R < 2:
                continue
            for h, c in cos:
                cc = np.minimum(c, D - c)
                rho = F(int((cc < R).sum()), len(H))
                m = int(cc.min())
                bnd = core.coset_crowding_bound(D, R, m)
                check(rho <= bnd, f"coset bound violated D={D} h={h} R={R}: {rho} > {bnd}")
                nco += 1
                if bnd > 0:
                    worst = max(worst, float(rho / bnd))
                if den == 16 and len(H) >= 20:
                    most.append((float(rho), D, h, m, len(H), float(bnd)))
    most.sort(reverse=True)
    print(f"  {nco} (coset, R) pairs with D < 1500 coprime to 6, R = floor(D/4), D/8, D/16, D/32: bound holds "
          f"everywhere (max rho/bound = {worst:.3f})")
    print("  most crowded cosets at R = D/16 with |H| >= 20 (rho, D, h, least |element| m, |H|, bound):")
    for row in most[:8]:
        print(f"    rho={row[0]:.3f}  D={row[1]:5d}  h={row[2]:4d}  m={row[3]:3d}  |H|={row[4]:4d}  bound={row[5]:.3f}")
    # explicit consequence: rho >= theta forces m <= ...
    for theta in (0.25, 0.5):
        print(f"  consequence (R = D/16): rho(hH) >= {theta} forces the least element m of hH to satisfy "
              f"sigma(R/u)/sigma(D/2u) >= {theta} for some 6-free u >= m, e.g. for D = 10^6 + 3: m <= "
              f"{max_m_for_theta(1000003, 1000003 // 16, theta)}")


def max_m_for_theta(D, R, theta):
    """largest m such that the PROVED bound (max over 6-free u >= m) is still >= theta."""
    best = 0
    vals = []
    for u in range(1, R):
        if u % 2 and u % 3:
            vals.append((u, core.sigma3(F(R, u)) / max(1, core.sigma3(F(D, 2 * u)))))
    for u, v in vals:
        if v >= theta:
            best = u
    return best


def b4_census():
    hdr("B4. Crowd census in boxes: theta-crowded times h vs best numerator u*(h) = min_box |r|")
    print("  near fraction rho_R(h) = #{(j,k) in box : |r| < R}/n;  background mean = 2R/D")
    targets = [("prime", 10007), ("prime", 100003), ("prime", 1000003), ("gate (20,12)", 517135), ("gate (21,13)", abs(2 ** 21 - 3 ** 13))]
    results = {}
    for name, D in targets:
        for c in (0.5, 1.0, 2.0):
            J = max(1, round(c * math.log(D) / core.LOG2))
            K = max(1, round(c * math.log(D) / core.LOG3))
            n = J * K
            Rs = (D // 8, D // 16, D // 32)
            best, cnt = core.scan_box(D, J, K, Rs)
            best = best[1:]
            for R in Rs:
                fr = cnt[R][1:] / n
                bg = 2 * R / D
                line = []
                for dth in (0.10, 0.20, 0.30):
                    th = bg + dth
                    sel = fr >= th
                    if sel.any():
                        um = int(best[sel].max())
                        line.append(f"theta=bg+{dth:.2f}: #h={int(sel.sum()):6d} max u*={um:6d} (eta={math.log(um) / math.log(D) if um > 1 else 0:.2f})")
                        results.setdefault((name, c, R * 1.0 / D, dth), um)
                    else:
                        line.append(f"theta=bg+{dth:.2f}: none")
                print(f"  {name:13s} D={D:8d} c={c:3.1f} box({J:2d},{K:2d}) R=D/{round(D / R):2d}: max rho={fr.max():.3f}; " + "; ".join(line))
            # single-triangle law for h = u (apex at the corner)
            R = D // 16
            law = []
            for u in (1, 5, 25, 125, 625):
                pred = sum(1 for j in range(J) for k in range(K) if (2 ** j) * (3 ** k) * u < R)
                law.append(f"u={u}: rho={cnt[R][u] / n:.3f} tri={pred / n:.3f}")
            print(f"      single-triangle law at R=D/16 (rho(h=u) vs #{{2^j3^k u < R}}/n): " + ", ".join(law))
            sys.stdout.flush()
    # growth of max u* among crowded times with D (box c = 1): refutes a uniform U(theta, delta)
    print("\n  uniform-U test (box c = 1, R = D/16, theta = bg + 0.10 / 0.20):")
    for dth in (0.10, 0.20):
        seq = [(D, results.get((name, 1.0, (D // 16) * 1.0 / D, dth))) for name, D in targets[:3]]
        print(f"    theta=bg+{dth:.2f}: max u* over crowded h: " + ", ".join(f"D={D}: {u}" for D, u in seq))
        vals = [u for _, u in seq if u is not None]
        if len(vals) == 3:
            check(vals[2] > vals[0], "expected max u* among crowded h to grow with D")
    print("    => in boxes matched to log D the max numerator of a theta-crowded time grows with D (power-like):")
    print("       a uniform bound U(theta, delta) independent of D FAILS in the box setting (EMPIRICAL REFUTATION)")


def main():
    core.lean_malloc()
    b1_triangle()
    b2_uniqueness()
    b3_coset_bound()
    b4_census()
    print(f"\n[Part B done in {time.time() - T0:.0f}s, peak RSS {core.mem_mb():.0f} MB]")
    check(core.mem_mb() < 500, "memory budget exceeded")
    print("PART B CHECKS PASSED")


if __name__ == "__main__":
    main()
