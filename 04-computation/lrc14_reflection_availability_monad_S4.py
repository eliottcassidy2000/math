#!/usr/bin/env python3
r"""
lrc14_reflection_availability_monad_S4.py   (monad-explorer-2026-07-09-S4, HYP-5737/THM-669)

THE REFLECTION AVAILABILITY LEMMA -- verification + the explicit floor tables.

  (i)  Av(E, r) = (1-r^2)^{-1} int sum_gaps (g(1-r) - 2r a - 1/7)_+ dx
       >= intW_{1/7+r}(E)/(1-r^2)                      [reflection + subadditivity]
  (ii) intW_{th'} >= int_{th'} mu_s ds                  [layer cake]
  (iii) mu_th(E_n) >= 1 - 2(n-1)(n th - 1)/n for 1/n < th <= 2/n   [parametric tent]

PARTS:
  A. exact Av(E, r) engine (cell decomposition; gap lengths AND positions affine per
     cell) + verification of (i) across the zoo x r-grid (both sides exact rationals).
  B. parametric tent floor vs exact mu_th (spot checks; THM-651 recovery at n=8, th=1/7).
  C. the explicit availability floor table over (n, r): pointwise (n <= 6), tent-
     integrated (n <= 10-ish), and where the parametric moment floors are needed.
  D. end-to-end honesty: proven floor vs empirical Av vs needed |P|/7 on the S2/S3
     battery splits.
"""
import sys
from fractions import Fraction as F
from itertools import combinations
from math import gcd

exec(open('/home/bigo/math/04-computation/lrc14_cubic_moment_gate_monad_S11.py')
     .read().split("if __name__")[0])
THETA = F(1, 7)

def sorted_phase_forms(E, x0, x1):
    """Sorted phase affine forms (alpha, beta) on the cell (x0,x1), anchored at 0.
    Returns list of (alpha_i, beta_i) with values increasing on the cell."""
    E = sorted(set(int(e) for e in E))
    xm = (x0 + x1) / 2
    forms = []
    for e in E:
        c = (e * xm).__floor__()
        forms.append((e, F(-c)))
    forms.sort(key=lambda ab: ab[0] * xm + ab[1])
    return forms

def av_exact(E, r):
    """Exact Av(E, r) = (1-r^2)^{-1} int sum_gaps (g(1-r) - 2 r a - 1/7)_+ dx."""
    E = sorted(set(int(e) for e in E))
    assert 0 in E, "anchor 0 required (observer co-offset)"
    bps = breakpoints(E)
    n = len(E)
    tot = F(0)
    one = F(1)
    for ci in range(len(bps) - 1):
        x0, x1 = bps[ci], bps[ci + 1]
        if x0 == x1:
            continue
        forms = sorted_phase_forms(E, x0, x1)
        # per gap i: a = p_i, g = p_{i+1} - p_i  (last: g = 1 - p_{n-1}, a = p_{n-1})
        terms = []
        for i in range(n):
            ai = forms[i]
            if i < n - 1:
                gi = (forms[i + 1][0] - forms[i][0], forms[i + 1][1] - forms[i][1])
            else:
                gi = (-forms[i][0], one - forms[i][1])
            # arg(x) = g(x)(1-r) - 2 r a(x) - 1/7   (affine)
            A = gi[0] * (1 - r) - 2 * r * ai[0]
            B = gi[1] * (1 - r) - 2 * r * ai[1] - F(1, 7)
            terms.append((A, B))
        # integrate sum of (Ax+B)_+ with kink subdivision
        cuts = {x0, x1}
        for (A, B) in terms:
            if A != 0:
                xc = -B / A
                if x0 < xc < x1:
                    cuts.add(xc)
        cuts = sorted(cuts)
        for si in range(len(cuts) - 1):
            u0, u1 = cuts[si], cuts[si + 1]
            um = (u0 + u1) / 2
            for (A, B) in terms:
                if A * um + B > 0:
                    tot += (A * (u0 + u1) / 2 + B) * (u1 - u0)
    return tot / (1 - r * r)

def intW_exact(E, theta):
    aV, _, _, _ = excess_moments(E, [theta], max_pow=1) if False else (None,)*4
    # direct: reuse excess machinery at threshold theta via the W-pieces integral
    E = sorted(set(int(e) for e in E))
    bps = breakpoints(E)
    tot = F(0)
    for ci in range(len(bps) - 1):
        x0, x1 = bps[ci], bps[ci + 1]
        if x0 == x1:
            continue
        gaps = cell_gaps(E, x0, x1)
        cuts = {x0, x1}
        for (A, B) in gaps:
            if A != 0:
                xc = (theta - B) / A
                if x0 < xc < x1:
                    cuts.add(xc)
        cuts = sorted(cuts)
        for si in range(len(cuts) - 1):
            u0, u1 = cuts[si], cuts[si + 1]
            um = (u0 + u1) / 2
            for (A, B) in gaps:
                if A * um + B > theta:
                    tot += (A * (u0 + u1) / 2 + B - theta) * (u1 - u0)
    return tot

def mu_exact_theta(E, theta):
    E = sorted(set(int(e) for e in E))
    bps = breakpoints(E)
    tot = F(0)
    for ci in range(len(bps) - 1):
        x0, x1 = bps[ci], bps[ci + 1]
        if x0 == x1:
            continue
        gaps = cell_gaps(E, x0, x1)
        cuts = {x0, x1}
        for (A, B) in gaps:
            if A != 0:
                xc = (theta - B) / A
                if x0 < xc < x1:
                    cuts.add(xc)
        cuts = sorted(cuts)
        for si in range(len(cuts) - 1):
            u0, u1 = cuts[si], cuts[si + 1]
            um = (u0 + u1) / 2
            if any(A * um + B > theta for (A, B) in gaps):
                tot += u1 - u0
    return tot

def tent_floor(n, theta):
    """mu_theta(E_n) >= 1 - 2(n-1)(n theta - 1)/n for 1/n < theta <= 2/n (else None)."""
    if theta <= F(1, n) or theta > F(2, n):
        return None
    return 1 - 2 * (n - 1) * (n * theta - 1) / F(n)

if __name__ == "__main__":
    print("=" * 100)
    print("PART A -- THM-669(i) verification: Av(E, r) >= intW_{1/7+r}/(1-r^2), exact")
    print("=" * 100)
    ZOO = [
        ("AP_8 co-offsets {0..7}", list(range(8))),
        ("hard-core P={9..13} cluster", [0, 1, 2, 3, 4]),
        ("kps dissoc 13", [0,1,4,9,11,16,20,23,25,28,30,33,35]),
        ("covering worst {17-v}", sorted(17 - v for v in [1,2,3,4,7,8,9,10,11,12,13,14,17])),
        ("tight AP {0..12}", list(range(13))),
        ("6-cluster", [0, 2, 5, 9, 11, 14]),
    ]
    viol = 0
    for name, E in ZOO:
        for r in (F(0), F(1, 20), F(1, 10), F(1, 7), F(1, 5), F(1, 3)):
            av = av_exact(E, r)
            rhs = intW_exact(E, F(1, 7) + r) / (1 - r * r)
            ok = av >= rhs
            if not ok:
                viol += 1
            if r in (F(0), F(1, 10), F(1, 3)):
                print(f"  {name:28s} r={str(r):>5s}: Av = {float(av):.5f} >= "
                      f"intW_lift/(1-r^2) = {float(rhs):.5f}  {'OK' if ok else 'VIOLATION'}")
        sys.stdout.flush()
    print(f"  total violations over zoo x r-grid: {viol}")

    print()
    print("=" * 100)
    print("PART B -- parametric tent floor vs exact mu_theta")
    print("=" * 100)
    for n, E in [(8, list(range(8))), (8, [0,1,2,3,4,5,6,8]), (9, list(range(9))),
                 (10, [0,1,2,4,7,9,12,14,17,19])]:
        for theta in (F(1, 7), F(8, 49), F(1, 6)):
            fl = tent_floor(n, theta)
            if fl is None:
                continue
            mu = mu_exact_theta(E, theta)
            print(f"  n={n:2d} theta={str(theta):>5s}: tent floor {float(fl):.4f} <= "
                  f"exact mu = {float(mu):.4f}  {'OK' if mu >= fl else 'VIOLATION'}")
        sys.stdout.flush()

    print()
    print("=" * 100)
    print("PART C -- the explicit availability floor table (proved, structure-free)")
    print("=" * 100)
    print("  n <= 6 (pointwise): Av >= (1/n - 1/7 - r)_+/(1-r^2)")
    for n in (4, 5, 6):
        rmax = F(1, n) - F(1, 7)
        print(f"    n={n}: positive for r < {rmax} = {float(rmax):.4f}")
    print("  7 <= n <= 10 (tent-integrated): Av >= (1-r^2)^-1 int_{1/7+r}^{2/n} tentfloor(n,s) ds")
    for n in (8, 9, 10):
        for r in (F(0), F(1, 50), F(1, 25)):
            th0 = F(1, 7) + r
            if th0 >= F(2, n):
                print(f"    n={n}, r={str(r):>5s}: window empty (1/7+r >= 2/n)")
                continue
            # integral of 1 - 2(n-1)(ns-1)/n ds from th0 to s1 where floor hits 0
            # floor = 0 at s* = (1 + n/(2(n-1)))/n
            sstar = (1 + F(n, 2 * (n - 1))) / n
            s1 = min(F(2, n), sstar)
            if s1 <= th0:
                print(f"    n={n}, r={str(r):>5s}: tent window empty")
                continue
            # antiderivative: s - 2(n-1)(n s^2/2 - s)/n = s - (n-1)(n s^2 - 2 s)/n... compute exactly
            def anti(s):
                return s - 2 * (n - 1) * (n * s * s / 2 - s) / F(n)
            val = (anti(s1) - anti(th0)) / (1 - r * r)
            print(f"    n={n}, r={str(r):>5s}: Av >= {val} = {float(val):.6f}")
        sys.stdout.flush()

    print()
    print("=" * 100)
    print("PART D -- proven floor vs empirical Av on battery splits (honesty table)")
    print("=" * 100)
    for name, E_L, r in [
        ("S2 split L={0,14,28,42}/182", [0, 14, 28, 42], F(45, 182)),
        ("hard-core L 5-teeth", [0, 1, 2, 3, 4], F(1, 10)),
        ("8-teeth cluster", list(range(8)), F(1, 25)),
    ]:
        # dilation-normalize E_L (Av is dilation-invariant: phases identical)
        g = 0
        for e in E_L[1:]:
            g = gcd(g, e)
        En = [e // g if g > 1 else e for e in E_L]
        av = av_exact(En, r)
        n = len(En)
        pw = max(F(0), (F(1, n) - F(1, 7) - r)) / (1 - r * r) if n <= 6 else F(0)
        th0 = F(1, 7) + r
        tent_av = F(0)
        if F(1, n) < th0 <= F(2, n):
            sstar = (1 + F(n, 2 * (n - 1))) / n
            s1 = min(F(2, n), sstar)
            if s1 > th0:
                def anti(s):
                    return s - 2 * (n - 1) * (n * s * s / 2 - s) / F(n)
                tent_av = (anti(s1) - anti(th0)) / (1 - r * r)
        floor = max(pw, tent_av)
        print(f"  {name:30s} n={n} r={float(r):.3f}: PROVEN floor = {float(floor):.5f}, "
              f"empirical Av = {float(av):.5f}, ratio = {float(av/floor) if floor > 0 else float('inf'):.1f}x")
        sys.stdout.flush()
