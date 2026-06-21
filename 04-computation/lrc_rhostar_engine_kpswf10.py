#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_rhostar_engine_kpswf10.py   (kind-pasteur 2026-06-21, THM-527 Thread A)

EXACT-rational engine for the lonely-density

    rho*(P, E) = meas{ x in [0,1) :
                       ||p x|| >= 1/14  for all p in P        (x in G_P)
                       AND  the k phases {frac(e_i x): e_i in E}
                            have circular max-gap > 2/7 }.

This is THE single-variable object of THM-527. Two ingredients:

  (A) G_P = { x : ||p x|| >= 1/14, all p in P }.
      For each p, ||p x|| >= 1/14  <=>  frac(p x) in [1/14, 13/14].
      Boundaries of frac(p x) at thresholds 1/14, 13/14 occur at
          x = (14 j + 1)/(14 p)  and  x = (14 j + 13)/(14 p),  j = 0..p-1.
      So all G_P breakpoints live on the grid  t/(14 p),  t in Z.

  (B) GOOD(E) = { x : max circular gap of {frac(e_i x)} > 2/7 }.
      maxgap(x) is piecewise-constant-pattern; it changes (i) when two
      phases collide: (e_i - e_j) x in Z, i.e. x = t/(e_i - e_j); and the
      INDICATOR [maxgap>2/7] additionally flips where some gap = 2/7 exactly,
      i.e. frac(e_j x) - frac(e_i x) = +-2/7  =>  (e_j-e_i) x = +-2/7 + Z
      =>  x = (7 t +- 2)/(7 (e_j - e_i)).
      All such breakpoints are rational with denominator dividing 7*lcm-diffs.

  COMBINED breakpoint set B = { all G_P breaks } U { all GOOD breaks } U {0,1}.
  On each open cell (a,b) between consecutive breakpoints, BOTH the G_P
  membership and the [maxgap>2/7] indicator are constant; evaluate at the
  midpoint (an exact Fraction) and accumulate (b-a) when both hold.

  rho* is then an EXACT Fraction.

Cross-checks built in:
  * measS7 complement:  meas{maxgap <= 2/7} should equal measS7? NO --
    measS7 = meas{all 7 sectors hit}.  "all sectors hit" => no gap of a
    full empty sector, but maxgap<=2/7 is the GAP condition, a different
    (finer) event.  We instead check GOOD against an independent brute grid.
  * pure consec mu_k = meas(GOOD({0..k-1}))  must hit the canon table
    (mu_4=19/21, mu_13=829/4620).
  * the 1/84 consec floor: min_P meas(GOOD({0..8}) cap G_P) over P subset
    {1..13}, |P| = 13-9 = 4, attained P={1,2,3,12}.
"""
import itertools
from fractions import Fraction as Fr
from math import gcd


# ============================================================ circular max-gap
def circ_maxgap_at(E, x):
    """Exact circular max-gap of the multiset {frac(e*x): e in E}, x a Fraction.
       Returns a Fraction in [0,1]. (If all phases equal, gap = 1.)"""
    phases = sorted(set((Fr(e) * x) % 1 for e in E))
    if len(phases) == 1:
        return Fr(1)
    g = Fr(0)
    for a, b in zip(phases, phases[1:]):
        if b - a > g:
            g = b - a
    # wrap-around gap
    wrap = (phases[0] + 1) - phases[-1]
    if wrap > g:
        g = wrap
    return g


def in_GP_at(P, x, thr=Fr(1, 14)):
    """||p x|| >= thr for all p in P ?  x a Fraction."""
    for p in P:
        f = (Fr(p) * x) % 1
        d = f if f <= Fr(1, 2) else 1 - f   # ||p x||
        if d < thr:
            return False
    return True


# ============================================================ breakpoint sets
def gp_breaks(P):
    """G_P breakpoints: x = (14 j +- 1)/(14 p) in (0,1)."""
    bps = set()
    for p in P:
        if p == 0:
            continue
        # frac(p x) = 1/14 or 13/14  =>  p x = m + 1/14 or m + 13/14
        # x = (14 m + 1)/(14 p), (14 m + 13)/(14 p)
        for m in range(0, p):
            for r in (1, 13):
                v = Fr(14 * m + r, 14 * p)
                if 0 < v < 1:
                    bps.add(v)
    return bps


def good_breaks(E):
    """GOOD(E) breakpoints in (0,1): phase collisions x=t/d and gap=+-2/7
       crossings x=(7t+-2)/(7 d), for every pairwise difference d=|e_i-e_j|."""
    bps = set()
    diffs = set()
    Elist = list(E)
    for i in range(len(Elist)):
        for j in range(i + 1, len(Elist)):
            d = abs(Elist[i] - Elist[j])
            if d != 0:
                diffs.add(d)
    for d in diffs:
        # collisions: x = t/d
        for t in range(1, d):
            bps.add(Fr(t, d))
        # gap = 2/7 crossings: (e_j - e_i) x = +-2/7 + integer
        #   x = (7 m +- 2)/(7 d)
        hi = 7 * d
        for m in range(0, hi + 1):
            for s in (2, -2):
                num = 7 * m + s
                v = Fr(num, 7 * d)
                if 0 < v < 1:
                    bps.add(v)
    return bps


# ============================================================ exact rho*
def rho_star(P, E, thr=Fr(1, 14), gapthr=Fr(2, 7), return_arcs=False):
    """Exact meas{ x in G_P AND maxgap{frac(e x)} > gapthr }."""
    P = [int(p) for p in P]
    E = [int(e) for e in E]
    bps = {Fr(0), Fr(1)}
    bps |= gp_breaks(P)
    bps |= good_breaks(E)
    pts = sorted(bps)
    total = Fr(0)
    arcs = []
    for a, b in zip(pts, pts[1:]):
        mid = (a + b) / 2
        if in_GP_at(P, mid, thr) and circ_maxgap_at(E, mid) > gapthr:
            total += (b - a)
            if return_arcs:
                arcs.append((a, b))
    if return_arcs:
        return total, arcs
    return total


def mu_pure(E, gapthr=Fr(2, 7)):
    """Pure cluster measure meas{maxgap{frac(e x)} > gapthr}, no G_P."""
    E = [int(e) for e in E]
    bps = {Fr(0), Fr(1)} | good_breaks(E)
    pts = sorted(bps)
    total = Fr(0)
    for a, b in zip(pts, pts[1:]):
        mid = (a + b) / 2
        if circ_maxgap_at(E, mid) > gapthr:
            total += (b - a)
    return total


def meas_GP(P, thr=Fr(1, 14)):
    """Exact meas(G_P)."""
    P = [int(p) for p in P]
    bps = {Fr(0), Fr(1)} | gp_breaks(P)
    pts = sorted(bps)
    total = Fr(0)
    for a, b in zip(pts, pts[1:]):
        mid = (a + b) / 2
        if in_GP_at(P, mid, thr):
            total += (b - a)
    return total


# ============================================================ validation suite
def brute_mu(E, N=20160, gapthr=2.0 / 7.0):
    """Numeric maxgap measure on a uniform grid, for cross-check."""
    cnt = 0
    for t in range(N):
        x = (t + 0.5) / N
        ph = sorted(((e * x) % 1.0) for e in E)
        ph = sorted(set(round(v, 12) for v in ph))
        if len(ph) == 1:
            g = 1.0
        else:
            g = max(b - a for a, b in zip(ph, ph[1:]))
            g = max(g, ph[0] + 1 - ph[-1])
        if g > gapthr + 1e-12:
            cnt += 1
    return cnt / N


def main():
    print("=" * 78)
    print("THM-527 Thread A: EXACT rho*(P,E) engine  (kind-pasteur kpswf10)")
    print("=" * 78)

    # ---- (1) pure consec mu_k vs canon table ----
    print("\n--- (1) pure consecutive mu_k = meas(maxgap{j x: j<k} > 2/7) ---")
    canon = {3: Fr(1), 4: Fr(19, 21), 5: Fr(9, 14), 13: Fr(829, 4620)}
    for k in [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]:
        E = list(range(k))
        m = mu_pure(E)
        tag = ""
        if k in canon:
            tag = f"  canon={canon[k]}  {'OK' if m == canon[k] else 'MISMATCH!!'}"
        print(f"  k={k:2d}: mu = {m} = {float(m):.6f}{tag}")

    # ---- (2) brute cross-check on a few shapes ----
    print("\n--- (2) exact vs brute (numeric grid) cross-check ---")
    shapes = [
        [0, 1, 2, 3, 4, 5, 6, 7, 8],          # consec k=9
        [0, 2, 3, 4, 5, 6, 8],                # k=7 perforated {0..8}\{1,7}
        [0, 1, 3, 5, 6, 8, 10, 11, 12, 14],   # k=10 perforated near-AP
        [0, 5, 7, 8, 9, 10, 11, 13, 18],      # k=9 spread-18
    ]
    for E in shapes:
        ex = mu_pure(E)
        br = brute_mu(E)
        ok = abs(float(ex) - br) < 2e-3
        print(f"  E={E}\n      exact={float(ex):.6f}  brute={br:.6f}  {'OK' if ok else 'CHECK'}")

    # ---- (3) meas(G_P) caps cross-check ----
    print("\n--- (3) meas(G_P) for canonical worst P (cap cross-check) ---")
    # cap_k = min over |P|=13-k of meas(G_P); canon worst examples
    for P in [[1, 2, 3, 12], [1, 2, 3], [1, 2], [1]]:
        m = meas_GP(P)
        print(f"  P={P}: meas(G_P) = {m} = {float(m):.6f}")

    # ---- (4) the 1/84 consecutive floor ----
    print("\n--- (4) consec k=9 floor: min_P meas(GOOD({0..8}) cap G_P), |P|=4 ---")
    E = list(range(9))
    best = None
    best_P = None
    for P in itertools.combinations(range(1, 14), 4):
        r = rho_star(list(P), E)
        if best is None or r < best:
            best = r
            best_P = P
    print(f"  min rho* = {best} = {float(best):.6f}  at P={best_P}")
    print(f"  target 1/84 = {Fr(1,84)} = {float(Fr(1,84)):.6f}  "
          f"{'OK' if best == Fr(1,84) else 'DIFFERS'}")

    # ---- (5) k=3 unconditional: rho* = meas(G_P) (every x good) ----
    print("\n--- (5) k=3 check: GOOD({0,1,2}) = all of [0,1)? rho*=meas(G_P)? ---")
    E3 = [0, 1, 2]
    print(f"  mu_pure({E3}) = {mu_pure(E3)}  (expect 1)")
    for P in [[1, 2, 3, 12], [4, 6, 9]]:
        r = rho_star(P, E3)
        g = meas_GP(P)
        print(f"  P={P}: rho*={float(r):.6f}  meas(G_P)={float(g):.6f}  "
              f"{'EQUAL (all-good)' if r == g else 'NOT equal'}")

    print("\nDONE.")


if __name__ == "__main__":
    main()
