#!/usr/bin/env python3
r"""
lrc14_sharp_aliasing_tv_monad_S1.py   (monad-explorer-2026-07-09-S1, HYP-5707 / THM-665)

THE SHARP ALIASING BOUND + THE EXACT TV(W') LEDGER -- the fleet's named next brick
(kps-S108/S112 NEXT lines) for the sole remaining Part-A node.

THEOREM (THM-665, proof in the canon file):
  W(x) = sum_i (g_i(x) - 1/7)_+  is 1-periodic, continuous, piecewise linear.  For the
  Vmax-ruler grid average E_grid[W](V) = (1/V) sum_{j=0}^{V-1} W(j/V):
    (i)   E_grid[W](V) = sum_{m in Z} What(mV)                    [Poisson aliasing, EXACT]
    (ii)  |What(r)| <= TV(W') / (4 pi^2 r^2)   for r != 0         [BV coefficient bound]
    (iii) |R_grid| = |E_grid[W](V) - int_0^1 W| <= TV(W') * zeta(2) * 2 / (4 pi^2 V^2)
                   = TV(W') / (12 V^2).
  This is 4pi^2 ~ 39.5x sharper than the circulated pi^2/3 envelope (kps-S108), so the
  a-priori good-period-existence window  V > V_0 = sqrt(TV(W')/(12 * int W))  shrinks by
  2pi ~ 6.28x relative to V_0^env = sqrt(TV(W') pi^2/(3 int W)).

THIS SCRIPT:
  1. EXACT TV(W') per cluster (order-cell engine + theta-kink subdivision; W' is piecewise
     CONSTANT INTEGER, TV = sum |jumps| including the periodic seam).  Settles the true
     scale of TV(W') (kps's 2.8*spread window implicitly assumed TV ~ 0.3 spread^2; the
     naive collision-count bound is ~ 2 sum (e_i-e_j)^2 ~ 100x bigger).
  2. VALIDATION: for kps-S108's two clusters, measure R_grid at V = 200..6400 and check
     |R_grid| <= TV/(12 V^2) (and compare to the pi^2/3 envelope).
  3. ALIASING IDENTITY CHECK: What(mV) by exact per-cell Fourier integrals of the linear
     pieces; partial sums vs measured R_grid.
  4. THE LEDGER: per-shape TV(W'), int W, V_0 (sharp vs envelope), V_0/spread over the
     covering zoo (klein-S206 worst-margin covering set, mac-mini's 7-structured @91,
     kps-S108 clusters, tight AP reference).
"""
import sys
from fractions import Fraction as F
from itertools import combinations
from math import gcd, pi, sin, cos, sqrt
import numpy as np

exec(open('/home/bigo/math/04-computation/lrc14_cubic_moment_gate_monad_S11.py')
     .read().split("if __name__")[0])
THETA = F(1, 7)

# ---------------------------------------------------------------------------
# exact subcell decomposition: [0,1] -> intervals on which W is affine
# returns list of (x0, x1, slope:int(F), intercept:F) for W restricted there
# ---------------------------------------------------------------------------

def W_pieces(E, theta=THETA):
    E = sorted(set(int(e) for e in E))
    bps = breakpoints(E)
    pieces = []
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
            P, Q = F(0), F(0)
            for (A, B) in gaps:
                if A * um + B > theta:
                    P += A
                    Q += B - theta
            pieces.append((u0, u1, P, Q))
    # merge zero-length / consecutive same-slope pieces are fine to keep separate
    return pieces


def tv_wprime(pieces):
    """Total variation of W' over the circle: sum of |slope jumps| at internal
    boundaries + the periodic seam."""
    slopes = [p[2] for p in pieces]
    tv = F(0)
    for i in range(len(slopes)):
        tv += abs(slopes[i] - slopes[i - 1])   # i-1 = -1 wraps to the last piece = seam
    return tv


def intW(pieces):
    tot = F(0)
    for (u0, u1, P, Q) in pieces:
        tot += (P * (u0 + u1) / 2 + Q) * (u1 - u0)
    return tot


def What(pieces, r):
    """Exact-formula Fourier coefficient of the piecewise-linear W at integer r != 0
    (float evaluation of the closed form)."""
    s = 0.0 + 0.0j
    w = 2 * pi * r
    for (u0, u1, P, Q) in pieces:
        a, b = float(P), float(Q)
        x0, x1 = float(u0), float(u1)
        # int (a x + b) e^{-i w x} dx over [x0,x1]
        e0 = complex(cos(w * x0), -sin(w * x0))
        e1 = complex(cos(w * x1), -sin(w * x1))
        term = ((a * x1 + b) * e1 - (a * x0 + b) * e0) / (-1j * w) - a * (e1 - e0) / (-1j * w) ** 2
        s += term
    return s


def grid_avg_W(E, V):
    Ea = np.array(E, float)
    j = np.arange(0, V)
    ph = np.mod(np.outer(j, Ea), V) / V
    ph.sort(axis=1)
    g = np.concatenate([np.diff(ph, axis=1), (ph[:, 0] + 1 - ph[:, -1])[:, None]], axis=1)
    return float(np.maximum(g - 1/7, 0).sum(axis=1).mean())


if __name__ == "__main__":
    ZOO = [
        ("kps-S108 dissoc (spr 35)",  [0,1,4,9,11,16,20,23,25,28,30,33,35]),
        ("kps-S108 7-struct (spr 82)",[0,7,14,21,26,29,37,44,51,58,67,75,82]),
        ("klein-S206 worst covering", sorted(17 - v for v in [1,2,3,4,7,8,9,10,11,12,13,14,17])),
        ("mac-mini @91 7-struct",     sorted(91 - v for v in [91,84,77,70,65,62,54,47,40,33,24,16,9])),
        ("tight AP (non-cov ref)",    list(range(13))),
        ("GW co-offsets (24-v)",      sorted(24 - v for v in [1,2,3,4,5,6,7,8,9,10,11,13,24])),
    ]

    print("=" * 104)
    print("PART 1 -- THE EXACT TV(W') LEDGER  (+ int W, V0 sharp vs envelope)")
    print("=" * 104)
    print(f"{'cluster':30s} {'spread':>6} {'TV(Wp)':>10} {'TV/spr^2':>9} {'intW':>8} "
          f"{'V0sharp':>8} {'V0env':>8} {'V0/spr':>7}")
    ledger = {}
    for name, E in ZOO:
        pieces = W_pieces(E)
        tv = tv_wprime(pieces)
        iw = intW(pieces)
        spread = max(E) - min(E)
        v0s = sqrt(float(tv) / (12 * float(iw))) if iw > 0 else float('inf')
        v0e = sqrt(float(tv) * pi * pi / (3 * float(iw))) if iw > 0 else float('inf')
        ledger[name] = (E, pieces, tv, iw, spread, v0s)
        print(f"{name:30s} {spread:>6} {float(tv):>10.1f} {float(tv)/spread**2:>9.3f} "
              f"{float(iw):>8.5f} {v0s:>8.1f} {v0e:>8.1f} {v0s/spread:>7.2f}")
        sys.stdout.flush()

    print()
    print("=" * 104)
    print("PART 2 -- BOUND VALIDATION on kps-S108's clusters: |R_grid| vs TV/(12V^2)")
    print("=" * 104)
    for name in list(ledger)[:2]:
        E, pieces, tv, iw, spread, _ = ledger[name]
        print(f"[{name}] int W = {float(iw):.6f}, TV(W') = {float(tv):.1f}")
        print(f"  {'V':>6} {'R_grid':>12} {'TV/(12V^2)':>12} {'ok?':>4} {'env pi^2/3':>12} {'sharp gain':>10}")
        for V in (200, 400, 800, 1600, 3200, 6400):
            R = grid_avg_W(E, V) - float(iw)
            bnd = float(tv) / (12 * V * V)
            env = float(tv) * pi * pi / (3 * V * V)
            print(f"  {V:>6} {R:>+12.2e} {bnd:>12.2e} {'OK' if abs(R) <= bnd else 'FAIL':>4} "
                  f"{env:>12.2e} {bnd/abs(R) if R else float('inf'):>10.1f}x")
        sys.stdout.flush()

    print()
    print("=" * 104)
    print("PART 3 -- ALIASING IDENTITY CHECK: sum_m What(mV) vs measured R_grid (V=400)")
    print("=" * 104)
    for name in list(ledger)[:2]:
        E, pieces, tv, iw, spread, _ = ledger[name]
        V = 400
        R = grid_avg_W(E, V) - float(iw)
        S = 0.0
        for m in range(1, 60):
            S += 2 * What(pieces, m * V).real
        print(f"[{name}] V={V}: measured R_grid = {R:+.3e}; sum_(|m|<=59) What(mV) = {S:+.3e}; "
              f"agree to {abs(R - S):.1e}")
        sys.stdout.flush()

    print()
    print("=" * 104)
    print("PART 4 -- THE WINDOW: V0_sharp vs the drift-embed threshold 1.41*spread per cluster")
    print("=" * 104)
    for name, (E, pieces, tv, iw, spread, v0s) in ledger.items():
        if float(iw) <= 0:
            continue
        drift = 1.41 * spread
        verdict = ("EXISTENCE window INSIDE drift threshold (embed governs)" if v0s <= drift
                   else f"existence needs V > {v0s:.0f} = {v0s/spread:.2f}*spread (above drift 1.41)")
        print(f"  {name:30s} V0_sharp = {v0s:>7.1f}  1.41*spread = {drift:>7.1f}  -> {verdict}")
    sys.stdout.flush()
