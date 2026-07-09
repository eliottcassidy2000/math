#!/usr/bin/env python3
r"""
lrc14_parametric_D3_rerun_monad_S5.py   (monad-explorer-2026-07-09-S5, HYP-5747/THM-670)

THE PARAMETRIC-THETA D3/E[W] RERUN for n = 11..13 -- the one remaining input of
THM-669's assembly -- plus THM-670 (the Lipschitz transfer).

  THM-670 (proved in canon file): for theta >= 1/7, at most floor(1/theta) <= 7 gaps
  exceed theta, so dW_theta/dtheta >= -7 pointwise ==>
      E[W_{theta2}] >= E[W_{theta1}] - 7 (theta2 - theta1)      (1/7 <= th1 <= th2).

  D3 (opus-S148 / THM-661, threshold-agnostic; W = W_theta', M = max W):
      mu_theta'(E) >= D3_theta'(E) = E[W]/M + (E[W] - E[W^2]/M)^2 / (E[W^2] - E[W^3]/M).

PARTS:
  A. calibration: D3 at theta = 1/7 vs THM-661's exact compact minima (k=12: 0.355876
     at (0,2,4,5,6..12,14); k=13: 0.308844 at (0,2,3,4..12,14)).
  B. THE LEDGER: compact scans (prim-diam <= 15, float-first grid, exact confirm of
     the minima) of E[W_theta'] and D3_theta' at theta' = 1/7 + r,
     r in {0, 1/50, 1/25, 1/14}, n = 11, 12, 13.
  C. Lipschitz cross-check: ledger min E[W_theta'] vs [min E[W_1/7] - 7r] (THM-670
     never violated; how conservative).
  D. dilated-AP ladders: D3_theta' vs diameter (tail rises at lifted theta too).
  E. assembly margins: |L| in {11,12,13} splits, |P| <= 2 exact G_P^c, the
     THM-669 + ledger criterion margin per battery family.
"""
import sys
from fractions import Fraction as F
from itertools import combinations
from math import gcd
import numpy as np

exec(open('/home/bigo/math/04-computation/lrc14_cubic_moment_gate_monad_S11.py')
     .read().split("if __name__")[0])

def wmoments_numeric(E, theta, ngrid=40_000):
    Ea = np.array(sorted(set(E)), float)
    xs = (np.arange(ngrid) + 0.5) / ngrid
    ph = np.mod(np.outer(xs, Ea), 1.0)
    ph.sort(axis=1)
    g = np.concatenate([np.diff(ph, axis=1), (ph[:, 0] + 1 - ph[:, -1])[:, None]], axis=1)
    W = np.clip(g - theta, 0, None).sum(axis=1)
    return W.mean(), (W**2).mean(), (W**3).mean(), W.max()

def d3_value(m1, m2, m3, M):
    if m1 <= 0 or M <= 0:
        return 0.0
    base = m1 / M
    num = m1 - m2 / M
    den = m2 - m3 / M
    if den <= 0:
        return base
    return base + num * num / den

def wmoments_exact(E, theta):
    """Exact E[W_theta^i], i=1..3, and exact max W_theta, via the cell engine."""
    E = sorted(set(int(e) for e in E))
    bps = breakpoints(E)
    m1 = m2 = m3 = F(0)
    wmax = F(0)
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
            # integrate (Px+Q)^i on (u0,u1)
            w0, w1 = P * u0 + Q, P * u1 + Q
            wmax = max(wmax, w0, w1)
            d = u1 - u0
            m1 += (w0 + w1) * d / 2
            m2 += (w0*w0 + w0*w1 + w1*w1) * d / 3
            m3 += (w0**3 + w0*w0*w1 + w0*w1*w1 + w1**3) * d / 4
    return m1, m2, m3, wmax

def compact_shapes(n, dmax):
    """primitive 0-anchored n-clusters with diameter <= dmax (diam >= n-1)."""
    for D in range(n - 1, dmax + 1):
        for mid in combinations(range(1, D), n - 2):
            E = (0,) + mid + (D,)
            g = 0
            for e in E[1:]:
                g = gcd(g, e)
            if g == 1:
                yield E

if __name__ == "__main__":
    TH0 = 1/7
    print("=" * 100)
    print("PART A -- CALIBRATION at theta = 1/7 vs THM-661")
    print("=" * 100)
    for n, Emin, ref in [(12, (0,2,4,5,6,7,8,9,10,11,12,14), 0.355876),
                         (13, (0,2,3,4,5,6,7,8,9,10,11,12,14), 0.308844)]:
        m1, m2, m3, M = wmoments_numeric(Emin, TH0, 120_000)
        d3 = d3_value(m1, m2, m3, M)
        print(f"  n={n}: D3(THM-661 minimizer) numeric = {d3:.6f} (ref {ref}) "
              f"{'OK' if abs(d3-ref) < 2e-3 else 'CHECK'}")
    sys.stdout.flush()

    print()
    print("=" * 100)
    print("PART B -- THE PARAMETRIC LEDGER (compact scan prim-diam <= 15, float-first)")
    print("=" * 100)
    RS = [(0, F(0)), ('1/50', F(1, 50)), ('1/25', F(1, 25)), ('1/14', F(1, 14))]
    ledger = {}
    for n in (11, 12, 13):
        shapes = list(compact_shapes(n, 15))
        print(f"  n={n}: {len(shapes)} primitive compact shapes (diam <= 15)")
        for rname, r in RS:
            th = TH0 + float(r)
            best_ew = (1e9, None)
            best_d3 = (1e9, None)
            for E in shapes:
                m1, m2, m3, M = wmoments_numeric(E, th, 20_000)
                d3 = d3_value(m1, m2, m3, M)
                if m1 < best_ew[0]:
                    best_ew = (m1, E)
                if d3 < best_d3[0]:
                    best_d3 = (d3, E)
            # exact confirm the minima
            e1, e2, e3, eM = wmoments_exact(list(best_ew[1]), F(1,7) + r)
            ew_exact = e1
            f1, f2, f3, fM = wmoments_exact(list(best_d3[1]), F(1,7) + r)
            d3_exact = d3_value(float(f1), float(f2), float(f3), float(fM))
            ledger[(n, str(r))] = (float(ew_exact), d3_exact)
            print(f"    r={rname:>5}: min E[W] = {float(ew_exact):.6f} at {best_ew[1]}")
            print(f"             min D3   = {d3_exact:.6f} at {best_d3[1]}")
            sys.stdout.flush()

    print()
    print("=" * 100)
    print("PART C -- THM-670 LIPSCHITZ CROSS-CHECK: min E[W_th'] vs min E[W_1/7] - 7r")
    print("=" * 100)
    for n in (11, 12, 13):
        base = ledger[(n, '0')][0]
        for rname, r in RS[1:]:
            got = ledger[(n, str(r))][0]
            trans = base - 7 * float(r)
            print(f"  n={n}, r={rname:>5}: ledger {got:.5f} >= transfer {trans:.5f} "
                  f"{'OK' if got >= trans - 1e-9 else 'VIOLATION'} "
                  f"(headroom {got - trans:+.5f})")
    sys.stdout.flush()

    print()
    print("=" * 100)
    print("PART D -- dilated-AP ladders: D3_theta' rises with diameter (tail behavior)")
    print("=" * 100)
    for n in (12, 13):
        for rname, r in [('0', F(0)), ('1/25', F(1, 25))]:
            th = TH0 + float(r)
            row = []
            for c in (1, 2, 3, 5, 8):
                E = [c * i for i in range(n)]
                m1, m2, m3, M = wmoments_numeric(E, th, 30_000)
                row.append(d3_value(m1, m2, m3, M))
            # dilation-invariant? (should be exactly invariant) -- use stretched non-dilated ladder
            row2 = []
            for D in (n - 1, 2 * n, 4 * n, 8 * n):
                E = sorted(set([0] + [round(i * D / (n - 1)) for i in range(1, n)]))
                if len(E) != n:
                    row2.append(float('nan')); continue
                m1, m2, m3, M = wmoments_numeric(E, th, 30_000)
                row2.append(d3_value(m1, m2, m3, M))
            print(f"  n={n}, r={rname:>4}: dilates {['%.4f' % v for v in row]} (invariance check); "
                  f"stretched-AP diam {n-1}/{2*n}/{4*n}/{8*n}: {['%.4f' % v for v in row2]}")
            sys.stdout.flush()

    print()
    print("=" * 100)
    print("PART E -- ASSEMBLY MARGINS for |L| in {11,12,13} splits (|P| <= 2)")
    print("  criterion: Av_floor(n, r) = ledger-min E[W_{1/7+r}]/(1-r^2)  vs  meas(G_P^c) exact")
    print("=" * 100)
    for name, v in [("14*(AP\\6)+83 full cluster", [14,28,42,56,70,83,98,112,126,140,154,168,182]),
                    ("covering worst 17-set", [1,2,3,4,7,8,9,10,11,12,13,14,17]),
                    ("@91 7-struct", [9,16,24,33,40,47,54,62,65,70,77,84,91])]:
        v = sorted(v)
        Vmax = v[-1]
        for nL in (13, 12, 11):
            # split: L = nL largest speeds (smallest co-offsets)
            L = v[-nL:]
            P = v[:-nL]
            S_L = Vmax - L[0]
            r = S_L / Vmax
            if r >= 0.999:
                continue
            rkey = min(RS[1:], key=lambda t: abs(float(t[1]) - r)) if r > 0.005 else RS[0]
            # use the ledger at the NEAREST r >= actual (conservative if ledger-r >= r)
            usable = [(rn, rv) for rn, rv in RS if float(rv) >= r - 1e-12]
            if not usable:
                print(f"  [{name}] |L|={nL}: r = {r:.3f} beyond ledger range -- needs direct")
                continue
            rn, rv = usable[0]
            ewfloor = ledger[(nL, str(rv))][0]
            av_floor = ewfloor / (1 - float(rv)**2)
            gpc = sum(min(1.0, vp / 7) for vp in P) if P else 0.0   # crude cap
            # exact meas(G_P^c) for small P
            if P:
                bads = []
                for vp in P:
                    for m in range(vp + 1):
                        bads.append((max(0.0, (m - 1/14) / vp), min(1.0, (m + 1/14) / vp)))
                bads.sort()
                tot = 0.0
                cur = 0.0
                for lo, hi in bads:
                    if hi <= cur:
                        continue
                    tot += hi - max(lo, cur)
                    cur = max(cur, hi)
                gpc = tot
            print(f"  [{name}] |L|={nL} (r = {r:.3f} -> ledger r = {rn}): "
                  f"Av_floor = {av_floor:.5f} vs meas(G_P^c) = {gpc:.5f} "
                  f"{'-> CLOSES (mass criterion, modulo grid port)' if av_floor > gpc else '-> short (needs sharper split/exact Av)'}")
        sys.stdout.flush()
