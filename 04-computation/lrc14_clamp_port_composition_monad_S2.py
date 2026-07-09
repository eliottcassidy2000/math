#!/usr/bin/env python3
r"""
lrc14_clamp_port_composition_monad_S2.py   (monad-explorer-2026-07-09-S2, HYP-5717/THM-666)

THE CLAMPED GRID-PORT + THE LOW-ANCHOR DRIFT COMPOSITION -- working hrefl.

PIECES:
  THM-666 (the clamp port; proof = THM-665's with C in place of W):
    C(x) = clamp((maxgap_E(x) - th2)/(th1 - th2), 0, 1),  th2 < th1.
    C is continuous piecewise-linear; 1[maxgap > th1] <= C <= 1[maxgap >= th2].
    Poisson aliasing + BV:  gridfrac{j : maxgap(j/V) >= th2}
                              >= E_grid[C] >= mu(th1) - TV(C')/(12 V^2).
    This ports the FULL closed-leg measure mu (0.3..0.94) to the ruler grid at 1/V^2.

  THE LOW-ANCHOR DRIFT CONDITION (from klein-S205's sharpened drift bound, drift ~ phi):
    at ruler index j, an L-tooth-free gap (a, a+g) certifies minReach >= 1/14 at
    tau = (j + a + g/2)/Vmax  iff   g(1 - r) > 1/7 + 2 r a,   r := S_L/Vmax
    (rearranged from  g > 1/7 + 2 (S_L/Vmax)(a + g/2); LOW gaps pay almost no drift).

  THE COMPOSITION CRITERION (per velocity set v, split parameter S_L):
    exists j in [0, Vmax):  [L-teeth at j admit a gap with g(1-r) > 1/7 + 2ra]
                        AND [tau_j-window inside G_P (all slow runners >= 1/14 safe)]
    ==> Mreach(v) >= 1/14, via klein-S205 Mreach_ge_of_driftGap + interval checks.
    All hypotheses exactly checkable (rationals).

THIS SCRIPT:
  A. THM-666 verification: exact TV(C') (extended cell engine: maxgap = max of affines;
     subdivide at argmax crossings + band crossings), grid-fraction vs mu(th1) - TV/(12V^2)
     on the zoo (multiple V).
  B. The composition checker over the 966 covering 13-subsets of [1,18] (which splits fire).
  C. THE r-MAP: synthetic single-scale covering families at r = spread/Vmax in a grid;
     where does the low-anchor composition stop firing (r*)?  The residual class, named.
"""
import sys, random
from fractions import Fraction as F
from itertools import combinations
from math import gcd, pi, sqrt
import numpy as np

exec(open('/home/bigo/math/04-computation/lrc14_cubic_moment_gate_monad_S11.py')
     .read().split("if __name__")[0])
THETA = F(1, 7)

# ---------------------------------------------------------------------------
# maxgap as exact piecewise-linear function; clamp C pieces; TV(C')
# ---------------------------------------------------------------------------

def maxgap_pieces(E):
    """[(x0, x1, A, B)]: maxgap(x) = A x + B on (x0, x1). Subdivides cells at
    argmax crossings of the gap lines."""
    E = sorted(set(int(e) for e in E))
    bps = breakpoints(E)
    out = []
    for ci in range(len(bps) - 1):
        x0, x1 = bps[ci], bps[ci + 1]
        if x0 == x1:
            continue
        gaps = cell_gaps(E, x0, x1)
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
            A, B = max(gaps, key=lambda ab: ab[0] * um + ab[1])
            out.append((u0, u1, A, B))
    return out


def clamp_pieces(mg_pieces, th2, th1):
    """C = clamp((maxgap - th2)/(th1 - th2), 0, 1) as exact PL pieces; returns
    (pieces [(x0,x1,slope,val0)], TV(C'), int C).  Subdivides at band crossings."""
    den = th1 - th2
    pieces = []
    for (x0, x1, A, B) in mg_pieces:
        cuts = {x0, x1}
        if A != 0:
            for t in (th1, th2):
                xc = (t - B) / A
                if x0 < xc < x1:
                    cuts.add(xc)
        cuts = sorted(cuts)
        for si in range(len(cuts) - 1):
            u0, u1 = cuts[si], cuts[si + 1]
            um = (u0 + u1) / 2
            val = A * um + B
            if val >= th1:
                pieces.append((u0, u1, F(0), F(1)))
            elif val <= th2:
                pieces.append((u0, u1, F(0), F(0)))
            else:
                pieces.append((u0, u1, A / den, (A * u0 + B - th2) / den))
    tv = F(0)
    n = len(pieces)
    for i in range(n):
        tv += abs(pieces[i][2] - pieces[i - 1][2])
    iC = F(0)
    for (u0, u1, s, v0) in pieces:
        iC += (v0 + s * (u1 - u0) / 2) * (u1 - u0)
    return pieces, tv, iC


def mu_exact(E, th):
    """meas{x : maxgap > th} exact."""
    tot = F(0)
    for (x0, x1, A, B) in maxgap_pieces(E):
        lo, hi = x0, x1
        if A == 0:
            if B <= th:
                continue
        elif A > 0:
            lo = max(lo, (th - B) / A)
        else:
            hi = min(hi, (th - B) / A)
        if hi > lo:
            tot += hi - lo
    return tot


# ---------------------------------------------------------------------------
# the composition checker
# ---------------------------------------------------------------------------

def good_p_intervals(vp, Vmax):
    """G_p = {tau in [0,1) : nearInt(vp*tau) >= 1/14} as a sorted interval list
    (exact rationals): safe between (m + 1/14)/vp and (m + 13/14)/vp."""
    out = []
    for m in range(vp):
        out.append(((F(m) + F(1, 14)) / vp, (F(m) + F(13, 14)) / vp))
    return out


def in_G_P(tau, Pspeeds):
    for vp in Pspeeds:
        y = (vp * tau) % 1
        if min(y, 1 - y) < F(1, 14):
            return False
    return True


def composition_fires(v, S_L, verbose=False):
    """v: sorted velocity list (13 ints). Split by co-offset e = Vmax - v_i <= S_L.
    Returns (fires: bool, j, gap, detail) -- checks all ruler indices j."""
    v = sorted(v)
    Vmax = v[-1]
    E_L = [Vmax - vi for vi in v if Vmax - vi <= S_L]
    P = [vi for vi in v if Vmax - vi > S_L]
    if len(E_L) < 1:
        return False, None, None, "empty L"
    r = F(S_L, Vmax)
    if r >= 1:
        return False, None, None, "r >= 1"
    for j in range(Vmax):
        # L-teeth at index j
        teeth = sorted(set(F((e * j) % Vmax, Vmax) for e in E_L))
        # circular gaps of the teeth; for each gap (a_rel, g) with LOW anchor:
        # positions measured as absolute in [0,1): gap from teeth[t] to teeth[t+1]
        nt = len(teeth)
        for t in range(nt):
            lo = teeth[t]
            hi = teeth[(t + 1) % nt] if t + 1 < nt else teeth[0] + 1
            g = hi - lo
            a = lo  # gap start (absolute phase); embed uses phi = a + g/2
            # drift condition (klein-S205 sharpened): g(1-r) > 1/7 + 2 r a  with a, g absolute
            # NOTE: the Lean takes the gap in [0,1] with a >= 0; wrap gap handled by a = teeth[last]
            if a + g > 1:
                # wrap gap: the free interval is (teeth[-1], 1) u (0, teeth[0]);
                # embed needs an interval inside [0,1]: use (teeth[-1], 1) part with the
                # observer tooth at 0 => actually e=0 is always a tooth at 0; skip wraps
                continue
            if g * (1 - r) > F(1, 7) + 2 * r * a:
                phi = a + g / 2
                tau = (j + phi) / Vmax
                if all(True for _ in [0]) and in_G_P(tau, P):
                    if verbose:
                        print(f"    fires at j={j}, gap=({float(a):.4f},{float(a+g):.4f}), "
                              f"tau={float(tau):.6f}")
                    return True, j, (a, g), f"|L|={len(E_L)}, |P|={len(P)}, r={float(r):.3f}"
    return False, None, None, f"|L|={len(E_L)}, |P|={len(P)}, r={float(r):.3f}"


def is_covering(v, n=14):
    return all(any(vi % q == 0 for vi in v) for q in range(2, n + 1))


if __name__ == "__main__":
    print("=" * 100)
    print("PART A -- THM-666 CLAMP PORT verification")
    print("=" * 100)
    ZOO = [
        ("kps dissoc (spr 35)", [0,1,4,9,11,16,20,23,25,28,30,33,35]),
        ("covering worst (spr 16)", sorted(17 - v for v in [1,2,3,4,7,8,9,10,11,12,13,14,17])),
        ("tight AP", list(range(13))),
    ]
    th1, th2 = F(1, 7) + F(1, 49), F(1, 7)   # th1 = 8/49, th2 = 1/7
    for name, E in ZOO:
        mg = maxgap_pieces(E)
        cp, tvC, iC = clamp_pieces(mg, th2, th1)
        mu1 = mu_exact(E, th1)
        print(f"[{name}] mu(8/49) = {float(mu1):.5f}, int C = {float(iC):.5f}, "
              f"TV(C') = {float(tvC):.1f}")
        for V in (400, 1600, 6400):
            # grid fraction with maxgap >= th2
            Ea = np.array(sorted(set(E)), float)
            jj = np.arange(V)
            ph = np.mod(np.outer(jj, Ea), V) / V
            ph.sort(axis=1)
            g = np.concatenate([np.diff(ph, axis=1), (ph[:, 0] + 1 - ph[:, -1])[:, None]], axis=1)
            frac = float((g.max(axis=1) >= float(th2) - 1e-15).mean())
            bound = float(mu1) - float(tvC) / (12 * V * V)
            print(f"    V={V:>5}: gridfrac(maxgap>=1/7) = {frac:.5f} >= "
                  f"mu(th1) - TV/(12V^2) = {bound:.5f}  {'OK' if frac >= bound - 1e-12 else 'FAIL'}")
        sys.stdout.flush()

    print()
    print("=" * 100)
    print("PART B -- COMPOSITION CHECKER over covering 13-subsets of [1,18]")
    print("=" * 100)
    allsets = [c for c in combinations(range(1, 19), 13) if is_covering(c)]
    print(f"covering sets: {len(allsets)}")
    fired = 0
    fails = []
    for v in allsets:
        Vmax = v[-1]
        ok = False
        for S_L in sorted(set([Vmax // 3, Vmax // 2, 2 * Vmax // 3, Vmax - 1])):
            f, j, gapinfo, det = composition_fires(list(v), S_L)
            if f:
                ok = True
                break
        if ok:
            fired += 1
        else:
            fails.append(v)
    print(f"composition fires on {fired}/{len(allsets)}; residual {len(fails)}")
    for v in fails[:8]:
        print(f"   residual: {v}")
    sys.stdout.flush()

    print()
    print("=" * 100)
    print("PART C -- THE r-MAP: synthetic single-scale covering families, r = spread/Vmax grid")
    print("=" * 100)
    rng = random.Random(97170709)
    print(f"{'r-target':>9} {'n':>4} {'fired':>6} {'frac':>6}   (low-anchor drift composition, best split)")
    for rt in (0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85):
        tried = fire = 0
        for trial in range(40):
            Vmax = rng.randint(120, 400)
            spread = int(rt * Vmax)
            # build a covering 13-set inside [Vmax - spread, Vmax]
            v = set([Vmax])
            need = list(range(2, 15))
            rng.shuffle(need)
            okbuild = True
            for q in need:
                lo = Vmax - spread
                mults = [m * q for m in range(max(1, lo // q + (1 if lo % q else 0)), Vmax // q + 1)]
                mults = [m for m in mults if m not in v]
                if not mults:
                    okbuild = False
                    break
                v.add(rng.choice(mults))
            if not okbuild:
                continue
            while len(v) < 13:
                c = rng.randint(Vmax - spread, Vmax)
                v.add(c)
            v = sorted(v)[:13]
            if len(v) != 13 or not is_covering(v):
                continue
            tried += 1
            Vm = max(v)
            got = False
            for S_L in sorted(set([Vm // 4, Vm // 3, Vm // 2, 2 * Vm // 3, int(0.8 * Vm)])):
                f, _, _, _ = composition_fires(list(v), S_L)
                if f:
                    got = True
                    break
            if got:
                fire += 1
        print(f"{rt:>9.2f} {tried:>4} {fire:>6} {fire/max(tried,1):>6.2f}")
        sys.stdout.flush()
