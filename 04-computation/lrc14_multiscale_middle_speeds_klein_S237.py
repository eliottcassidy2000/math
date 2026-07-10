#!/usr/bin/env python3
"""
THE MULTI-SCALE EXTENSION FOR MIDDLE SPEEDS (klein-2026-07-10-S237,
HYP-5910, THM-688). Completes THM-687 to arbitrary scale structure.

(I) SEPARATED SCALES: S = P u C_1 u ... u C_r, C_i = {V_i - e : e in E_i},
    ratios rho_i = V_i/V_{i+1} -> inf. Each cluster gets its own fast
    variable beta_i = frac(V_i alpha); jointly independent in the limit:

      mu(S) -> mu_inf = Int_{G_P} PROD_i m_{E_i}(alpha) d alpha,
      rate C * (1/V_r + Sum_i V_{i+1}/V_i)     (iterated two-scale slicing).

    Exact evaluator: common breakpoint refinement; per cell each m_{E_i} is
    linear (fit from interior thirds), product = exact polynomial, integrated
    exactly (Fraction coefficients, Sum c_j/(j+1)).

(II) SAME SCALE / RATIONAL RAY (the true "middle speed" case): W ~ (a/b)V
    couples through the b-fold cover. Demo b = 2 (V_1 = 2 V_2 EXACTLY):
    frac((2V_2 - e')alpha) = frac(gamma - e' alpha) with gamma = 2 beta_2 on
    the circumference-2 circle; the E_1 arcs appear TWICE (width 1/7 each,
    at e'alpha and e'alpha + 1) while the E_2 arcs double in width (2/7 at
    2e alpha). Joint fiber m_joint = meas(avoid all)/2. The PRODUCT limit is
    WRONG here; the b-cover limit is the truth -- both computed exactly and
    tested against exact mu(S).

(III) floors: all k_i <= 6 => mu_inf >= PROD(1 - k_i/7) * m_P (arcs cannot
    cover, per cluster, independently); at most one cluster can have
    k_i >= 7, leaving a single finite criterion per class.
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd

BAND_LO, BAND_HI = F(1, 14), F(13, 14)


def lcm(a, b):
    return a * b // gcd(a, b)


def mu_exact(S):
    L = 1
    for v in S:
        L = lcm(L, v)
    D2 = 28 * L
    pts = {0, D2}
    for v in S:
        s = 2 * L // v
        for k in range(v):
            pts.add((14 * k + 1) * s)
            pts.add((14 * k + 13) * s)
    pts = sorted(pts)
    good = 0
    for x, y in zip(pts, pts[1:]):
        m = x + y
        ok = True
        for v in S:
            t = v * m % (2 * D2)
            if not (2 * D2 <= 14 * t <= 26 * D2):
                ok = False
                break
        if ok:
            good += y - x
    return F(good, D2)


def in_band_frac(v, a):
    r = (v * a) % 1
    return BAND_LO <= r <= BAND_HI


def m_arcs(centers_widths, circ=1):
    """Union-complement measure of arcs (center, halfwidth) on circle of
    circumference circ; returns circ - |union| (NOT normalized)."""
    arcs = []
    for c, h in centers_widths:
        c = c % circ
        lo, hi = c - h, c + h
        if lo < 0:
            arcs.append((lo + circ, F(circ)))
            arcs.append((F(0), hi))
        elif hi > circ:
            arcs.append((lo, F(circ)))
            arcs.append((F(0), hi - circ))
        else:
            arcs.append((lo, hi))
    arcs.sort()
    tot = F(0)
    cl, ch = None, None
    for lo, hi in arcs:
        if ch is None or lo > ch:
            if ch is not None:
                tot += ch - cl
            cl, ch = lo, hi
        else:
            ch = max(ch, hi)
    if ch is not None:
        tot += ch - cl
    return circ - tot


def m_E(E, a):
    return m_arcs([((e * a) % 1, F(1, 14)) for e in E])


def cluster_breaks(E):
    pts = set()
    for e, f in combinations(E, 2):
        d = abs(e - f)
        if d:
            for k in range(d):
                pts.add(F(k, d))
                pts.add(F(7 * k + 1, 7 * d))
                pts.add(F(7 * k + 6, 7 * d))
    return pts


def P_breaks(P):
    pts = set()
    for p in P:
        for k in range(p):
            pts.add(F(14 * k + 1, 14 * p))
            pts.add(F(14 * k + 13, 14 * p))
    return pts


def poly_mul(A, B):
    out = [F(0)] * (len(A) + len(B) - 1)
    for i, x in enumerate(A):
        for j, y in enumerate(B):
            out[i + j] += x * y
    return out


def mu_inf_multi(P, Es):
    """Exact Int_{G_P} PROD_i m_{E_i}: product of per-cell linears."""
    pts = {F(0), F(1)} | P_breaks(P)
    for E in Es:
        pts |= cluster_breaks(E)
    pts = sorted(pts)
    total = F(0)
    for x, y in zip(pts, pts[1:]):
        mid = (x + y) / 2
        if not all(in_band_frac(p, mid) for p in P):
            continue
        ln = y - x
        poly = [F(1)]
        for E in Es:
            f1 = m_E(E, x + ln / 3)
            f2 = m_E(E, x + 2 * ln / 3)
            poly = poly_mul(poly, [2 * f1 - f2, 3 * (f2 - f1)])  # a + b s
        total += ln * sum(c / (j + 1) for j, c in enumerate(poly))
    return total


def m_joint_2cover(E2, E1, a):
    """b=2 ray fiber: gamma in [0,2); E2 arcs (2e*a, 1/7 halfwidth... width
    2/7 => halfwidth 1/7); E1 arcs at (e'a, 1/28*2=1/14 halfwidth) twice."""
    cw = [(((2 * e * a) % 2), F(1, 7)) for e in E2]
    for e in E1:
        c = (e * a) % 1
        cw.append((c, F(1, 14)))
        cw.append((c + 1, F(1, 14)))
    return m_arcs(cw, circ=2) / 2


def mu_inf_2cover(P, E2, E1):
    """Exact Int_{G_P} m_joint for the V1 = 2 V2 ray coupling."""
    slopes = [(2 * e, F(1, 7)) for e in E2] + [(e, F(1, 14)) for e in E1]
    pts = {F(0), F(1)} | P_breaks(P)
    gens = []
    for s, h in slopes:
        gens.append((s, -h))
        gens.append((s, h))
        if h == F(1, 14):  # E1 arcs live twice (offset +1)
            gens.append((s, 1 - h))
            gens.append((s, 1 + h))
    for (s1, w1), (s2, w2) in combinations(gens, 2):
        d = s1 - s2
        if d == 0:
            continue
        w = w2 - w1
        k0 = int(-(w) / d) - 3 if d else 0
        for k in range(k0, k0 + 3 * abs(d) + 7):
            alpha = (w + 2 * k) / d
            if 0 < alpha < 1:
                pts.add(F(alpha))
    pts = sorted(pts)
    total = F(0)
    for x, y in zip(pts, pts[1:]):
        mid = (x + y) / 2
        if not all(in_band_frac(p, mid) for p in P):
            continue
        ln = y - x
        f1 = m_joint_2cover(E2, E1, x + ln / 3)
        f2 = m_joint_2cover(E2, E1, x + 2 * ln / 3)
        total += ln * (f1 + f2) / 2
    return total


P = [1, 2, 3, 4, 5]
E2 = [0, 1]
E1 = [0, 1, 2, 3, 4, 5]
print("(I) SEPARATED SCALES: P={1..5}, E2={0,1}@V2, E1={0..5}@V1")
mi = mu_inf_multi(P, [E2, E1])
mP = mu_inf_multi(P, [])
floor = F(5, 7) * F(1, 7) * mP
print(f"    mu_inf (product) = {mi} = {float(mi):.6f}; m_P = {float(mP):.6f}; "
      f"floor (5/7)(1/7)m_P = {float(floor):.6f} "
      f"[{'OK' if mi >= floor else 'VIOLATION'}]")
print("    ratio study at fixed V2=60 (isolates the V2/V1 term):")
for V2, V1 in [(60, 240), (60, 600), (60, 1200), (60, 2400),
               (40, 1600), (100, 10000)]:
    S = sorted(P + [V2 - e for e in E2] + [V1 - e for e in E1])
    assert len(set(S)) == 13
    m = mu_exact(S)
    err = float(m - mi)
    pred = 1 / V2 + V2 / V1
    print(f"      V2={V2:4d} V1={V1:6d}: mu = {float(m):.6f} "
          f"err = {err:+.6f}  err/(1/V2+V2/V1) = {err / pred:+.3f}")

print("\n(II) RATIONAL RAY b=2 (V1 = 2*V2 EXACTLY -- the middle-speed coupling):")
mj = mu_inf_2cover(P, E2, E1)
print(f"    product limit (WRONG here) = {float(mi):.6f}; "
      f"2-cover joint limit = {mj} = {float(mj):.6f}")
for V2 in (300, 600, 1200):
    S = sorted(P + [V2 - e for e in E2] + [2 * V2 - e for e in E1])
    assert len(set(S)) == 13
    m = mu_exact(S)
    print(f"      V2={V2:5d}: mu = {float(m):.6f}  "
          f"err_vs_product = {float(m - mi):+.6f}  "
          f"err_vs_2cover = {float(m - mj):+.6f}  "
          f"|err2cover|*V2 = {abs(float(m - mj)) * V2:.2f}")

print("\n(III) three-scale sanity: P={1..3}, {0,1}@V3, {0,1}@V2, {0..5}@V1:")
P3 = [1, 2, 3]
Ea, Eb, Ec = [0, 1], [0, 1], [0, 1, 2, 3, 4, 5]
mi3 = mu_inf_multi(P3, [Ea, Eb, Ec])
mP3 = mu_inf_multi(P3, [])
print(f"    mu_inf = {float(mi3):.6f}; floor (5/7)^2(1/7)m_P = "
      f"{float(F(5,7)**2 * F(1,7) * mP3):.6f}")
for V3, V2, V1 in [(40, 400, 4000), (60, 900, 13500)]:
    S = sorted(P3 + [V3 - e for e in Ea] + [V2 - e for e in Eb]
               + [V1 - e for e in Ec])
    assert len(set(S)) == 13
    m = mu_exact(S)
    print(f"      ({V3},{V2},{V1}): mu = {float(m):.6f} "
          f"err = {float(m - mi3):+.6f}")
