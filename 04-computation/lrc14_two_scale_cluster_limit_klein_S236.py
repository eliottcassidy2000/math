#!/usr/bin/env python3
"""
THE TWO-SCALE CLUSTER LIMIT (klein-2026-07-10-S236, HYP-5905, THM-687).

For the unbounded residual families S_V = P u {V - e : e in E} (P subset
[1,13] the small speeds, E the fixed co-offsets, V -> infinity):

  frac((V-e)alpha) = frac(beta - e*alpha) with beta = frac(V*alpha), so the
  cluster runner is safe iff beta avoids the BAD ARC (e*alpha - 1/14,
  e*alpha + 1/14). As V -> infinity, (alpha, beta) equidistributes on T^2 and

    mu(S_V) -> mu_inf(P,E) := Int_{G_P} m_E(alpha) d alpha,
    m_E(alpha) = 1 - meas( Union_{e in E} (e*alpha - 1/14, e*alpha + 1/14) ),

  with rate C(P,E)/V (the two-scale transfer -- slice + crossing counting,
  the THM-685 toolbox one level up).

  THE k <= 6 FLOOR: k = |E| arcs of length 1/7 cannot cover the circle, so
  m_E >= 1 - k/7 pointwise and mu_inf >= (1 - k/7) * m_P; m_P > 0 is FREE
  from LRC(<=13)-strict (THM-667 Lemma A: m_P >= 1/(91 maxP)).  With THM-685
  (transfer at 14-coprime q) + kps LRCStrictRuler: THE WALL holds on every
  bounded-E slice with k <= 6 beyond an explicit V_0.

  k >= 7: mu_inf(P,E) > 0 is a FINITE exact criterion per class (this file's
  evaluator); the apex-7 dead zone = {alpha : the E-arcs cover} meeting G_P.

Everything exact (Fractions / integer-scaled sweeps; no floats in results).
The deep well {1..12, 182} is the V = 182 sample of (P,E) = ({1..12},{0}).
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations

HALF = F(1, 14)   # arc radius
BAND_LO, BAND_HI = F(1, 14), F(13, 14)


def lcm(a, b):
    return a * b // gcd(a, b)


def mu_exact(S):
    """Exact mu by integer-scaled sweep (S235 machinery)."""
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


def frac(x):
    return x - int(x) if x >= 0 else x - int(x) + (1 if x != int(x) else 0)


def m_E(E, a):
    """Exact fiber measure at rational alpha: 1 - |union of bad arcs|."""
    arcs = []
    for e in E:
        c = (e * a) % 1
        lo, hi = c - HALF, c + HALF
        if lo < 0:
            arcs.append((lo + 1, F(1)))
            arcs.append((F(0), hi))
        elif hi > 1:
            arcs.append((lo, F(1)))
            arcs.append((F(0), hi - 1))
        else:
            arcs.append((lo, hi))
    arcs.sort()
    tot = F(0)
    cur_lo, cur_hi = None, None
    for lo, hi in arcs:
        if cur_hi is None or lo > cur_hi:
            if cur_hi is not None:
                tot += cur_hi - cur_lo
            cur_lo, cur_hi = lo, hi
        else:
            cur_hi = max(cur_hi, hi)
    if cur_hi is not None:
        tot += cur_hi - cur_lo
    return 1 - tot


def in_band_frac(v, a):
    r = (v * a) % 1
    return BAND_LO <= r <= BAND_HI


def mu_inf(P, E):
    """Exact mu_inf(P,E) = Int_{G_P} m_E: cell decomposition + trapezoid
    (m_E is piecewise linear; cells cut at P-band crossings and E-pair arc
    topology changes frac(d*alpha) in {0, 1/7, 6/7})."""
    pts = {F(0), F(1)}
    for p in P:
        for k in range(p):
            pts.add(F(14 * k + 1, 14 * p))
            pts.add(F(14 * k + 13, 14 * p))
    ds = {abs(e - f) for e, f in combinations(E, 2) if e != f}
    for d in ds:
        for k in range(d):
            pts.add(F(k, d))
            pts.add(F(7 * k + 1, 7 * d))
            pts.add(F(7 * k + 6, 7 * d))
    pts = sorted(pts)
    total = F(0)
    comp_mP = F(0)
    for x, y in zip(pts, pts[1:]):
        mid = (x + y) / 2
        if all(in_band_frac(p, mid) for p in P):
            comp_mP += y - x
            p1 = x + (y - x) / 3
            p2 = x + 2 * (y - x) / 3
            total += (y - x) * (m_E(E, p1) + m_E(E, p2)) / 2
    return total, comp_mP


def study(P, E, name, Vs):
    k = len(E)
    mi, mP = mu_inf(P, E)
    floor = (1 - F(k, 7)) * mP
    print(f"\n{name}: P={P}, E={E} (k={k})")
    print(f"  m_P = {mP} = {float(mP):.6f}   mu_inf = {mi} = {float(mi):.6f}")
    print(f"  k<=6 floor (1-k/7)*m_P = {float(floor):.6f}: "
          f"mu_inf >= floor {'OK' if mi >= floor else 'VIOLATION'}")
    Cmax = 0
    for V in Vs:
        S = sorted(P + [V - e for e in E])
        assert len(set(S)) == 13, (name, V)
        m = mu_exact(S)
        err = m - mi
        Cmax = max(Cmax, abs(float(err)) * V)
        print(f"    V={V:5d}: mu(S_V) = {float(m):.6f}  "
              f"err = {float(err):+.6f}  |err|*V = {abs(float(err)) * V:7.2f}")
    print(f"  measured C = max|err|*V = {Cmax:.2f}  "
          f"(2(sum P + sum E) = {2 * (sum(P) + sum(E))})")
    return mi, mP


# (a) k=1, the deep-well family
study(list(range(1, 13)), [0], "k=1 deep-well family {1..12, V}",
      [50, 100, 182, 300, 500, 1000, 2000])

# (b) k=2
study(list(range(1, 12)), [0, 1], "k=2 family {1..11, V-1, V}",
      [100, 300, 1000, 2000])

# (c) k=6 boundary of the free zone
study(list(range(1, 8)), [0, 1, 2, 3, 4, 5], "k=6 family {1..7, V-5..V}",
      [100, 300, 1000])

# (d) k=7 apex demo -- the first class where the arcs CAN cover
study(list(range(1, 7)), [0, 1, 2, 3, 4, 5, 6], "k=7 APEX {1..6, V-6..V}",
      [100, 300, 1000])

# (e) k=7 with spread co-offsets (arcs spread out -- do they still choke?)
study(list(range(1, 7)), [0, 2, 5, 9, 14, 20, 27],
      "k=7 spread co-offsets {1..6, V-27..V}", [300, 1000])

print("\nREADING: (C) verified where mu_inf >= (1-k/7) m_P; the two-scale "
      "error is O(1/V) with measured C; k=7 classes show whether the apex "
      "dead-zone kills mu_inf or the fiber survives on G_P.")
