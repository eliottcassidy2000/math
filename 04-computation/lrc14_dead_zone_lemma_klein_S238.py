#!/usr/bin/env python3
"""
THE DEAD-ZONE LEMMA (klein-2026-07-10-S238, HYP-5915, THM-689).

Question (THM-687/688's named structural gap): is mu_inf(P,E) > 0 for EVERY
class -- equivalently, can the dead zone D(E) = {alpha : the k arcs
(e*alpha +- 1/14) COVER the circle} contain the P-safe set G_P?

(1) k = 7 RIGIDITY (proved; confirmed here): m_E(alpha) = 0 with k = 7 needs
    meas(union) = 1 = total mass => pairwise-disjoint open arcs => a PERFECT
    1/7-net => (e-e')alpha in (1/7)Z for all pairs => alpha in a finite set.
    So meas(D) = 0 at k = 7: mu_inf > 0 whenever m_P > 0. k <= 7 CLOSED.

(2) MEASURE CRITERION at k >= 8: mu_inf = 0 requires G_P (ess.) subset D(E),
    hence meas(D) >= m_P. Exact meas(D) evaluator (cells where the linear
    m_E vanishes identically); census over shapes; min m_P census per |P|.

(3) THE FIRST-WINDOW HANDLE (proved): W_P := [1/(14 pmin), 13/(14 pmax)] is
    ALWAYS inside G_P (for alpha in W_P: p*alpha in [p/(14 pmin),
    13p/(14 pmax)] subset [1/14, 13/14], no wrap), nonempty iff
    pmax <= 13*pmin -- degenerating EXACTLY at the ratio-13 tight locus.
    The lemma reduces to: can D(E) cover W_P?

(4) THE HUNT: mu_inf over a hard grid -- min-m_P P-shapes x big-D E-shapes
    (consecutive/AP/detuned/random/DESIGNED: centers ~ j/(7 alpha_0) equally
    spaced at alpha_0 in W_P) -- looking for any zero (the wall's boundary).
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd
import random

random.seed(238)
BAND_LO, BAND_HI = F(1, 14), F(13, 14)


def in_band_frac(v, a):
    r = (v * a) % 1
    return BAND_LO <= r <= BAND_HI


def m_arcs(cw, circ=1):
    arcs = []
    for c, h in cw:
        c = c % circ
        lo, hi = c - h, c + h
        if lo < 0:
            arcs += [(lo + circ, F(circ)), (F(0), hi)]
        elif hi > circ:
            arcs += [(lo, F(circ)), (F(0), hi - circ)]
        else:
            arcs.append((lo, hi))
    arcs.sort()
    tot = F(0)
    cl = ch = None
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


def P_breaks(P):
    pts = set()
    for p in P:
        for k in range(p):
            pts.add(F(14 * k + 1, 14 * p))
            pts.add(F(14 * k + 13, 14 * p))
    return pts


def E_breaks(E):
    pts = set()
    for e, f in combinations(E, 2):
        d = abs(e - f)
        if d:
            for k in range(d):
                pts.add(F(k, d))
                pts.add(F(7 * k + 1, 7 * d))
                pts.add(F(7 * k + 6, 7 * d))
    return pts


def cells(E, extra=()):
    pts = {F(0), F(1)} | E_breaks(E)
    for x in extra:
        pts |= x
    return sorted(pts)


def meas_D(E, window=None):
    """Exact measure of the dead zone (m_E == 0 cells), optionally within a
    window [a,b]; also returns the dead interval list."""
    pts = cells(E)
    if window:
        pts = sorted(set(pts) | {window[0], window[1]})
    tot = F(0)
    dead = []
    for x, y in zip(pts, pts[1:]):
        if window and (y <= window[0] or x >= window[1]):
            continue
        ln = y - x
        if m_E(E, x + ln / 3) == 0 and m_E(E, x + 2 * ln / 3) == 0:
            tot += ln
            dead.append((x, y))
    return tot, dead


def mu_inf(P, E):
    pts = cells(E, [P_breaks(P)])
    total = F(0)
    mP = F(0)
    for x, y in zip(pts, pts[1:]):
        mid = (x + y) / 2
        if all(in_band_frac(p, mid) for p in P):
            ln = y - x
            mP += ln
            total += ln * (m_E(E, x + ln / 3) + m_E(E, x + 2 * ln / 3)) / 2
    return total, mP


print("(1) k=7 RIGIDITY: meas(D) for k=7 shapes (must all be 0):")
for name, E in [("consec", list(range(7))), ("AP2", [0, 2, 4, 6, 8, 10, 12]),
                ("AP3", [0, 3, 6, 9, 12, 15, 18]),
                ("rand", sorted(random.sample(range(25), 7)))]:
    d, _ = meas_D(E)
    print(f"    {name:8s} E={E}: meas(D) = {d}")

print("\n(2) meas(D) census at k = 8..13 (consecutive/AP/detuned/random):")
maxD = {}
for k in range(8, 14):
    shapes = [("consec", list(range(k))), ("AP2", [2 * j for j in range(k)]),
              ("AP3", [3 * j for j in range(k)]),
              ("detun", list(range(k - 1)) + [k + 3])]
    for i in range(6):
        shapes.append((f"rnd{i}", sorted(random.sample(range(0, 3 * k), k))))
    mx = F(0)
    for name, E in shapes:
        d, _ = meas_D(E)
        mx = max(mx, d)
        if name in ("consec", "AP2"):
            print(f"    k={k} {name:6s}: meas(D) = {float(d):.6f}")
    maxD[k] = mx
    print(f"    k={k} MAX over {len(shapes)} shapes: {float(mx):.6f}")

print("\n    min m_P per |P| (exact, exhaustive over P subset [1,13]):")
minmP = {}
for size in (5, 4, 3, 2, 1):
    best = None
    for P in combinations(range(1, 14), size):
        _, mP = mu_inf(list(P), [0])
        if best is None or mP < best[0]:
            best = (mP, P)
    minmP[size] = best
    k = 13 - size
    verdict = ("criterion CLOSES k=%d" % k) if maxD.get(k, F(1)) < best[0] \
        else ("criterion OPEN at k=%d (containment needed)" % k)
    print(f"    |P|={size}: min m_P = {float(best[0]):.6f} at {best[1]}; "
          f"max meas(D) = {float(maxD.get(k, F(0))):.6f} -> {verdict}")

print("\n(3)+(4) THE HUNT: mu_inf over hard (P,E) grid (min-m_P P's x "
      "adversarial E's incl. window-designed):")
hunt_min = None
zeros = []
for size in (5, 4, 3, 2, 1):
    k = 13 - size
    Ps = [list(minmP[size][1])]
    # ratio-13 (degenerate-window) P's and near-AP P's of this size
    if size >= 2:
        Ps.append([1] + list(range(15 - size, 14)))     # {1, high block} ratio 13
    Ps.append(list(range(1, size + 1)))                  # {1..size}
    for P in Ps:
        pmin, pmax = min(P), max(P)
        W = (F(1, 14 * pmin), F(13, 14 * pmax))
        Es = [list(range(k)), [2 * j for j in range(k)],
              list(range(k - 1)) + [k + 3]]
        if W[1] > W[0]:
            for t in (F(1, 3), F(1, 2), F(2, 3)):
                a0 = W[0] + t * (W[1] - W[0])
                step = 1 / (7 * a0)
                Ed = sorted({int(j * step + F(1, 2)) for j in range(k)})
                while len(Ed) < k:
                    Ed.append(max(Ed) + 1)
                Es.append(Ed)
        for i in range(4):
            Es.append(sorted(random.sample(range(0, 40), k)))
        for E in Es:
            mi, mP = mu_inf(P, E)
            if mi == 0:
                zeros.append((P, E))
            if hunt_min is None or mi < hunt_min[0]:
                hunt_min = (mi, P, E)
        # window coverage check for the designed adversaries
        if W[1] > W[0]:
            dW, _ = meas_D(Es[3], window=W)
            print(f"    |P|={size} P={P}: W_P len = {float(W[1] - W[0]):.5f}; "
                  f"designed-E dead mass in W_P = {float(dW):.6f}")
print(f"\n    HUNT RESULT: {len(zeros)} zero classes found; "
      f"min mu_inf = {float(hunt_min[0]):.6f} at P={hunt_min[1]}, E={hunt_min[2]}")
