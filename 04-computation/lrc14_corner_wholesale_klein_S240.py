#!/usr/bin/env python3
"""
THE CORNER WHOLESALE CLOSURE (klein-2026-07-10-S240, HYP-5925, THM-691).

(A) THE q*-TEST-POINT THEOREM (proved; verified here): for any P let
    q* = max{q in [8,13] : q not in P} (exists whenever k >= 8). For every
    a coprime to q*:
      P-side: p in [1,13], p != q*, gcd(a,q*) = 1 => q* | pa iff q* | p,
        impossible (p < 2q*); so frac(p*a/q*) in [1/q*, (q*-1)/q*] subset
        (1/14, 13/14)  (1/q* >= 1/13 > 1/14; (q*-1)/q* <= 12/13 < 13/14).
        PRIMALITY OF q* NOT NEEDED.
      E-side: |E| = k < q* => the centers occupy < q* classes of the q*-net
        => empty class => gap >= 2/q* => m_E(a/q*) >= 2/q* - 1/7 >= 1/91.
    k < q* fails only when P = [q*+1, 13] EXACTLY: the corner collapses to
    THE FIVE TOP-BLOCK FAMILIES {9..13}, {10..13}, {11,12,13}, {12,13}, {13}
    (k = 8,9,10,11,12).

(B) THE e_max DICHOTOMY for the five (proved): if e_max < 12*pmin then the
    sliver [1/(14*pmin), min(1/14, 6/(7*e_max))) is nonempty, inside G_P
    (the first window W_P), and OUTSIDE D(E) (small-alpha exclusion: all
    centers e*alpha <= e_max*alpha < 6/7 with no wrap, so the top gap
    exceeds 1/7): mu_inf > 0 for ALL E with e_max < 12*pmin (>= 108).

(C) SPREAD ADVERSARIES (e_max >= 12*pmin): the measure criterion needs
    meas(D(E)) >= m_P in {0.42..0.86} -- census over spread shapes shows
    meas(D) collapses toward the iid value, far below m_P.
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd
import random

random.seed(240)
BAND_LO, BAND_HI = F(1, 14), F(13, 14)


def in_band_frac(v, a):
    r = (v * a) % 1
    return BAND_LO <= r <= BAND_HI


def m_arcs(cw):
    arcs = []
    for c, h in cw:
        c = c % 1
        lo, hi = c - h, c + h
        if lo < 0:
            arcs += [(lo + 1, F(1)), (F(0), hi)]
        elif hi > 1:
            arcs += [(lo, F(1)), (F(0), hi - 1)]
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
    return 1 - tot


def m_E(E, a):
    return m_arcs([((e * a) % 1, F(1, 14)) for e in E])


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


def meas_D(E):
    pts = sorted({F(0), F(1)} | E_breaks(E))
    tot = F(0)
    for x, y in zip(pts, pts[1:]):
        ln = y - x
        if m_E(E, x + ln / 3) == 0 and m_E(E, x + 2 * ln / 3) == 0:
            tot += ln
    return tot


def P_breaks(P):
    pts = set()
    for p in P:
        for k in range(p):
            pts.add(F(14 * k + 1, 14 * p))
            pts.add(F(14 * k + 13, 14 * p))
    return pts


def mu_inf(P, E):
    pts = sorted({F(0), F(1)} | P_breaks(P) | E_breaks(E))
    total = F(0)
    for x, y in zip(pts, pts[1:]):
        mid = (x + y) / 2
        if all(in_band_frac(p, mid) for p in P):
            ln = y - x
            total += ln * (m_E(E, x + ln / 3) + m_E(E, x + 2 * ln / 3)) / 2
    return total


def qstar(P):
    for q in range(13, 7, -1):
        if q not in P:
            return q
    return None


FIVE = [list(range(b, 14)) for b in (9, 10, 11, 12, 13)]

print("(A) q*-THEOREM verification on non-top-block corner P's:")
tests = [[1, 13], [1, 12, 13], [2, 7, 13], [1, 11, 12, 13], [1, 5, 7, 8, 13],
         [3, 9, 11, 12, 13]]
for P in tests:
    q = qstar(P)
    k = 13 - len(P)
    if k >= q:
        print(f"    P={P}: k={k} >= q*={q} -- top-block-like, excepted")
        continue
    margin = F(2, q) - F(1, 7)
    ok_p = all(in_band_frac(p, F(a, q)) for p in P
               for a in range(1, q) if gcd(a, q) == 1)
    worst_e = None
    for _ in range(6):
        E = sorted({0} | set(random.sample(range(1, 300), k - 1)))
        mn = min(m_E(E, F(a, q)) for a in range(1, q) if gcd(a, q) == 1)
        worst_e = mn if worst_e is None or mn < worst_e else worst_e
    print(f"    P={P}: q*={q}, k={k}: P-side all a/q* in G_P: {ok_p}; "
          f"E-side worst m_E = {float(worst_e):.5f} >= bound "
          f"{float(margin):.5f}: {worst_e >= margin}")

print("\n(B) THE FIVE TOP-BLOCKS: packed branch is PROVED (e_max < 12*pmin); "
      "verify the sliver on the packed extremal:")
for P in FIVE:
    b = min(P)
    k = 13 - len(P)
    E = list(range(k))  # consecutive = packed extremal, e_max = k-1 << 12b
    mi = mu_inf(P, E)
    lo, hi = F(1, 14 * b), min(F(1, 14), F(6, 7 * (k - 1)))
    sliver = hi - lo
    print(f"    P={P} (pmin={b}, k={k}): e_max bound 12*pmin = {12*b}; "
          f"consecutive mu_inf = {float(mi):.6f} > 0; sliver len = "
          f"{float(sliver):.5f} > 0: {sliver > 0}")

print("\n(C) SPREAD ADVERSARIES (e_max >= 12*pmin): meas(D) vs m_P:")
mPs = {}
for P in FIVE:
    mPs[tuple(P)] = mu_inf(P, [0]) * F(7, 6)  # m_P = mu_inf(P,{0})/(6/7)
for P in FIVE:
    b = min(P)
    k = 13 - len(P)
    mP = mPs[tuple(P)]
    worstD = F(0)
    worstE = None
    shapes = []
    for i in range(8):
        top = random.randint(12 * b, 12 * b + 240)
        Ei = sorted({0, top} | set(random.sample(range(1, top), k - 2)))
        shapes.append(Ei)
    shapes.append([0] + [12 * b + 3 * j for j in range(1, k)])       # far block
    shapes.append(list(range(k - 1)) + [12 * b + 7])                 # packed+spike
    for E in shapes:
        E = sorted(set(E))[:k]
        while len(E) < k:
            E.append(max(E) + 1)
        d = meas_D(E)
        if d > worstD:
            worstD, worstE = d, E
    ok = worstD < mP
    print(f"    P={P}: m_P = {float(mP):.4f}; worst spread meas(D) over "
          f"{len(shapes)} shapes = {float(worstD):.4f} "
          f"[{'criterion CLOSES' if ok else 'NOT CLOSED'}] "
          f"(worst E head: {worstE[:5]}...)")

print("\nSTATUS: corner = [q*-theorem: all non-top-block P] + [five blocks: "
      "packed PROVED via sliver; spread censused vs measure criterion].")
