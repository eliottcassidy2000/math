#!/usr/bin/env python3
"""Kepler-polygon reformulation of the lonely runner conjecture (THM-4737 checks).
Loop = regular n-gon of perimeter 1 (n = number of runners including the stationary observer at a vertex), sides 1/n.
(A) Kepler II: constant speed along the sides sweeps area about the CENTRE at the constant rate a v / 2 (a = apothem);
    about the OBSERVER's vertex the rate on side s is d_s v / 2 with d_s = R (cos(pi/n) - cos((2s+1) pi/n)), zero exactly on
    the two incident sides s = 0, n-1; per lap the swept area about either point is the polygon area P.
(B) Kepler criterion: observer lonely at t (all ||v_i t|| >= 1/n) <=> no runner on an open incident side <=> all observer
    angular momenta positive (boundary: a runner exactly at the far vertex of an incident side has distance 1/n).
(C) Exact kappa(D) = max_t min_d ||t d|| over candidate times m/(2d), m/(d_i +- d_j); lonely measure m(D); minimal polygon
    order N_min(D) = least N with a lonely vertex time j/N (N = n iff no speed divisible by n); LRC check kappa >= 1/(k+1).
(D) Tight family {1..k}: N_min = n, m = 0."""
import math, itertools, random, sys
from fractions import Fraction as F
def polygon(n):
    R = 1 / (2 * n * math.sin(math.pi / n)); a = R * math.cos(math.pi / n); P = n * (1 / n) * a / 2
    return R, a, P
def point_on_polygon(n, x):
    """position at arc-length fraction x in [0,1) along the n-gon (vertex 0 at angle 0, counter-clockwise)."""
    R, a, P = polygon(n); s = int(math.floor(x * n)) % n; u = x * n - math.floor(x * n)
    A = (R * math.cos(2 * math.pi * s / n), R * math.sin(2 * math.pi * s / n)); B = (R * math.cos(2 * math.pi * (s + 1) / n), R * math.sin(2 * math.pi * (s + 1) / n))
    return (A[0] + u * (B[0] - A[0]), A[1] + u * (B[1] - A[1])), s
def area_rate(n, x, v, about):
    """instantaneous area sweep rate (1/2)|r x v| about the point 'about' for a runner at arc-fraction x with speed v (laps/time)."""
    R, a, P = polygon(n); s = int(math.floor(x * n)) % n
    A = (R * math.cos(2 * math.pi * s / n), R * math.sin(2 * math.pi * s / n)); B = (R * math.cos(2 * math.pi * (s + 1) / n), R * math.sin(2 * math.pi * (s + 1) / n))
    pos, _ = point_on_polygon(n, x); L = math.hypot(B[0] - A[0], B[1] - A[1]); vel = (v * (B[0] - A[0]) / L * 1.0, v * (B[1] - A[1]) / L)
    # speed along the perimeter is v (perimeter 1 per lap): |vel| = v * (side length)/(side length)... perimeter 1 so speed = v
    rx, ry = pos[0] - about[0], pos[1] - about[1]
    return 0.5 * abs(rx * vel[1] - ry * vel[0])
def dist(x):  # x Fraction mod 1 distance to Z
    x = x - (x.numerator // x.denominator); return min(x, 1 - x)
def kappa(D):
    cands = set()
    for d in D:
        for m in range(0, 2 * d + 1): cands.add(F(m, 2 * d))
    for a, b in itertools.combinations(D, 2):
        for s in (a + b, abs(a - b)):
            if s:
                for m in range(0, s + 1): cands.add(F(m, s))
    best = (F(0), None)
    for t in cands:
        val = min(dist(t * d) for d in D)
        if val > best[0]: best = (val, t)
    return best
def lonely_measure(D, n):
    # measure of {t in [0,1): all ||d t|| >= 1/n}: complement = union over d of intervals |d t - m| < 1/n
    bad = []
    for d in D:
        for m in range(0, d + 1):
            lo = (F(m) - F(1, n)) / d; hi = (F(m) + F(1, n)) / d; bad.append((max(lo, F(0)), min(hi, F(1))))
    bad.sort(); tot = F(0); cur_lo, cur_hi = None, None
    for lo, hi in bad:
        if lo >= hi: continue
        if cur_lo is None: cur_lo, cur_hi = lo, hi
        elif lo <= cur_hi: cur_hi = max(cur_hi, hi)
        else: tot += cur_hi - cur_lo; cur_lo, cur_hi = lo, hi
    if cur_lo is not None: tot += cur_hi - cur_lo
    return 1 - tot
def N_min(D, n, Nmax=400):
    for N in range(2, Nmax + 1):
        for j in range(1, N):
            if all(dist(F(d * j, N)) >= F(1, n) for d in D): return N, j
    return None, None

def main():
    print("(A) Kepler II on the regular n-gon (perimeter 1)")
    for n in (3, 4, 5, 6, 7, 10, 14):
        R, a, P = polygon(n); centre = (0.0, 0.0); obs = (R, 0.0)
        rates_c = [area_rate(n, x, 1.0, centre) for x in [i / 1000 + 1e-7 for i in range(1000)]]
        ds = [R * (math.cos(math.pi / n) - math.cos((2 * s + 1) * math.pi / n)) for s in range(n)]
        rates_o = [area_rate(n, (s + 0.37) / n, 1.0, obs) for s in range(n)]
        lap_c = sum(rates_c) / 1000; lap_o = sum(ds) * (1 / n) / 2
        print(f"  n={n}: centre rate a/2 = {a/2:.6f} (max dev {max(abs(r - a/2) for r in rates_c):.1e}); observer rates d_s/2 vs formula max dev {max(abs(rates_o[s] - ds[s]/2) for s in range(n)):.1e}; d_0 = {ds[0]:.1e}, d_(n-1) = {ds[-1]:.1e}; lap area centre {lap_c:.6f} observer {lap_o:.6f} P {P:.6f}")
    print("(B) Kepler criterion: lonely <=> all observer angular momenta > 0 (random speeds, random times)")
    random.seed(1); mism = 0; tested = 0
    for _ in range(20000):
        k = random.randint(2, 7); n = k + 1; D = random.sample(range(1, 60), k); t = random.random()
        lonely = all(min(v * t % 1, 1 - v * t % 1) >= 1 / n - 1e-12 for v in D)
        R, a, P = polygon(n); obs = (R, 0.0)
        momenta = [area_rate(n, (v * t) % 1, v, obs) for v in D]
        sides = [int(math.floor((v * t % 1) * n)) % n for v in D]
        crit = all(s not in (0, n - 1) for s in sides)
        # allow the measure-zero boundary (position exactly at a vertex)
        if lonely != crit and all(abs((v * t % 1) * n - round((v * t % 1) * n)) > 1e-9 for v in D): mism += 1
        tested += 1
    print(f"  tested {tested} random (D,t): mismatches between ||v t|| >= 1/n and 'no runner on an incident side' = {mism}")
    print("(C) exact kappa, lonely measure, minimal polygon order; all k-subsets of [1..M]")
    from collections import Counter
    for k, M in ((2, 20), (3, 16), (4, 13), (5, 11)):
        n = k + 1; cnt = 0; viol = 0; tight = []; nmin_hist = Counter(); zero_measure = 0; nmin_n = 0; divisible = 0
        for D in itertools.combinations(range(1, M + 1), k):
            if math.gcd(*D) != 1: continue
            cnt += 1; kap, t = kappa(D)
            if kap < F(1, n): viol += 1
            if kap == F(1, n): tight.append(D)
            mD = lonely_measure(D, n)
            if mD == 0: zero_measure += 1
            N, j = N_min(D, n); nmin_hist[N] += 1
            anydiv = any(d % n == 0 for d in D)
            if anydiv: divisible += 1
            if N == n: nmin_n += 1
            if (N == n) != (not anydiv): print("   !! N_min = n iff no speed divisible by n FAILS for", D, N, anydiv)
        print(f"  k={k}, M={M}: primitive sets {cnt}, kappa < 1/{n}: {viol}, tight (kappa = 1/{n}): {len(tight)} e.g. {tight[:6]}; zero lonely measure {zero_measure}; N_min = n: {nmin_n} (= sets with no speed divisible by n: {cnt - divisible}); N_min histogram {dict(sorted(nmin_hist.items()))}")
    print("(D) tight family {1..k}")
    for k in range(2, 8):
        n = k + 1; D = tuple(range(1, k + 1)); kap, t = kappa(D); print(f"  k={k}: kappa = {kap} at t = {t}, lonely measure = {lonely_measure(D, n)}, N_min = {N_min(D, n)}")

if __name__ == '__main__':
    main()
