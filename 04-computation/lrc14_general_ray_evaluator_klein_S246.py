#!/usr/bin/env python3
"""
THE GENERAL-(a,b) RAY EVALUATOR (klein-2026-07-10-S246, HYP-5955):
completing THM-688(B) -- the exact limit measure for rational-ray coupled
clusters at ANY ray a/b, gcd(a,b) = 1 (S237 implemented only b = 2).

Setup: top cluster {V - e : e in E_top} at scale V; RAY cluster
{V_mid - f : f in E_mid} with V_mid = (a*V - c)/b.  Writing V*alpha = j +
beta:  V_mid*alpha = (a(j + beta) - c*alpha)/b, so with s' = a*j mod b
(uniform over Z_b, gcd(a,b) = 1) the mid-cluster phase is

    phi = frac( (s' + a*beta)/b - (c/b)*alpha ),

and the (V_mid - f)-runner is safe iff phi - f*alpha is in the band, i.e.

    a*beta  in  (c + b*f)*alpha - s' + b*Z  +/-  b/14      (mod nothing:
    a*beta ranges over [0, a) -- enumerate the integer copies).

m_joint(alpha) = (1/b) Sum_{s'} meas{beta in [0,1): beta avoids the top
arcs (e*alpha +- 1/14)  AND  a*beta avoids every mid copy};
mu_inf = Int_{G_P} m_joint  (piecewise linear in alpha: top endpoints move
at rates e, mid beta-endpoints at rates (c + b*f)/a; trapezoid on interior
thirds of the breakpoint cells).

Verified here: (i) regression vs the S237 b=2 evaluator; (ii) new rays
(2,3), (1,3), (3,4) against exact finite-V sweeps; (iii) the merged floor;
(iv) a ray-class positivity census.
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd
import random

random.seed(246)
BAND_LO, BAND_HI = F(1, 14), F(13, 14)


def lcm(x, y):
    return x * y // gcd(x, y)


def in_band_frac(v, al):
    r = (v * al) % 1
    return BAND_LO <= r <= BAND_HI


def mu_exact(S):
    L = 1
    for v in S:
        L = lcm(L, v)
    D2 = 28 * L
    pts = {0, D2}
    for v in S:
        st = 2 * L // v
        for k in range(v):
            pts.add((14 * k + 1) * st)
            pts.add((14 * k + 13) * st)
    pts = sorted(pts)
    good = 0
    for x, y in zip(pts, pts[1:]):
        mch = x + y
        ok = True
        for v in S:
            t = v * mch % (2 * D2)
            if not (2 * D2 <= 14 * t <= 26 * D2):
                ok = False
                break
        if ok:
            good += y - x
    return F(good, D2)


def merge_meas(arcs):
    """Total length of a union of intervals within [0,1) (already clipped)."""
    arcs = sorted(a for a in arcs if a[0] < a[1])
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
    return tot


def clip_mod1(lo, hi):
    """Split (lo, hi) into pieces inside [0,1) modulo 1."""
    out = []
    span = hi - lo
    if span >= 1:
        return [(F(0), F(1))]
    lo = lo % 1
    hi = lo + span
    if hi <= 1:
        out.append((lo, hi))
    else:
        out.append((lo, F(1)))
        out.append((F(0), hi - 1))
    return out


def m_joint(E_top, E_mid, A, B, C, al):
    """Exact joint fiber measure at rational alpha for ray a/b = A/B, shift C."""
    total = F(0)
    top = []
    for e in E_top:
        for piece in clip_mod1(e * al - F(1, 14), e * al + F(1, 14)):
            top.append(piece)
    for sp in range(B):
        bad = list(top)
        for f in E_mid:
            center = (C + B * f) * al - sp
            # a*beta in (center - B/14, center + B/14) + B*k, a*beta in [0, A)
            lo0 = center - F(B, 14)
            k0 = int((-lo0) / B) - 2
            for k in range(k0, k0 + int(A / B) + 5):
                lo = lo0 + B * k
                hi = center + F(B, 14) + B * k
                # intersect with [0, A), then divide by A into beta-space
                l2, h2 = max(lo, F(0)), min(hi, F(A))
                if l2 < h2:
                    bad.append((l2 / A, h2 / A))
        total += 1 - merge_meas(bad)
    return total / B


def P_breaks(P):
    pts = set()
    for p in P:
        for k in range(p):
            pts.add(F(14 * k + 1, 14 * p))
            pts.add(F(14 * k + 13, 14 * p))
    return pts


def mu_inf_ray(P, E_top, E_mid, A, B, C, extra_denom=None):
    """Exact Int_{G_P} m_joint: breakpoints from all pairwise endpoint
    collisions (rates e for top, (C + Bf)/A for mid, per s' and copy)."""
    rates = [F(e) for e in E_top] + [F(C + B * f, A) for f in E_mid]
    pts = {F(0), F(1)} | P_breaks(P)
    # endpoint collisions: (r1 - r2)*alpha = rational offsets; enumerate on a
    # common grid: denominators divide 14*A*B*lcm(rate denominators) -- use a
    # conservative refinement grid from pairwise rate differences
    offs = [F(1, 14), F(-1, 14), F(B, 14 * A), F(-B, 14 * A)]
    for r1, r2 in combinations(set(rates), 2):
        d = r1 - r2
        if d == 0:
            continue
        for o1 in offs:
            for o2 in offs:
                base = (o2 - o1) / d
                step = abs(F(1) / d)
                k0 = int((0 - base) / step) - 2
                for k in range(k0, k0 + int(1 / step) + 5):
                    al = base + k * step
                    if 0 < al < 1:
                        pts.add(al)
    # also per-object wrap events (endpoint crossing 0/1): rate*alpha + off in Z
    for r in set(rates):
        if r == 0:
            continue
        for o in offs + [F(sp) for sp in range(B)]:
            step = abs(F(1) / r)
            base = (F(0) - o) / r
            k0 = int((0 - base) / step) - 2
            for k in range(k0, k0 + int(1 / step) + 5):
                al = base + k * step
                if 0 < al < 1:
                    pts.add(al)
    pts = sorted(pts)
    total = F(0)
    for x, y in zip(pts, pts[1:]):
        mid = (x + y) / 2
        if all(in_band_frac(p, mid) for p in P):
            ln = y - x
            f1 = m_joint(E_top, E_mid, A, B, C, x + ln / 3)
            f2 = m_joint(E_top, E_mid, A, B, C, x + 2 * ln / 3)
            total += ln * (f1 + f2) / 2
    return total


P = [1, 2, 3, 4, 5]
E_mid = [0, 1]          # the S237 demo had E2 = {0,1} as the mid cluster
E_top = [0, 1, 2, 3, 4, 5]

print("(i) REGRESSION vs S237's b = 2 evaluator (ray 1/2, c = 0):")
mi = mu_inf_ray(P, E_top, E_mid, 1, 2, 0)
print(f"    general evaluator (a,b,c) = (1,2,0): mu_inf = {mi} = {float(mi):.6f}")
print(f"    S237 2-cover value: 49379/470400 = {float(F(49379, 470400)):.6f} "
      f"-> match: {mi == F(49379, 470400)}")

print("\n(ii) NEW RAYS vs exact finite-V sweeps:")
for A, B, C in [(2, 3, 0), (1, 3, 0), (3, 4, 1)]:
    mi = mu_inf_ray(P, E_top, E_mid, A, B, C)
    print(f"    ray (a,b,c) = ({A},{B},{C}): mu_inf = {float(mi):.6f}")
    for Vtop in (600, 1200, 2400):
        if (A * Vtop - C) % B != 0:
            Vtop += B - ((A * Vtop - C) % B) * pow(A, -1, B) % B  # adjust
            while (A * Vtop - C) % B != 0:
                Vtop += 1
        Vmid = (A * Vtop - C) // B
        S = sorted(P + [Vmid - f for f in E_mid] + [Vtop - e for e in E_top])
        if len(set(S)) != 13:
            continue
        m = mu_exact(S)
        print(f"      V_top={Vtop:5d} (V_mid={Vmid}): mu = {float(m):.6f}  "
              f"err = {float(m - mi):+.6f}  |err|*V = {abs(float(m - mi)) * Vtop:.2f}")

print("\n(iii) merged floor m_joint >= 1 - (k_top + k_mid)/7 spot checks:")
floor = 1 - F(len(E_top) + len(E_mid), 7)
worst = None
for _ in range(60):
    al = F(random.randint(1, 999), 1000)
    v = m_joint(E_top, E_mid, 2, 3, 0, al)
    if worst is None or v < worst:
        worst = v
print(f"    min m_joint over 60 random alpha = {float(worst):.6f} "
      f">= floor {float(floor):.6f}: {worst >= floor}")

print("\n(iv) ray-class positivity census (mu_inf > 0):")
zeros = 0
cases = 0
for A, B in [(1, 2), (2, 3), (1, 3), (3, 4), (1, 4), (4, 5)]:
    for C in (0, 1):
        if gcd(A, B) != 1:
            continue
        for shapes in [([1, 2, 3, 4, 5], [0, 1], [0, 1, 2, 3, 4, 5]),
                       ([1, 2, 3], [0, 1, 2, 3], [0, 1, 2, 3, 4, 5]),
                       ([9, 10, 11, 12, 13], [0, 1], [0, 1, 2, 3, 4, 5])]:
            Pc, Em, Et = shapes
            mi = mu_inf_ray(Pc, Et, Em, A, B, C)
            cases += 1
            if mi <= 0:
                zeros += 1
                print(f"    ZERO at ray ({A},{B},{C}), P={Pc}")
print(f"    {cases} ray classes evaluated: zeros = {zeros}; all positive: "
      f"{zeros == 0}")
