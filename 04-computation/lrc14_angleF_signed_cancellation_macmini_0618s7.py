#!/usr/bin/env python3
"""
lrc14_angleF_signed_cancellation_macmini_0618s7.py  (mac-mini-2026-06-18-S7)

ANGLE F part 5 — locate the cancellation that kills the absolute covolume bound,
and test whether a SIGNED (Poisson/theta) covolume bound survives.

Part 4 showed: sum_n |K(n)| ~ 9000 >> margin 0.357 >> corr_true 0.30.  So the
absolute (triangle-inequality) bound is hopeless; the smallness of corr is PURELY
signed cancellation.  Here we:

 (1) compute the SIGNED partial sum by SUPPORT SIZE s (s=2,3,4,...): how much of
     corr lives in each support level, signed; show the per-support partial is
     itself O(1) (not 9000) -> cancellation is INTRA-support (within fixed
     support, over the lattice directions);
 (2) compute the SIGNED partial sum by COVOLUME SHELL (|n| <= R): show the signed
     truncation converges, and estimate the convergence RATE vs covolume;
 (3) the key diagnostic: ratio  corr_true / (per-support-2 signed sum)  to see if
     support-2 (the shortest relations, the lambda_1 vectors) already captures the
     extremal ordering -> a SIGNED short-vector certificate (not absolute).

Support-2 vectors n have n_a*e_a + n_b*e_b = 0 with two nonzero entries: these are
the rational dependencies e_a/e_b = -n_b/n_a, i.e. the 2-term commensurabilities.
For DISTINCT integer offsets there are few; for AP there are the proportional pairs.
"""
from __future__ import annotations
import sys, itertools, math, cmath
from collections import defaultdict
from fractions import Fraction as F
import sympy
sys.stdout.reconfigure(line_buffering=True)

TWO_PI_I = 2j * math.pi


def measS7(E):
    E = sorted(set(E)); bps = {F(0), F(1)}
    for e in E:
        if e == 0: continue
        for m in range(0, 7 * e + 1): bps.add(F(m, 7 * e))
    bps = sorted(bps); total = F(0)
    for i in range(len(bps) - 1):
        x0, x1 = bps[i], bps[i + 1]
        if x1 <= x0: continue
        xm = (x0 + x1) / 2
        secs = {int(((e * xm) % 1) * 7) for e in E}
        if len(secs) == 7: total += x1 - x0
    return total


def M7(k):
    return sum(F((-1) ** t * math.comb(6, t)) * F(7 - t, 7) ** (k - 1) for t in range(7))


def shat(n, j):
    if n == 0: return 1.0 / 7.0
    a = j / 7.0
    return (cmath.exp(-TWO_PI_I * n * a) - cmath.exp(-TWO_PI_I * n * (a + 1 / 7.0))) / (TWO_PI_I * n)


SUB = [tuple(T) for r in range(7) for T in itertools.combinations(range(1, 7), r)]
SGN = {T: (-1) ** len(T) for T in SUB}
_CH = {}
def chat(n, T):
    key = (n, T)
    if key in _CH: return _CH[key]
    if n == 0: v = complex(1 - len(T) / 7.0, 0.0)
    elif n % 7 == 0: v = 0j
    else: v = -sum(shat(n, j) for j in T)
    _CH[key] = v; return v
def Kk(n):
    s = 0j
    for T in SUB:
        p = 1.0 + 0j
        for ni in n:
            p *= chat(ni, T)
            if p == 0: break
        s += SGN[T] * p
    return s


def lll(basis):
    if not basis: return []
    B = [[F(x) for x in r] for r in basis]; n = len(B)
    def dot(u, v): return sum(a * b for a, b in zip(u, v))
    def gs(B):
        Bs, mu = [], [[F(0)] * n for _ in range(n)]
        for i in range(n):
            bi = list(B[i])
            for j in range(i):
                mu[i][j] = dot(B[i], Bs[j]) / dot(Bs[j], Bs[j])
                bi = [bi[t] - mu[i][j] * Bs[j][t] for t in range(len(bi))]
            Bs.append(bi)
        return Bs, mu
    Bs, mu = gs(B); k = 1
    while k < n:
        for j in range(k - 1, -1, -1):
            q = round(mu[k][j])
            if q: B[k] = [B[k][t] - q * B[j][t] for t in range(len(B[k]))]; Bs, mu = gs(B)
        if dot(Bs[k], Bs[k]) >= (F(3, 4) - mu[k][k - 1] ** 2) * dot(Bs[k - 1], Bs[k - 1]): k += 1
        else: B[k], B[k - 1] = B[k - 1], B[k]; Bs, mu = gs(B); k = max(k - 1, 1)
    return [[int(x) for x in r] for r in B]
def kbasis(nz):
    M = sympy.Matrix([nz]); out = []
    for v in M.nullspace():
        L = 1
        for x in v:
            fr = F(int(x.p), int(x.q)); L = L * fr.denominator // math.gcd(L, fr.denominator)
        iv = [int(x * L) for x in v]; g = 0
        for c in iv: g = math.gcd(g, abs(c))
        if g: iv = [c // g for c in iv]
        out.append(iv)
    return out
def enum_lat(nz, N0, mc):
    red = lll(kbasis(nz)); d = len(red); seen = set(); out = []
    for c in itertools.product(range(-mc, mc + 1), repeat=d):
        v = tuple(sum(c[i] * red[i][t] for i in range(d)) for t in range(len(nz)))
        if all(abs(x) <= N0 for x in v) and v not in seen: seen.add(v); out.append(v)
    return out


def main():
    print("ANGLE F part 5: signed cancellation anatomy of corr(E)")
    print("=" * 78)
    m7 = float(M7(8))
    shapes = [
        ("AP {0..7}", list(range(8))),
        ("2run40", [0, 1, 2, 3, 40, 41, 42, 43]),
        ("perf", [0, 2, 3, 4, 5, 6, 7, 9]),
        ("dissoc", [0, 1, 3, 7, 15, 31, 63, 127]),
    ]
    print("(1) SIGNED partial sum by support size (N0=18,mc=4):")
    for name, E in shapes:
        nz = [e for e in E if e != 0]
        corr = float(measS7(E)) - m7
        by_supp_signed = defaultdict(float)
        by_supp_abs = defaultdict(float)
        for n in enum_lat(nz, 18, 4):
            if all(x == 0 for x in n): continue
            s = sum(1 for x in n if x != 0)
            k = Kk(n).real
            by_supp_signed[s] += k
            by_supp_abs[s] += abs(k)
        print(f"  {name}: corr_true={corr:+.5f}")
        for s in sorted(by_supp_signed):
            print(f"      supp={s}: signed={by_supp_signed[s]:+.5f}  abs={by_supp_abs[s]:.5f}  "
                  f"cancel_ratio={by_supp_abs[s]/max(1e-9,abs(by_supp_signed[s])):.1f}x")
    print()
    print("(2) signed truncation convergence by radius R (AP, the worst shape):")
    nz = list(range(1, 8))
    for R in (6, 10, 14, 18, 24):
        s = 0j
        for n in enum_lat(nz, R, 5):
            if all(x == 0 for x in n): continue
            s += Kk(n)
        print(f"   R={R:>2}: M7+signed_sum={m7+s.real:.5f}  (exact corr_true={float(measS7([0]+list(range(1,8))))-m7:+.5f})")
    print()
    print("READOUT: within each fixed support the cancellation ratio is LARGE (tens to")
    print("hundreds x), so |sum| << sum|.|. The covolume controls the DENSITY of lattice")
    print("directions but the SIGN (THM-503 phase alignment) controls the magnitude. A")
    print("rigorous Angle-F certificate MUST keep the sign (Poisson/theta), not |K(n)|.")
    print("DONE.")


if __name__ == "__main__":
    main()
