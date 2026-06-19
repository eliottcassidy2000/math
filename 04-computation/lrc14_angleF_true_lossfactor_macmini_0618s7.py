#!/usr/bin/env python3
"""
lrc14_angleF_true_lossfactor_macmini_0618s7.py  (mac-mini-2026-06-18-S7)

ANGLE F part 6 — the TRUE loss factor of the absolute relation-lattice bound.

Part 4's Hbound (9000) was a grossly loose per-vector bound (it bounded
|K(n)| <= sum_t C(6,t)(t s7/pi)^supp / height, which explodes in supp).
Part 5 showed the EXACT per-vector |K(n)| summed is only O(1-2).  Here we compute
   AbsTrue(E) = sum_{0!=n in Lambda^o, truncated} |K(n)|   (EXACT kernel)
the best-possible absolute (triangle-inequality) bound, and the loss factor
   AbsTrue(E) / corr_true(E),
to report HONESTLY how far an |corr|-type covolume certificate is from closing.

We also report AbsTrue vs margin_8 = 0.35698:  if AbsTrue(AP) <= margin the
absolute certificate CLOSES; the data will show it does NOT (AbsTrue ~ 1.8 for AP)
but is within a single-digit factor -> a partial-sign-saving argument (keep the
signs of the SHORT vectors, bound only the support>=4 tail absolutely) is the
realistic finish.  We test that hybrid: signed short (supp<=3) + absolute tail.
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
        if len({int(((e * xm) % 1) * 7) for e in E}) == 7: total += x1 - x0
    return total
def M7(k): return sum(F((-1) ** t * math.comb(6, t)) * F(7 - t, 7) ** (k - 1) for t in range(7))
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
def covol_sat(nz):
    g = 0
    for x in nz: g = math.gcd(g, x)
    return math.sqrt(sum((x // g) ** 2 for x in nz))


def main():
    print("ANGLE F part 6: TRUE absolute loss factor & signed-short/absolute-tail hybrid")
    print("=" * 86)
    m7 = float(M7(8)); cap = float(F(2243, 5880)); margin = cap - m7
    print(f"k=8: M7={m7:.5f} cap={cap:.5f} margin={margin:.5f}")
    print(f"{'shape':<28}{'corr':>9}{'covol':>8}{'AbsTrue':>9}{'loss':>7}{'hybrid':>9}{'hyb<=marg':>10}")
    shapes = [
        ("AP {0..7}", list(range(8))),
        ("3AP{..6,8}", [0, 1, 2, 3, 4, 5, 6, 8]),
        ("2run40", [0, 1, 2, 3, 40, 41, 42, 43]),
        ("perf", [0, 2, 3, 4, 5, 6, 7, 9]),
        ("Sidon", [0, 1, 3, 7, 12, 20, 30, 44]),
        ("dissoc", [0, 1, 3, 7, 15, 31, 63, 127]),
    ]
    for name, E in shapes:
        nz = [e for e in E if e != 0]
        corr = float(measS7(E)) - m7
        cov = covol_sat(nz)
        abs_true = 0.0
        signed_short = 0.0
        abs_tail = 0.0
        for n in enum_lat(nz, 18, 4):
            if all(x == 0 for x in n): continue
            kv = Kk(n).real
            abs_true += abs(kv)
            supp = sum(1 for x in n if x != 0)
            if supp <= 3:
                signed_short += kv
            else:
                abs_tail += abs(kv)
        hybrid = abs(signed_short) + abs_tail
        loss = abs_true / corr if abs(corr) > 1e-6 else float("nan")
        ok = "YES" if hybrid <= margin else "no"
        print(f"{name:<28}{corr:>9.5f}{cov:>8.1f}{abs_true:>9.4f}{loss:>7.1f}{hybrid:>9.4f}{ok:>10}")
    print("-" * 86)
    print("AbsTrue = sum|K(n)| (exact); loss=AbsTrue/corr; hybrid=|signed supp<=3|+sum|supp>=4|.")
    print("If hybrid<=margin for ALL shapes incl AP, a SIGNED-SHORT certificate closes k=8.")
    print("DONE.")


if __name__ == "__main__":
    main()
