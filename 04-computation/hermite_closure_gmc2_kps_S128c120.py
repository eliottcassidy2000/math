#!/usr/bin/env python3
"""hermite_closure_gmc2_kps_S128c120.py -- kind-pasteur-2026-07-20-S128c120

THE {-1,0,1} / M=1 STRATUM CLOSES, AND THE OBJECT IS HERMITE.

THM-1585 showed the Gamma-domination bridge is false: the mass of E_r[psi_m] sits at an
INTERIOR k, not at the top.  That is a saddle, and the right tool for a saddle is not an
asymptotic estimate -- it is the orthogonal polynomial that owns it.

STEP 1 (Lagrange-Buermann, exact).  psi_m = [t^m] log(phi(u(t))/phi(0)) where u = t phi(u)
and phi(u) = R(r,u) = g_{-1} + g_0 u + g_1 u^2.  With H = log phi - log phi(0),
H' phi^m = phi' phi^{m-1} = (phi^m)'/m, so
    [t^m]H(u) = (1/m)[u^{m-1}]{H'(u) phi(u)^m} = (1/m^2) * m * [u^m] phi^m,
                                   ==>   psi_m = (1/m) [u^m] phi(u)^m.
Expanding phi = rho a + b u + rho c u^2, the u^m coefficient forces #(rho a) = #(rho c),
so every surviving term carries rho^{2k} = r^k -- RHO-FREE, and
    psi_m = (1/m) sum_k  m!/(k!^2 (m-2k)!) * W^k * B^{m-2k},   B = b,  W = r a c.

STEP 2 (the Gamma average, constant coefficients).  E_r[r^k] = k! cancels one k!:
    E_r[psi_m] = (1/m) sum_k  m!/(k! (m-2k)!) * w^k * b^{m-2k},   w = a c.
Write T_m(b,w) for that sum.  Setting w = -s^2/2, i.e. s = sqrt(-2w):
    T_m(b,w) = s^m * He_m(b/s),        He = PROBABILISTS' HERMITE,
because He_m(x) = sum_k m!/(k!(m-2k)!) (-1/2)^k x^{m-2k} exactly.

STEP 3 (the closure).  E_r[psi_m] = 0  <=>  He_m(b/s) = 0.
  * w = ac > 0: s is imaginary, so b/s is purely imaginary; every root of He_m is REAL,
    so He_m(b/s) != 0 unless b = 0, and at b = 0, He_m(0) != 0 for EVEN m.  Never all-zero.
  * w = ac < 0: s real, b/s real.  He_m(b/s) can vanish for some m -- but NOT for all m,
    because CONSECUTIVE HERMITE POLYNOMIALS HAVE NO COMMON ROOT: from the three-term
    recurrence He_{m+1} = x He_m - m He_{m-1}, a common root of He_m and He_{m+1} forces
    He_{m-1}(x) = 0, and downward induction reaches He_0 = 1 = 0, absurd.
So for ac != 0 there is ALWAYS an m with E_r[psi_m] != 0: such P is NOT in the nullcone.
And ac != 0 is exactly "both extreme charges present" = TWO-SIDED.  That is the one-sided
conjecture on this stratum, PROVED, with no asymptotics and no domination.

SCOPE, STATED UP FRONT: step 2 uses E_r[W^k B^{m-2k}] = k! w^k b^{m-2k}, which needs
a, b, c CONSTANT.  Non-constant a(r), b(r), c(r) are NOT covered by the Hermite closure.
"""
import sys
from fractions import Fraction as Fr
import numpy as np

MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 16


def fact(n):
    r = 1
    for i in range(2, n + 1):
        r *= i
    return r


def pmul(A, B):
    if not A or not B:
        return []
    C = [Fr(0)] * (len(A) + len(B) - 1)
    for i, a in enumerate(A):
        if a:
            for j, b in enumerate(B):
                C[i + j] += a * b
    while C and C[-1] == 0:
        C.pop()
    return C


def padd(A, B):
    C = [Fr(0)] * max(len(A), len(B))
    for i, a in enumerate(A):
        C[i] += a
    for i, b in enumerate(B):
        C[i] += b
    while C and C[-1] == 0:
        C.pop()
    return C


def series_v(B, W, N):
    v = [[Fr(1)]] + [[] for _ in range(N)]
    for _ in range(N + 1):
        v2 = [[] for _ in range(N + 1)]
        for i, vi in enumerate(v):
            if vi:
                for j, vj in enumerate(v):
                    if i + j <= N and vj:
                        v2[i + j] = padd(v2[i + j], pmul(vi, vj))
        new = [[Fr(1)]] + [[] for _ in range(N)]
        for k in range(N + 1):
            if k >= 1 and v[k - 1]:
                new[k] = padd(new[k], pmul(B, v[k - 1]))
            if k >= 2 and v2[k - 2]:
                new[k] = padd(new[k], pmul(W, v2[k - 2]))
        v = new
    return v


def series_log(v, N):
    w = [padd(v[i], [Fr(-1)] if i == 0 else []) for i in range(N + 1)]
    out = [[] for _ in range(N + 1)]
    powk = [[Fr(1)]] + [[] for _ in range(N)]
    for k in range(1, N + 1):
        nxt = [[] for _ in range(N + 1)]
        for i, pi in enumerate(powk):
            if pi:
                for j, wj in enumerate(w):
                    if i + j <= N and wj:
                        nxt[i + j] = padd(nxt[i + j], pmul(pi, wj))
        powk = nxt
        s = Fr((-1) ** (k + 1), k)
        for m in range(N + 1):
            if powk[m]:
                out[m] = padd(out[m], [x * s for x in powk[m]])
    return out


def Er(A):
    return sum(c * fact(k) for k, c in enumerate(A))


def psi_closed(a, b, c, m):
    """(1/m) sum_k m!/(k!^2 (m-2k)!) W^k B^{m-2k} as a poly in r; W = r*a*c, B = b."""
    out = [Fr(0)] * (m // 2 + 1)
    for k in range(m // 2 + 1):
        out[k] += (Fr(fact(m), fact(k) ** 2 * fact(m - 2 * k))
                   * Fr(a * c) ** k * Fr(b) ** (m - 2 * k)) / m
    while out and out[-1] == 0:
        out.pop()
    return out


def T(m, b, w):
    return sum(Fr(fact(m), fact(k) * fact(m - 2 * k)) * Fr(w) ** k * Fr(b) ** (m - 2 * k)
               for k in range(m // 2 + 1))


def He(m, x):
    """Probabilists Hermite via the three-term recurrence; exact on Fraction input."""
    h0, h1 = Fr(1), x
    if m == 0:
        return h0
    for k in range(1, m):
        h0, h1 = h1, x * h1 - k * h0
    return h1


print("=" * 88)
print("STEP 1 CHECK -- Lagrange-Buermann closed form vs the log-series computation")
print("=" * 88)
ok_all = True
for (a, b, c) in [(1, 1, 1), (1, 0, 1), (1, 3, 1), (1, 1, -1), (2, 5, -3), (1, -2, 4)]:
    v = series_v([Fr(b)] if b else [], [Fr(0), Fr(a * c)], MMAX)
    lg = series_log(v, MMAX)
    agree = all(lg[m] == psi_closed(a, b, c, m) for m in range(1, MMAX + 1))
    ok_all &= agree
    print("  (a,b,c) = %-14s  psi_m agrees for m = 1..%d : %s" % (str((a, b, c)), MMAX, agree))
print("  ALL AGREE: %s  -> psi_m = (1/m)[u^m] phi(u)^m confirmed exactly." % ok_all)

print()
print("=" * 88)
print("STEP 2 CHECK -- m * E_r[psi_m] = T_m(b, w) with w = ac, and T_m(b,-1/2) = He_m(b)")
print("=" * 88)
allok = True
for (a, b, c) in [(1, 1, 1), (1, 3, 1), (1, 1, -1), (2, 5, -3), (1, -2, 4), (1, 0, 1)]:
    w = a * c
    for m in (2, 5, 8, 11):
        lhs = m * Er(psi_closed(a, b, c, m))
        allok &= (lhs == T(m, b, w))
    print("  (a,b,c) = %-14s w = %-4d  m*E_r[psi_m] == T_m(b,w) for m = 2,5,8,11 : %s"
          % (str((a, b, c)), w, all(m * Er(psi_closed(a, b, c, m)) == T(m, b, w)
                                    for m in (2, 5, 8, 11))))
print("  ALL AGREE: %s" % allok)
print()
for b in (Fr(1), Fr(3), Fr(-2), Fr(1, 2)):
    print("     b = %-6s  T_m(b,-1/2) == He_m(b) for m = 1..8 : %s"
          % (b, all(T(m, b, Fr(-1, 2)) == He(m, b) for m in range(1, 9))))

print()
print("=" * 88)
print("STEP 3 CHECK -- consecutive Hermite polynomials share no root (the closure)")
print("=" * 88)
for m in range(2, 12):
    rm = np.polynomial.hermite_e.hermeroots(np.array([0] * m + [1], dtype=float))
    rm1 = np.polynomial.hermite_e.hermeroots(np.array([0] * (m + 1) + [1], dtype=float))
    mind = min(abs(x - y) for x in rm for y in rm1)
    print("  m = %-3d  min |root(He_m) - root(He_{m+1})| = %.6f" % (m, mind))

print()
print("=" * 88)
print("CONSEQUENCE -- the one-sided conjecture on this stratum (a, b, c CONSTANT)")
print("=" * 88)
tested, bad, worst = 0, [], 0
for a in range(-3, 4):
    for b in range(-3, 4):
        for c in range(-3, 4):
            if a * c == 0:
                continue
            tested += 1
            first = None
            for m in range(1, 25):
                if Er(psi_closed(a, b, c, m)) != 0:
                    first = m
                    break
            if first is None:
                bad.append((a, b, c))
            else:
                worst = max(worst, first)
print("  (a,b,c) with ac != 0 tested                  : %d" % tested)
print("  cases with NO nonzero E_r[psi_m] for m <= 24 : %d %s" % (len(bad), bad if bad else ""))
print("  largest 'first nonzero m' encountered        : %d" % worst)
print()
print("  So on the constant-coefficient {-1,0,1} stratum at M = 1, the nullcone contains")
print("  ONLY the one-sided P (ac = 0).  No asymptotics, no domination, no ESV saddle:")
print("  the answer is a Hermite polynomial, and orthogonality does the work.")
