#!/usr/bin/env python3
"""
lrc14_angleF_covolume_formula_macmini_0618s7.py  (mac-mini-2026-06-18-S7)

ANGLE F part 3 — the EXACT covolume formula and the transference structure.

CLAIM (standard geometry of numbers).  For a primitive integer row vector
a=(e_1,...,e_{k-1}) (gcd=1), the kernel lattice Lambda^o = {n in Z^{k-1}: a.n=0}
is a rank (k-2) sublattice of Z^{k-1} with
        covolume(Lambda^o) = |a|_2 = sqrt(e_1^2+...+e_{k-1}^2).
(Reason: Z^{k-1}/Lambda^o = a.Z^{k-1} = Z (primitive) embeds via a as the lattice
 |a|_2 * Z along direction a/|a|; index/covolume bookkeeping gives covol = |a|_2.)
Dually, the SHORTEST relation lambda_1(Lambda^o) is the shortest nonzero integer
solution of a.n=0 — the shortest affine relation among the offsets.

CONSEQUENCE (Angle F covolume law for the seven-sector correction).  We proved
   |corr(E)| <= sum_{0!=n in Lambda^o} |K(n)|,  |K(n)| <= sum_T prod_j |chat_T(n_j)|.
Combined with |chat_T(n)| <= |T| s7 /(pi|n|), s7=sin(pi/7), the dominant support-3
piece gives a 1/height sum.  The MINKOWSKI / counting transference:
  #{n in Lambda^o : |n|_inf <= R}  <=  C_k (R^{k-2}/covol + lower order)  (for R >= lambda_1),
so SMALL covolume => MANY short relations => large corr bound; and
  covol(Lambda^o) = |e|_2  is SCALE-equivariant: covol(dE)=d*|e_prim|_2 but the
  LATTICE Lambda^o(dE)=Lambda^o(E) is UNCHANGED (meas(S7) scale-invariance, THM-532).
So the right invariant is the PRIMITIVE covolume |e/gcd(e)|_2 — equivalently the
shape of the lattice, not its scale.  This matches meas(S7) scale-invariance.

This script:
 (A) VERIFIES covol(Lambda^o) = |e|_2 exactly (vs Gram-determinant computation);
 (B) computes lambda_1(Lambda^o) = shortest relation, and tabulates it vs corr;
 (C) tests the AP-extremal statement: AP minimises BOTH covol-per-scale and
     lambda_1 among bounded shapes (=> densest relation lattice => max corr).
"""
from __future__ import annotations
import sys, itertools, math
from fractions import Fraction as F
import sympy
sys.stdout.reconfigure(line_buffering=True)


def measS7(E):
    E = sorted(set(E))
    bps = {F(0), F(1)}
    for e in E:
        if e == 0:
            continue
        for m in range(0, 7 * e + 1):
            bps.add(F(m, 7 * e))
    bps = sorted(bps)
    total = F(0)
    for i in range(len(bps) - 1):
        x0, x1 = bps[i], bps[i + 1]
        if x1 <= x0:
            continue
        xm = (x0 + x1) / 2
        secs = {int(((e * xm) % 1) * 7) for e in E}
        if len(secs) == 7:
            total += x1 - x0
    return total


def M7(k):
    return sum(F((-1) ** t * math.comb(6, t)) * F(7 - t, 7) ** (k - 1) for t in range(7))


def lll(basis):
    if not basis:
        return []
    B = [[F(x) for x in r] for r in basis]
    n = len(B)

    def dot(u, v):
        return sum(a * b for a, b in zip(u, v))

    def gs(B):
        Bs, mu = [], [[F(0)] * n for _ in range(n)]
        for i in range(n):
            bi = list(B[i])
            for j in range(i):
                mu[i][j] = dot(B[i], Bs[j]) / dot(Bs[j], Bs[j])
                bi = [bi[t] - mu[i][j] * Bs[j][t] for t in range(len(bi))]
            Bs.append(bi)
        return Bs, mu

    Bs, mu = gs(B)
    k = 1
    while k < n:
        for j in range(k - 1, -1, -1):
            q = round(mu[k][j])
            if q:
                B[k] = [B[k][t] - q * B[j][t] for t in range(len(B[k]))]
                Bs, mu = gs(B)
        if dot(Bs[k], Bs[k]) >= (F(3, 4) - mu[k][k - 1] ** 2) * dot(Bs[k - 1], Bs[k - 1]):
            k += 1
        else:
            B[k], B[k - 1] = B[k - 1], B[k]
            Bs, mu = gs(B)
            k = max(k - 1, 1)
    return [[int(x) for x in r] for r in B]


def kernel_basis(nz):
    M = sympy.Matrix([nz])
    out = []
    for v in M.nullspace():
        L = 1
        for x in v:
            fr = F(int(x.p), int(x.q))
            L = L * fr.denominator // math.gcd(L, fr.denominator)
        iv = [int(x * L) for x in v]
        g = 0
        for c in iv:
            g = math.gcd(g, abs(c))
        if g:
            iv = [c // g for c in iv]
        out.append(iv)
    return out


def covol(basis):
    if not basis:
        return 0.0
    G = sympy.Matrix([[sum(a * b for a, b in zip(u, w)) for w in basis] for u in basis])
    return math.sqrt(float(G.det()))


def shortest_relation(nz, mc=3):
    red = lll(kernel_basis(nz))
    d = len(red)
    best = None
    for c in itertools.product(range(-mc, mc + 1), repeat=d):
        if all(x == 0 for x in c):
            continue
        v = [sum(c[i] * red[i][t] for i in range(d)) for t in range(len(nz))]
        nv = sum(x * x for x in v)
        if best is None or nv < best[0]:
            best = (nv, v)
    return math.sqrt(best[0]), best[1]


def main():
    print("=" * 90)
    print("ANGLE F part 3: exact covolume formula covol(Lambda^o)=|e|_2 + transference")
    print("=" * 90)
    print("(A) Verify covol(Lambda^o) = |e_prim|_2 :")
    print(f"{'shape':<34}{'covol(Gram)':>14}{'|e_prim|_2':>14}{'match':>8}")
    test = [
        [1, 2, 3, 4],
        [1, 2, 3, 4, 5, 6, 7],
        [2, 3, 4, 5, 6, 7, 9],
        [1, 3, 7, 12, 20, 30, 44],
        [1, 3, 7, 15, 31, 63, 127],
        [5, 13, 27, 41, 58, 79, 97],
    ]
    for nz in test:
        g = 0
        for x in nz:
            g = math.gcd(g, x)
        prim = [x // g for x in nz]
        cv = covol(kernel_basis(nz))
        norm = math.sqrt(sum(x * x for x in prim))
        ok = "YES" if abs(cv - norm) < 1e-6 else "NO"
        print(f"{str(nz):<34}{cv:>14.4f}{norm:>14.4f}{ok:>8}")
    print()
    print("(B) shortest relation lambda_1 vs corr(E), k=8:")
    print(f"{'shape':<34}{'corr':>10}{'covol=|e|2':>12}{'lambda1':>10}{'shortest rel vec':>22}")
    m7 = float(M7(8))
    shapes = [
        ("consec {0..7}", list(range(8))),
        ("3AP {..,6,8}", [0, 1, 2, 3, 4, 5, 6, 8]),
        ("2run40", [0, 1, 2, 3, 40, 41, 42, 43]),
        ("perf", [0, 2, 3, 4, 5, 6, 7, 9]),
        ("Sidon", [0, 1, 3, 7, 12, 20, 30, 44]),
        ("dissoc", [0, 1, 3, 7, 15, 31, 63, 127]),
    ]
    for name, E in shapes:
        nz = [e for e in E if e != 0]
        corr = float(measS7(E)) - m7
        cv = covol(kernel_basis(nz))
        l1, vec = shortest_relation(nz, 3)
        print(f"{name:<34}{corr:>10.5f}{cv:>12.2f}{l1:>10.3f}{str(vec):>22}")
    print()
    print("(C) AP minimises covol-per-scale & lambda_1 among bounded shapes (exhaustive):")
    for kk in (5, 6, 7, 8):
        cap = kk + 4
        best_cv = best_E = None
        best_l1 = best_l1E = None
        ap_cv = ap_l1 = None
        for combo in itertools.combinations(range(1, cap + 1), kk - 1):
            E = (0,) + combo
            g = 0
            for x in E:
                g = math.gcd(g, x)
            if g != 1:
                continue
            nz = list(combo)
            cv = math.sqrt(sum(x * x for x in nz))  # = |e|_2 (gcd=1)
            l1, _ = shortest_relation(nz, 2)
            if best_cv is None or cv < best_cv:
                best_cv, best_E = cv, E
            if best_l1 is None or l1 < best_l1:
                best_l1, best_l1E = l1, E
            if E == tuple(range(kk)):
                ap_cv, ap_l1 = cv, l1
        print(
            f"k={kk}: min covol={best_cv:.3f}@{best_E} (AP {ap_cv:.3f}) | "
            f"min lambda1={best_l1:.3f}@{best_l1E} (AP {ap_l1:.3f}) "
            f"{'AP=both-min' if best_E==tuple(range(kk)) and best_l1E==tuple(range(kk)) else ''}"
        )
    print()
    print("DONE.")


if __name__ == "__main__":
    main()
