#!/usr/bin/env python3
"""Audit proof-5128-AB: three extra checks.

1. (A) is attained: for r = 2, B the dot product, L != L' two lines through 0 in (Z/p)^2 and small t,
       D1 = U + t (1_{L minus 0} - 1_{L' minus 0}),   D2 = U + t (1_{Lperp minus 0} - 1_{L'perp minus 0})
   give Rhat(0) = ((p-1)/p) p^(r/2) |D1 - U|_2 |D2 - U|_2 exactly (rationals): the constant in (2) cannot be lowered
   for c = 0, r = 2.  For c != 0 the exact value of Rhat(c) is printed (it is smaller in modulus).
2. (B1) at m = 3, exactly: characteristic polynomial of T_1 over Q(zeta_3) (sympy) and its largest root modulus.
3. (B2): the cell formula.  For random cells (top word, bottom word, r free digits) the character sum over the cell
   equals (e_(0,0)^T T[top]) T_k^r (T[bottom] v_l), and its modulus is at most the l^infinity operator norm of T_k^r.
"""
import itertools
import sys
from fractions import Fraction

import numpy as np
import sympy as sp


def part1():
    ok = True
    for p in (3, 5, 7, 11):
        V = list(itertools.product(range(p), repeat=2))
        n = len(V)
        U = Fraction(1, n)
        t = Fraction(1, 10 * n)
        L = {(a % p, 0) for a in range(1, p)}                      # line spanned by (1, 0)
        L2 = {(a % p, a % p) for a in range(1, p)}                 # line spanned by (1, 1)
        Lp = {(0, a % p) for a in range(1, p)}                     # (1,0)^perp
        L2p = {(a % p, (-a) % p) for a in range(1, p)}             # (1,1)^perp
        D1 = {v: U + t * ((v in L) - (v in L2)) for v in V}
        D2 = {v: U + t * ((v in Lp) - (v in L2p)) for v in V}
        assert sum(D1.values()) == 1 and sum(D2.values()) == 1 and min(D1.values()) > 0 and min(D2.values()) > 0
        P = [Fraction(0)] * p
        for u in V:
            for v in V:
                P[(u[0] * v[0] + u[1] * v[1]) % p] += D1[u] * D2[v]
        pi1, pi2 = D1[(0, 0)], D2[(0, 0)]
        Rhat = [P[c] - Fraction(1, p) - (pi1 + pi2 - U) * ((c == 0) - Fraction(1, p)) for c in range(p)]
        n1 = sum((D1[v] - U) ** 2 for v in V)
        n2 = sum((D2[v] - U) ** 2 for v in V)
        bound_sq = Fraction(p - 1, p) ** 2 * p ** 2 * n1 * n2
        eq = Rhat[0] ** 2 == bound_sq
        ok &= eq and all(Rhat[c] ** 2 <= bound_sq for c in range(p))
        others = sorted({str(Rhat[c] / Rhat[0]) for c in range(1, p)})
        print(f"   p = {p:2d}: Rhat(0)^2 == bound^2 exactly: {eq};  Rhat(c)/Rhat(0) for c != 0: {others}")
    return ok


def part2():
    w = sp.Rational(-1, 2) + sp.sqrt(3) * sp.I / 2          # zeta_3
    m, k = 3, 1
    inv2 = 2
    T = sp.zeros(6, 6)
    for eta in range(m):
        for beta in range(2):
            src = 2 * eta + beta
            T[src, 2 * ((2 * eta) % m)] = 1
            T[src, 2 * ((2 * eta + 1) % m) + 1] = w ** ((k * ((eta + beta) * inv2 % m)) % m)
    lam = sp.symbols("lam")
    cpoly = sp.expand(T.charpoly(lam).as_expr())
    cp = sp.factor(cpoly, extension=[sp.sqrt(3), sp.I])
    print("   m = 3, k = 1: det(lam - T_1) =", cp)
    top = sp.Rational(3, 2) + sp.sqrt(3) * sp.I / 2
    is_root = sp.simplify(cpoly.subs(lam, top)) == 0
    print(f"   lam = 3/2 + sqrt(3) i/2 is an eigenvalue: {is_root}; its modulus squared = {sp.simplify(sp.Abs(top) ** 2)}")
    roots = np.roots([complex(c) for c in sp.Poly(cpoly, lam).all_coeffs()])
    mods = sorted((abs(complex(r)) for r in roots), reverse=True)
    print("   root moduli:", ", ".join(f"{v:.9f}" for v in mods), "  sqrt(3) =", f"{3 ** 0.5:.9f}")
    return bool(is_root) and abs(mods[0] - 3 ** 0.5) < 1e-6


def part3():
    rng = np.random.default_rng(5128)
    ok = True
    nb = 16
    size = 1 << nb
    q = np.zeros(size, dtype=np.int64)
    for kk in range(nb):
        lo, hi = 1 << kk, 1 << (kk + 1)
        idx = np.arange(lo, hi, dtype=np.int64)
        half = idx >> 1
        q[lo:hi] = q[half] + (idx & 1) * ((half + 1) >> 1)
    worst = 0.0
    worst_ratio = 0.0
    for m in (3, 5, 7, 9, 15):
        inv2 = (m + 1) // 2
        for trial in range(40):
            k = int(rng.integers(0, m))
            l = int(rng.integers(0, m))
            T0 = np.zeros((2 * m, 2 * m), dtype=np.complex128)
            T1 = np.zeros((2 * m, 2 * m), dtype=np.complex128)
            for eta in range(m):
                for beta in range(2):
                    src = 2 * eta + beta
                    T0[src, 2 * ((2 * eta) % m)] = 1.0
                    T1[src, 2 * ((2 * eta + 1) % m) + 1] = np.exp(2j * np.pi * ((k * ((eta + beta) * inv2 % m)) % m) / m)
            Tk = T0 + T1
            n = int(rng.integers(4, nb + 1))
            t = int(rng.integers(1, n))
            j = int(rng.integers(0, n - t + 1))
            r = n - t - j
            top = (1 << (t - 1)) | int(rng.integers(0, 1 << (t - 1)))
            bottom = int(rng.integers(0, 1 << j))
            free = np.arange(1 << r, dtype=np.int64)
            hs = (top << (n - t)) | (free << j) | bottom
            direct = np.sum(np.exp(2j * np.pi * ((k * q[hs] + l * hs) % m) / m))
            row = np.zeros(2 * m, dtype=np.complex128)
            row[0] = 1.0
            for i in range(t - 1, -1, -1):
                row = row @ (T1 if (top >> i) & 1 else T0)
            col = np.array([np.exp(2j * np.pi * ((l * eta) % m) / m) for eta in range(m) for _ in range(2)])
            for i in range(j):
                col = (T1 if (bottom >> i) & 1 else T0) @ col
            Tr = np.linalg.matrix_power(Tk, r)
            val = row @ Tr @ col
            worst = max(worst, abs(val - direct))
            opnorm = float(np.max(np.sum(np.abs(Tr), axis=1)))
            ok &= abs(direct) <= opnorm + 1e-6
            ok &= abs(np.max(np.abs(row)) - 1) < 1e-9 and int(np.sum(np.abs(row) > 1e-9)) == 1
            ok &= np.max(np.abs(col)) <= 1 + 1e-9
            if k:
                worst_ratio = max(worst_ratio, abs(direct) / opnorm)
    print(f"   cell formula: max |direct - (e^T T[top]) T_k^r (T[bottom] v_l)| = {worst:.2e}; "
          f"max |cell sum| / ||T_k^r||_inf = {worst_ratio:.3f}")
    return ok and worst < 1e-6


def main():
    print("1. the constant of (2) is attained (c = 0, r = 2)")
    ok1 = part1()
    print("2. exact characteristic polynomial at m = 3")
    ok2 = part2()
    print("3. cell character sums")
    ok3 = part3()
    print("ALL CHECKS PASSED" if (ok1 and ok2 and ok3) else "CHECK FAILED")
    sys.exit(0 if (ok1 and ok2 and ok3) else 1)


if __name__ == "__main__":
    main()
