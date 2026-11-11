#!/usr/bin/env python3
"""S91 rates_exact (A.1, A.3): exact arithmetic of the twisted transfer matrices T_k of twisted_spectrum.py.

Everything in sections 1-5 is exact arithmetic in Z[zeta_m] = Z[z] / Phi_m(z) (integer coordinates on the power
basis) or in Z[t] (python-flint).  Only the identification of the root of largest modulus (section 6) uses floating
point (mpmath, 60 digits).

Conventions (as in twisted_spectrum.py): states (eta, beta), eta in Z/m, beta in {0, 1}, index 2 eta + beta;
T_k[(eta, beta) -> (2 eta, 0)] = 1, T_k[(eta, beta) -> (2 eta + 1, 1)] = w^(eta + beta), w = zeta^(k/2).

Checked here:
 1. (complement symmetry)  C_k T_k = T_k C_k and C_k^2 = zeta^k, (C u)(eta, beta) = w^(-(eta - beta)) u(-1 - eta, 1 - beta),
    for every odd m <= MSYM and every k; consequently char(T_k) = char(T_k^+) char(T_k^-) with the m x m matrices
    (T_k^{+-} f)(eta) = f(2 eta) +- w^(-(eta + 1)) f(-2 eta - 2)   (product checked exactly for m <= MCHAR, all k).
 2. (norm) T_k T_k^* is the direct sum of the blocks [[2, 1 + w^-1], [1 + w, 2]]: singular values sqrt(2 +- |1 + w|).
    Hence rho(T_k) <= sqrt(2 + |1 + w|) and |det T_k| = |1 - w|^m.
 3. m = 3: u = (1, -omega, omega, -omega^2, omega, -omega^2) is an eigenvector of T_1 for 2 + omega, |2 + omega|^2 = 3
    = 2 + |1 + w|: the spectral radius is exactly sqrt(3) (the operator norm).  Characteristic polynomials printed.
 4. m = 3, 5, 7: the norm polynomials N^{+-}(t) = prod over Galois conjugates of char(T_1^{+-}) in Z[t], their
    factorisation over Q, cross-checked against the characteristic polynomial of the integer companion-block matrix.
 5. m = 5, 7: char(T_k^{+-}) are irreducible of degree m over Q(zeta_m) (their norms are irreducible over Q), and
    their Galois groups over Q(zeta_m) contain a transposition (Dedekind: a prime p = 1 mod m where the reduction has
    exactly one quadratic factor and otherwise factors of odd degree): the groups are S_m, not solvable.
 6. the eigenvalue of largest modulus and its minimal polynomial over Q.

Usage: python3 exact_charpoly.py            (about half a minute)
"""
import sys
from math import gcd

import mpmath as mp
from flint import fmpz_mat, fmpz_poly, nmod_poly

MSYM = 45
MCHAR = 13


# ----------------------------------------------------------------------------------------------- the ring Z[zeta_m]
class Ring:
    def __init__(self, m):
        self.m = m
        x = fmpz_poly([0, 1])
        phi = {}
        for d in range(1, m + 1):
            if m % d == 0:
                p = x ** d - 1
                for e in range(1, d):
                    if d % e == 0:
                        p = p // phi[e]
                phi[d] = p
        self.phi = [int(c) for c in phi[m].coeffs()]            # ascending, monic
        self.d = len(self.phi) - 1
        d = self.d
        # z^j mod Phi for j < max(2d - 1, m)
        self.red = []
        cur = [1] + [0] * (d - 1)
        for _ in range(max(2 * d - 1, m) + 1):
            self.red.append(tuple(cur))
            top = cur[d - 1]
            cur = [0] + cur[:d - 1]
            if top:
                cur = [c - top * self.phi[i] for i, c in enumerate(cur)]
        self.zero = tuple([0] * d)
        self.one = self.red[0]

    def zpow(self, e):
        return self.red[e % self.m]

    def add(self, a, b):
        return tuple(x + y for x, y in zip(a, b))

    def sub(self, a, b):
        return tuple(x - y for x, y in zip(a, b))

    def neg(self, a):
        return tuple(-x for x in a)

    def scal(self, c, a):
        return tuple(c * x for x in a)

    def mul(self, a, b):
        d = self.d
        conv = [0] * (2 * d - 1)
        for i, ai in enumerate(a):
            if ai:
                for j, bj in enumerate(b):
                    if bj:
                        conv[i + j] += ai * bj
        res = [0] * d
        for j, c in enumerate(conv):
            if c:
                r = self.red[j]
                for i in range(d):
                    res[i] += c * r[i]
        return tuple(res)

    def sigma(self, a, g):
        """Galois conjugation zeta -> zeta^g."""
        res = [0] * self.d
        for i, c in enumerate(a):
            if c:
                r = self.zpow(g * i)
                for j in range(self.d):
                    res[j] += c * r[j]
        return tuple(res)

    def divexact(self, a, k):
        if any(x % k for x in a):
            raise ArithmeticError("inexact division in Z[zeta]")
        return tuple(x // k for x in a)

    def show(self, a):
        terms = []
        for i, c in enumerate(a):
            if c:
                terms.append((f"{c}" if i == 0 else (("" if abs(c) == 1 else str(abs(c)) + "*") if c > 0 else ("-" if c == -1 else str(c) + "*")) + ("z" if i == 1 else f"z^{i}")))
        return "(" + " + ".join(terms).replace("+ -", "- ") + ")" if terms else "0"


# ------------------------------------------------------------------------------------------- matrices over the ring
def sparse_T(R, k):
    """rows of T_k as lists of (column, ring element)."""
    m = R.m
    inv2 = pow(2, -1, m)
    rows = []
    for eta in range(m):
        for beta in range(2):
            e = (k * inv2 * (eta + beta)) % m
            rows.append([(2 * ((2 * eta) % m), R.one), (2 * ((2 * eta + 1) % m) + 1, R.zpow(e))])
    return rows


def sparse_reduced(R, k, sign):
    m = R.m
    inv2 = pow(2, -1, m)
    rows = []
    for eta in range(m):
        ph = R.zpow(-(k * inv2 * (eta + 1)))
        rows.append([((2 * eta) % m, R.one), ((-2 * eta - 2) % m, ph if sign > 0 else R.neg(ph))])
    return rows


def sparse_C(R, k):
    m = R.m
    inv2 = pow(2, -1, m)
    rows = []
    for eta in range(m):
        for beta in range(2):
            rows.append([(2 * ((-1 - eta) % m) + (1 - beta), R.zpow(-(k * inv2 * (eta - beta))))])
    return rows


def dense(R, rows, n):
    M = [[R.zero] * n for _ in range(n)]
    for i, row in enumerate(rows):
        for c, a in row:
            M[i][c] = R.add(M[i][c], a)
    return M


def sp_times_dense(R, rows, M):
    n = len(M)
    out = []
    for row in rows:
        acc = [R.zero] * n
        for c, a in row:
            Mc = M[c]
            for j in range(n):
                if any(Mc[j]):
                    acc[j] = R.add(acc[j], R.mul(a, Mc[j]))
        out.append(acc)
    return out


def charpoly(R, rows):
    """Faddeev-LeVerrier; returns the coefficients c[0..n] (ring elements), c[n] = 1."""
    n = len(rows)
    M = [[R.one if i == j else R.zero for j in range(n)] for i in range(n)]
    c = [None] * (n + 1)
    c[n] = R.one
    for k in range(1, n + 1):
        AM = sp_times_dense(R, rows, M)
        tr = R.zero
        for i in range(n):
            tr = R.add(tr, AM[i][i])
        ck = R.divexact(R.neg(tr), k)
        c[n - k] = ck
        for i in range(n):
            AM[i][i] = R.add(AM[i][i], ck)
        M = AM
    return c


def pmul(R, a, b):
    out = [R.zero] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if any(x):
            for j, y in enumerate(b):
                if any(y):
                    out[i + j] = R.add(out[i + j], R.mul(x, y))
    return out


def peval(R, c, x):
    acc = R.zero
    for a in reversed(c):
        acc = R.add(R.mul(acc, x), a)
    return acc


def show_poly(R, c):
    return " + ".join(f"{R.show(a)} t^{i}" for i, a in reversed(list(enumerate(c))) if any(a))


def companion_block_charpoly(m, k):
    """characteristic polynomial in Z[t] of T_k with zeta replaced by the companion matrix of Phi_m."""
    R = Ring(m)
    d = R.d
    n = 2 * m
    B = [[0] * (n * d) for _ in range(n * d)]
    for i, row in enumerate(sparse_T(R, k)):
        for col, a in row:
            # multiplication by a on the basis 1, z, ..., z^(d-1): column b of the block is a * z^b
            for b in range(d):
                prod = R.mul(a, R.red[b])
                for r in range(d):
                    B[i * d + r][col * d + b] += prod[r]
    return fmpz_mat(B).charpoly()


# ---------------------------------------------------------------------------------------------------------- checks
def main():
    ok = True
    mp.mp.dps = 60

    # 1. complement symmetry, exact, all odd m <= MSYM, all k
    for m in range(3, MSYM + 1, 2):
        R = Ring(m)
        n = 2 * m
        for k in range(1, m):
            Trows, Crows = sparse_T(R, k), sparse_C(R, k)
            Td, Cd = dense(R, Trows, n), dense(R, Crows, n)
            ok &= sp_times_dense(R, Trows, Cd) == sp_times_dense(R, Crows, Td)
            C2 = sp_times_dense(R, Crows, Cd)
            zk = R.zpow(k)
            ok &= all(C2[i][j] == (zk if i == j else R.zero) for i in range(n) for j in range(n))
            # 2. Gram matrix T T^*
            w = R.zpow(k * pow(2, -1, m))
            winv = R.zpow(-k * pow(2, -1, m))
            for eta in range(m):
                r0 = dict(Trows[2 * eta])
                r1 = dict(Trows[2 * eta + 1])
                ok &= set(r0) == set(r1)
                g01 = R.zero
                for col in r0:
                    g01 = R.add(g01, R.mul(r0[col], R.sigma(r1[col], -1)))
                ok &= g01 == R.add(R.one, winv)
                ok &= R.sigma(g01, -1) == R.add(R.one, w)
            cols = [c for row in Trows[0::2] for c, _ in row]
            ok &= len(set(cols)) == n            # rows with different eta have disjoint supports
    print(f"1-2. C T = T C, C^2 = zeta^k, and T T^* = sum of blocks [[2, 1 + 1/w], [1 + w, 2]]: exact for all odd m <= {MSYM}, all k: {ok}")

    chars = {}
    for m in range(3, MCHAR + 1, 2):
        R = Ring(m)
        good = True
        for k in range(1, m):
            full = charpoly(R, sparse_T(R, k))
            cp = charpoly(R, sparse_reduced(R, k, +1))
            cm = charpoly(R, sparse_reduced(R, k, -1))
            good &= pmul(R, cp, cm) == full
            chars[(m, k)] = (cp, cm)
            if k > 1:           # Galois covariance sigma_k(char T_1) = char T_k (only meaningful when gcd(k, m) = 1)
                if gcd(k, m) == 1:
                    good &= [R.sigma(a, k) for a in chars[(m, 1)][0]] == cp and [R.sigma(a, k) for a in chars[(m, 1)][1]] == cm
        ok &= good
        print(f"   m = {m:2d}: char(T_k) = char(T_k^+) char(T_k^-) and sigma_k(char T_1^{{+-}}) = char T_k^{{+-}} for all k: {good}")

    # 3. m = 3
    R = Ring(3)
    om = R.zpow(1)
    lam = R.add(R.scal(2, R.one), om)                           # 2 + omega
    u = [R.one, R.neg(om), om, R.neg(R.zpow(2)), om, R.neg(R.zpow(2))]
    Tu = []
    for row in sparse_T(R, 1):
        acc = R.zero
        for c, a in row:
            acc = R.add(acc, R.mul(a, u[c]))
        Tu.append(acc)
    eig = Tu == [R.mul(lam, x) for x in u]
    nrm = R.mul(lam, R.sigma(lam, -1))
    ok &= eig and nrm == R.scal(3, R.one)
    w = R.zpow(pow(2, -1, 3))
    onepw = R.add(R.one, w)
    ok &= R.mul(onepw, R.sigma(onepw, -1)) == R.one             # |1 + w| = 1, so ||T_1||^2 = 2 + 1 = 3
    print(f"3. m = 3: T_1 u = (2 + omega) u for u = (1, -omega, omega, -omega^2, omega, -omega^2): {eig}; "
          f"|2 + omega|^2 = 3 and |1 + w|^2 = 1, so rho(T_1) = ||T_1|| = sqrt(3): {nrm == R.scal(3, R.one)}")
    cp, cm = chars[(3, 1)]
    print("   char(T_1^+) =", show_poly(R, cp))
    print("   char(T_1^-) =", show_poly(R, cm))
    rootp, rootm = any(peval(R, cp, lam)) is False, any(peval(R, cm, lam)) is False
    ok &= rootp != rootm
    print(f"   2 + omega is a root of char(T_1^+): {rootp}; of char(T_1^-): {rootm}")

    # 4.-6. norms, factorisation, Galois groups, dominant eigenvalue
    results = {}
    for m in (3, 5, 7):
        R = Ring(m)
        units = [g for g in range(1, m) if gcd(g, m) == 1]
        norms = {}
        for name, idx in (("+", 0), ("-", 1)):
            prod = [R.one]
            for g in units:
                prod = pmul(R, prod, [R.sigma(a, g) for a in chars[(m, 1)][idx]])
            good = all(not any(a[1:]) for a in prod)
            ok &= good
            norms[name] = fmpz_poly([a[0] for a in prod])
        full = companion_block_charpoly(m, 1)
        same = (norms["+"] * norms["-"]) == full
        ok &= same
        print(f"4. m = {m}: N^+ N^- equals the characteristic polynomial of the {2 * m * R.d} x {2 * m * R.d} integer block matrix: {same}")
        best = None
        for name in ("+", "-"):
            _, facs = norms[name].factor()
            print(f"   N^{name}: degrees of the irreducible factors over Q: {sorted(f.degree() for f, _ in facs)}")
            for f, mult in facs:
                co = [int(c) for c in f.coeffs()]
                roots = mp.polyroots(co[::-1], maxsteps=5000, extraprec=600)
                mods = sorted((abs(r) for r in roots), reverse=True)
                if m == 3:
                    print(f"      {f}   moduli of the roots: {[mp.nstr(v, 12) for v in mods]}")
                if best is None or mods[0] > best[0]:
                    second = best[0] if best else mp.mpf(0)
                    best = (mods[0], name, f, mods)
        top, name, f, mods = best
        # the second largest modulus over all roots of N^+ N^- that is not equal to the top (complex pairs share it)
        allmods = []
        for nm in ("+", "-"):
            _, facs = norms[nm].factor()
            for g, _ in facs:
                co = [int(c) for c in g.coeffs()]
                allmods += [abs(r) for r in mp.polyroots(co[::-1], maxsteps=5000, extraprec=600)]
        allmods.sort(reverse=True)
        ntop = sum(1 for v in allmods if abs(v - top) < mp.mpf(10) ** (-40))
        nxt = max(v for v in allmods if abs(v - top) >= mp.mpf(10) ** (-40))
        results[m] = (top, name, f)
        print(f"6. m = {m}: largest modulus of an eigenvalue of T_k, k != 0: {mp.nstr(top, 40)}  (rate {mp.nstr(top / 2, 20)}); "
              f"attained by {ntop} roots of N^{name} (one per Galois-conjugate block pair); next modulus {mp.nstr(nxt, 12)}")
        print(f"   minimal polynomial over Q of the dominant eigenvalue (degree {f.degree()}):")
        print("     ", f)
        if m > 3:
            irred = all(len(norms[nm].factor()[1]) == 1 and norms[nm].factor()[1][0][1] == 1 for nm in ("+", "-"))
            ok &= irred
            print(f"5. m = {m}: N^+ and N^- are irreducible over Q of degree {m * R.d}, so char(T_k^+-) are irreducible over Q(zeta_{m}): {irred}")
            for nm, idx in (("+", 0), ("-", 1)):
                c = chars[(m, 1)][idx]
                found = None
                types = {}
                p = 1
                while p < 4000:
                    p += m
                    if p % 2 == 0 or any(p % s == 0 for s in range(3, int(p ** 0.5) + 1, 2)):
                        continue
                    # a primitive m-th root of unity modulo p
                    g = next(pow(a, (p - 1) // m, p) for a in range(2, p) if pow(a, (p - 1) // m, p) != 1)
                    co = [sum(ci * pow(g, i, p) for i, ci in enumerate(a)) % p for a in c]
                    _, fac = nmod_poly(co, p).factor()
                    if any(e > 1 for _, e in fac):
                        continue
                    t = tuple(sorted(fp.degree() for fp, _ in fac))
                    types[t] = types.get(t, 0) + 1
                    if found is None and t.count(2) == 1 and all(d % 2 == 1 for d in t if d != 2):
                        found = (p, t)
                ok &= found is not None
                print(f"   char(T_1^{nm}) over Q(zeta_{m}): Frobenius cycle types at the split primes p = 1 mod {m}, p < 4000: "
                      f"{dict(sorted(types.items()))}")
                print(f"      transposition certificate: p = {found[0]}, factor degrees {found[1]} -> the Galois group is S_{m} (not solvable)")
            R7 = R
            print(f"   char(T_1^{name}) (the block of N^{name}; coefficients in Z[z]/Phi_{m}, z = zeta_{m}):")
            print("     ", show_poly(R7, chars[(m, 1)][0 if name == "+" else 1]))
    # exact values at m = 3
    ok &= str(results[3][2]) == "x^2 + (-3)*x + 3"
    ok &= abs(results[3][0] - mp.sqrt(3)) < mp.mpf(10) ** (-50)
    ok &= abs(results[5][0] / 2 - mp.mpf("0.741760072680922990")) < mp.mpf(10) ** (-15)
    ok &= abs(results[7][0] / 2 - mp.mpf("0.932038815554839891")) < mp.mpf(10) ** (-15)
    print("ALL CHECKS PASSED" if ok else "CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
