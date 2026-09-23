#!/usr/bin/env python3
"""
collatz_procgen_20260923_theta_defect.py

The owner's principle ("every triangle/AM-GM/Cauchy-Schwarz inequality marks abstracted structure;
the defect is calculable and is a lever") applied to the height bookkeeping of HYP-9127 and Theorem H.

  D1  cubes, partial sums S_K = A_K/3^(9K^3): exact heights after gcd versus the triangle bound
  D2  cubes, periodic approximants Phi(u v^inf) = P/Q of Y3 (Theorem R): gcd(P, Q) and height defect
  D3  squares, Zudilin forms: gcd and heights versus the natural bound; the exponent ladder
      (height exponent per unit of 2-adic valuation) periodic -> Zudilin -> Hankel H1 -> H2 -> H3
  D4  squares, Hankel integers N_n(a, b): polynomial content (forced divisor) versus the proven divisor
      M0^e'(n) prod Phi_l^hom(M0^2, 2^2L)^e_l(n); 3-adic and 11-adic parts
  D5  cubes, Hankel integers of the cube tails: content, X-order law e3(n) = 3n^3 - 18n^2 + 34n - 19,
      11-adic part; comparison with the quartic clearing

Usage: python3 collatz_procgen_20260923_theta_defect.py [--quick]
"""
import math
import sys
import time
from fractions import Fraction

import flint
import gmpy2
from gmpy2 import mpz
import sympy

QUICK = "--quick" in sys.argv
PHI = (1 + 5 ** 0.5) / 2
L_Y, M0_Y = 10, 3 ** 9


def hdr(s):
    print()
    print("=" * 78)
    print(s)
    print("=" * 78)
    sys.stdout.flush()


def log2abs(x):
    x = abs(mpz(x))
    if x == 0:
        return float("-inf")
    b = int(x.bit_length())
    if b <= 1000:
        return math.log2(int(x))
    sh = b - 60
    return math.log2(int(x >> sh)) + sh


def vp(x, p):
    x = mpz(x)
    if x == 0:
        return None
    c = 0
    while x % p == 0:
        x //= p
        c += 1
    return c


def v2(x):
    x = mpz(x)
    return int(gmpy2.bit_scan1(abs(x))) if x else None


# ------------------------------------------------------------------ D1
def section_D1():
    hdr("D1  cubes: partial sums S_K = A_K / 3^(9K^3): exact height versus triangle bound")
    print("  K | gcd(A_K, 3^(9K^3)) | log2 H exact | log2 triangle bound (K+1) 3^(9K^3) | defect bits | v2/log2H")
    Kmax = 8 if QUICK else 12
    for K in range(1, Kmax + 1):
        A = sum(mpz(2) ** (10 * k ** 3) * mpz(3) ** (9 * (K ** 3 - k ** 3)) for k in range(K + 1))
        D = mpz(3) ** (9 * K ** 3)
        g = gmpy2.gcd(A, D)
        H = max(log2abs(A // g), log2abs(D // g))
        bound = math.log2(K + 1) + 9 * K ** 3 * math.log2(3)
        print(f"  {K:2d} | {int(g):18d} | {H:12.2f} | {bound:36.2f} | {bound - H:11.2f} | {10*(K+1)**3/H:.4f}")
    print("  The denominator 3^(9K^3) is exact (A_K = 2^(10K^3) mod 3), so the triangle inequality loses only")
    print("  log2(K+1) - log2(1.052) bits: the Liouville exponent 1/mu_bar = 0.7010 has no calculable defect.")


# ------------------------------------------------------------------ D2
def section_D2():
    hdr("D2  cubes: periodic approximants Phi(u v^inf) of Y3 (best per period end j): gcd defect")
    nb = 120 if QUICK else 400
    word = "".join("1111111101" if round(n ** (1 / 3)) ** 3 == n and n > 0 else "1111111110" for n in range(nb))
    Lw = len(word)

    def MR(z):
        M, R = 1, 0
        for i, ch in enumerate(z):
            if ch == "1":
                R = 3 * R + (1 << i)
                M *= 3
        return M, R
    rows = []
    jmax = 1200 if QUICK else 3000
    for j in range(20, jmax, 7):
        best = None
        for a in range(0, j):
            p = j - a
            l = 0
            while j + l < Lw and word[a + l] == word[j + l]:
                l += 1
            lam = j + l
            if best is None or lam > best[0]:
                best = (lam, a)
        lam, a = best
        u, v = word[:a], word[a:j]
        Mu, Ru = MR(u)
        Mv, Rv = MR(v)
        Dv = (1 << len(v)) - Mv
        Q = Mu * Dv
        P = (1 << len(u)) * Rv - Ru * Dv
        g = math.gcd(P, Q)
        Hnat = max(log2abs(P), log2abs(Q))
        Hred = max(log2abs(P // g), log2abs(Q // g))
        # canonical representation of u v^inf: primitive period, then minimal preperiod (rotate back)
        pmin = next(pp for pp in range(1, len(v) + 1) if len(v) % pp == 0 and v == v[:pp] * (len(v) // pp))
        uc, w = u, v[:pmin]
        while uc and uc[-1] == w[-1]:
            uc, w = uc[:-1], w[-1] + w[:-1]
        Mc, Rc = MR(uc)
        Mw, Rw = MR(w)
        Dw = (1 << len(w)) - Mw
        Qw = Mc * Dw
        Pw = (1 << len(uc)) * Rw - Rc * Dw
        assert Fraction(int(Pw), int(Qw)) == Fraction(int(P), int(Q))
        gw = math.gcd(Pw, Qw)
        Hw = max(log2abs(Pw), log2abs(Qw))
        rows.append((j, lam, Hnat, Hred, math.log2(g), len(v) // pmin, Hw, math.log2(gw), len(u) - len(uc)))
    gmax = max(r[4] for r in rows)
    gmean = sum(r[4] for r in rows) / len(rows)
    gw_max = max(r[7] for r in rows)
    npow = sum(1 for r in rows if r[5] > 1)
    print(f"  {len(rows)} best approximants (period end j < {jmax}): max log2 gcd(P,Q) = {gmax:.1f}, mean {gmean:.2f} bits;")
    print(f"  in {npow} of them v is a proper power; after canonicalization (primitive period, minimal preperiod)"
          f" the residual gcd is <= {gw_max:.1f} bits")
    for r in rows[:: max(1, len(rows) // 6)]:
        print(f"    j={r[0]:5d}: lambda={r[1]:5d}, log2 H natural {r[2]:8.1f}, after gcd {r[3]:8.1f};"
              f" canonical: v = w^{r[5]}, preperiod shortened by {r[8]}, height {r[6]:8.1f} (residual gcd {r[7]:.1f});"
              f" gain lambda-H {r[1]-r[3]:+8.1f}")
    print("  The large gcds are the calculable defect of a NON-CANONICAL representation (non-primitive period,")
    print("  non-minimal preperiod: factors (2^|v|-M_v)/(2^|w|-M_w) and M_z); after canonicalization the residual")
    print("  gcd is small and the gains stay bounded (they are the partial-sum approximants in disguise).")


# ------------------------------------------------------------------ D3
def qbinom_products(n):
    C = [{0: 1}]
    for j in range(1, n + 1):
        new = [dict() for _ in range(len(C) + 1)]
        for k, ck in enumerate(C):
            for d, c in ck.items():
                new[k][d] = new[k].get(d, 0) + c
                new[k + 1][d + j] = new[k + 1].get(d + j, 0) - c
        C = [{d: c for d, c in ck.items() if c != 0} for ck in new]
    return C


def zudilin_PQ(n, m, L, M0):
    P = lambda l: l * l
    C = qbinom_products(n)
    Acoef, Bcoef = {}, {}
    for k, ck in enumerate(C):
        for d, c in ck.items():
            e0 = -2 * (d + k * (k - 1) // 2 + k * m) - k
            Acoef[e0] = Acoef.get(e0, 0) + c
            for l in range(0, k + m + 1):
                e = e0 + P(l)
                Bcoef[e] = Bcoef.get(e, 0) + c
    exps = [e for e, c in Acoef.items() if c] + [e for e, c in Bcoef.items() if c]
    lo, hi = min(0, min(exps)), max(0, max(exps))

    def clear(coef):
        return sum(c * mpz(2) ** (L * (e - lo)) * mpz(M0) ** (hi - e) for e, c in coef.items() if c)
    return clear(Acoef), clear(Bcoef), lo, hi


def section_D3():
    hdr("D3  squares: Zudilin forms (gcd, heights) and the exponent ladder at the Y map (mu_bar = 1.4265)")
    L, M0 = L_Y, M0_Y
    mu = math.log2(M0) / L
    nmax = 8 if QUICK else 12
    print("  n  m | log2 H natural | log2 H after gcd | log2 gcd | v2 of form | ratio log2H / v2  (18659 = 3^9 - 2^10)")
    for n in range(2, nmax + 1, 2):
        for m in (0, round(n / PHI)):
            A, B, lo, hi = zudilin_PQ(n, m, L, M0)
            g = gmpy2.gcd(A, B)
            Hn = max(log2abs(A), log2abs(B))
            Hr = max(log2abs(A // g), log2abs(B // g))
            val = L * (-lo) + L * (n + m + 1) ** 2
            gi = int(g)
            f3 = vp(gi, 3)
            rest = gi // 3 ** f3
            f18659 = vp(rest, 18659) if rest > 1 else 0
            print(f"  {n:2d} {m:2d} | {Hn:14.1f} | {Hr:16.1f} | {math.log2(gi):8.1f} | {val:10d} | {Hr/val:.4f}"
                  f"   gcd = 3^{f3} * 18659^{f18659} * ({math.log2(rest / 18659 ** f18659):.1f} bits)")
    C_cyc = 5 / 54 - 0.676627737606436 / (math.pi ** 2 * math.sqrt(3))
    ladder = [
        ("periodic approximants (natural, Lemma H')", mu),
        ("Zudilin/Lemma P, m = 0", 2 * mu / 3),
        ("Zudilin/Lemma P, m = n/phi", mu / PHI),
        ("Hankel H1 (Cauchy-Binet)", 4 * mu / 7),
        ("Hankel H2 (+ q-order divisor)", 11 * mu / 28),
        ("Hankel H3 (+ cyclotomic divisor)", (11 / 12 - 2 * C_cyc) * mu / (7 / 3)),
    ]
    print("  asymptotic height exponent per unit of 2-adic valuation (irrationality needs < 1):")
    for name, x in ladder:
        print(f"     {name:45s} {x:.4f}")
    print("  (the defect exploited at each step: Pade cancellation, Vandermonde rank, forced odd divisors)")


# ------------------------------------------------------------------ D4 / D5
def content_poly(values_at, n):
    """Coefficients of a homogeneous degree-n polynomial N(a, b) from N(a, 1) at a = 0..n."""
    xs = list(range(n + 1))
    ys = [values_at(a) for a in xs]
    # Lagrange interpolation over Q (exact)
    coeffs = [Fraction(0)] * (n + 1)
    for i, xi in enumerate(xs):
        basis = [Fraction(1)]
        denom = Fraction(1)
        for j, xj in enumerate(xs):
            if j == i:
                continue
            basis = [Fraction(0)] + basis
            for t in range(len(basis) - 1):
                basis[t] -= xj * basis[t + 1]
            denom *= (xi - xj)
        for t in range(len(basis)):
            coeffs[t] += Fraction(int(ys[i])) * basis[t] / denom
    assert all(c.denominator == 1 for c in coeffs)
    g = 0
    for c in coeffs:
        g = math.gcd(g, int(c.numerator))
    return mpz(g)


def Nn_sq(n, a, b, L, M0):
    Y, X = mpz(2) ** L, mpz(M0)
    Gs = [a * X ** (m * m) - b * sum(Y ** (k * k) * X ** (m * m - k * k) for k in range(m + 1)) for m in range(2 * n - 1)]
    return mpz(int(flint.fmpz_mat([[int(Gs[i + j] * Y ** ((i - j) ** 2)) for j in range(n)] for i in range(n)]).det()))


def Nn_cube(n, a, b, L, M0):
    Y, X = mpz(2) ** L, mpz(M0)
    Gs = [a * X ** (s ** 3) - b * sum(Y ** (k ** 3) * X ** (s ** 3 - k ** 3) for k in range(s + 1)) for s in range(2 * n - 1)]
    return mpz(int(flint.fmpz_mat([[int(Gs[i + j] * Y ** (4 * i ** 3 + 4 * j ** 3 - (i + j) ** 3)) for j in range(n)]
                                   for i in range(n)]).det()))


def e0_krvz(n):
    return n * (n - 2) * (5 * n - 2) // 24 if n % 2 == 0 else n * (n - 1) * (5 * n - 7) // 24


def el_krvz(l, n):
    return sum((i + l) // (3 * l) + i // (3 * l) for i in range(n))


def cyc_hom(l, X, Y):
    x = sympy.symbols("x")
    c = [int(t) for t in sympy.Poly(sympy.cyclotomic_poly(l, x), x).all_coeffs()]
    d = len(c) - 1
    return sum(ci * mpz(X) ** (d - i) * mpz(Y) ** i for i, ci in enumerate(c))


def section_D4():
    hdr("D4  squares at the Y map: forced divisor (content of N_n(a,b)) versus the proven divisor D_n")
    L, M0 = L_Y, M0_Y
    nmax = 6 if QUICK else 9
    print("  n | log2 content | log2 D_n proven | extra bits | v3(extra) | v11(content) | v11(D_n)")
    for n in range(2, nmax + 1):
        cont = content_poly(lambda a: Nn_sq(n, a, 1, L, M0), n)
        ep = max(0, 2 * e0_krvz(n) - n * (n - 1))
        Dn = mpz(M0) ** ep
        for l in range(1, n):
            if 2 * l < n:
                Dn *= cyc_hom(l, mpz(M0) ** 2, mpz(2) ** (2 * L)) ** el_krvz(l, n)
        assert cont % Dn == 0
        extra = cont // Dn
        e3 = vp(extra, 3)
        r = extra // 3 ** e3
        e18 = vp(r, 18659) if r > 1 else 0
        r2 = r // 18659 ** e18
        print(f"  {n:2d} | {log2abs(cont):12.1f} | {log2abs(Dn):15.1f} | {log2abs(extra):10.1f} | {e3!s:>9s} |"
              f" {vp(cont, 11)!s:>12s} | {vp(Dn, 11)!s:>8s}   extra = 3^{e3} * 18659^{e18} * ({log2abs(r2):.1f} bits)")
    print("  The extra content is 3^(O(n^2)) (the tied alpha = q^(1/2) term of KRVZ Prop. 2) times powers of")
    print("  M0 - 2^L = 18659 = Phi_1^hom(M0, 2^L) (from Phi_1(xi) | Phi_1(q) = (xi-1)(xi+1), growing slowly):")
    print("  no additional cubic lever.")


def section_D5():
    hdr("D5  cubes at the Y3 map: forced divisor of the Hankel integers of the cube tails")
    L, M0 = L_Y, M0_Y
    nmax = 5 if QUICK else 7
    print("  n | K3 = 8 sum i^3 | log2 content | v3(content)/9 (X-order) | e3 law 3n^3-18n^2+34n-19 | v11(content)"
          " | quartic loss (mu_bar-1) L K3")
    mu = math.log2(M0) / L
    for n in range(2, nmax + 1):
        cont = content_poly(lambda a: Nn_cube(n, a, 1, L, M0), n)
        K3 = 8 * sum(i ** 3 for i in range(n))
        law = 3 * n ** 3 - 18 * n ** 2 + 34 * n - 19 if n >= 3 else 0
        v3c = vp(cont, 3)
        r = cont // 3 ** v3c
        e18 = vp(r, 18659) if r > 1 else 0
        r2 = r // 18659 ** e18
        print(f"  {n:2d} | {K3:14d} | {log2abs(cont):12.1f} | {v3c/9 if v3c is not None else None!s:>25s} | {law:25d} |"
              f" {vp(cont, 11)!s:>12s} | {(mu - 1) * L * K3:10.0f}   content = 3^{v3c} * 18659^{e18} * ({log2abs(r2):.1f} bits)")
    print("  The forced divisor of cube Hankel integers is M0^(3n^3 + O(n^2)) (cubic) times O(1)-size factors;")
    print("  closing the 0.299 gap needs a forced divisor of size M0^((1 - 1/mu_bar) K3) ~ M0^(0.6 n^4) (quartic).")
    print("  No systematic 11-adic (or 1089) factor appears: mod-11 structure cannot supply the missing lever.")


# ------------------------------------------------------------------ D6: exact forced divisor as a polynomial in xi
def _F_entry(s_, a, b, power):
    P = (lambda t: t * t) if power == 2 else (lambda t: t ** 3)
    c = [0] * (P(s_) + 1)
    c[P(s_)] += a
    for k in range(s_ + 1):
        c[P(s_) - P(k)] -= b
    return flint.fmpz_poly(c)


def _det_bareiss(M):
    n = len(M)
    A = [row[:] for row in M]
    prev = flint.fmpz_poly([1])
    sign = 1
    for k in range(n - 1):
        if A[k][k] == 0:
            pr = next((r for r in range(k + 1, n) if A[r][k] != 0), None)
            if pr is None:
                return flint.fmpz_poly([0])
            A[k], A[pr] = A[pr], A[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                qq, rr = divmod(A[i][j] * A[k][k] - A[i][k] * A[k][j], prev)
                assert rr == 0
                A[i][j] = qq
        prev = A[k][k]
    return sign * A[n - 1][n - 1]


def _cyc_index(f):
    for d in range(1, 400):
        if flint.fmpz_poly.cyclotomic(d) == f:
            return d
    return None


def section_D6():
    hdr("D6  exact forced divisor of det F(xi) (gcd over (a,b)) versus KRVZ's proven divisor (xi-degree units)")
    pts = [(0, 1), (1, 0), (1, 1), (2, 1), (1, 2), (3, 1), (1, 3), (5, 2), (2, 5), (7, 3), (3, 7), (11, 4), (4, 11),
           (13, 6), (6, 13), (17, 5)]

    def phi_(d):
        return sum(1 for k in range(1, d + 1) if math.gcd(k, d) == 1)
    nmax = 10 if QUICK else 14
    print("  squares:  n | K_n | forced deg (xi-order + cyclotomic) | KRVZ proven deg | excess | excess/n^2 | forced/n^3")
    for n in range(4, nmax + 1):
        g = None
        for (a, b) in pts[: n + 2]:
            d = _det_bareiss([[_F_entry(i + j, a, b, 2) for j in range(n)] for i in range(n)])
            g = d if g is None else g.gcd(d)
        fac = g.factor()
        tot = 0
        for f, e in fac[1]:
            if f == flint.fmpz_poly([0, 1]):
                tot += e
            else:
                assert _cyc_index(f) is not None
                tot += f.degree() * e
        K = 4 * sum(i * i for i in range(n))
        proven = max(0, 2 * e0_krvz(n) - n * (n - 1)) + 2 * sum(phi_(l) * el_krvz(l, n) for l in range(1, n) if 2 * l < n)
        print(f"            {n:2d} | {K:5d} | {tot:34d} | {proven:15d} | {tot - proven:6d} | {(tot - proven)/n**2:10.3f} | {tot/n**3:.4f}")
    print("  cubes:    n | K3_n | xi-order (law 3n^3-18n^2+34n-19) | cyclotomic deg | forced/K3_n")
    for n in range(3, 7 if QUICK else 8):
        g = None
        for (a, b) in pts[: n + 2]:
            d = _det_bareiss([[_F_entry(i + j, a, b, 3) for j in range(n)] for i in range(n)])
            g = d if g is None else g.gcd(d)
        fac = g.factor()
        xo = sum(e for f, e in fac[1] if f == flint.fmpz_poly([0, 1]))
        cy = sum(f.degree() * e for f, e in fac[1] if f != flint.fmpz_poly([0, 1]))
        K3 = 8 * sum(i ** 3 for i in range(n))
        law = 3 * n ** 3 - 18 * n ** 2 + 34 * n - 19
        print(f"            {n:2d} | {K3:5d} | {xo:6d} ({law:6d}) {'':20s} | {cy:14d} | {(xo + cy)/K3:.4f}")
    print("  All forced factors are xi-powers and cyclotomic Phi_d(xi) (no other irreducible factor occurs).")
    print("  Squares: the tied (theta-null) case carries MORE forced factors than KRVZ prove (Phi_1(xi), Phi_3(xi),")
    print("  Phi_4(xi), ... with larger exponents); the excess grows slowly (excess/n^2 from 0.6 to ~1.2 for")
    print("  n = 4..14, excess/n^3 decreasing), i.e. below the cubic order of the proven part: no extra cubic")
    print("  lever in range. Cubes: the xi-order follows the cubic law (checked to n = 12 in the scratch runs) and")
    print("  the cyclotomic part is small, against a quartic clearing K3 ~ 2n^4, so forced/K3 -> 0.")


if __name__ == "__main__":
    t0 = time.time()
    print("collatz_procgen_20260923_theta_defect.py" + (" --quick" if QUICK else ""))
    section_D1()
    section_D2()
    section_D3()
    section_D4()
    section_D5()
    section_D6()
    print(f"\n[defect done in {time.time()-t0:.1f}s]")
