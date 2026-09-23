#!/usr/bin/env python3
"""
collatz_procgen_20260923_theta_euler.py

Euler's function P(rho) = (rho; rho)_inf = prod_(n>=1) (1 - rho^n) = sum_(k in Z) (-1)^k rho^(k(3k-1)/2)
in Q_2, rho = 2^L/M0: 2-adic transcription of the Hankel method for the q-exponential case of
KRVZ (lambda = 1):  P = F_q(alpha; lambda) = sum_n alpha^n / prod_(j<=n)(q^j - lambda),
with q = 1/rho = M0/2^L, alpha = -1, lambda = 1.

  E0  three representations of P agree mod 2^N (product, Euler's q-series, pentagonal bilateral sum)
  E1  v2(V_n) = L n(n+1)(2n+1)/6 for the Hankel determinant of the KRVZ tails v_n (several maps)
  E2  integer N_n = 2^(L D_n) b^n V_n(q = M0/2^L, alpha = -1, mu = a/b), D_n = n(n-1)(4n+1)/6:
      identity check (fake rationals), v2(N_n) = L n^3, archimedean bound
  E3  divisibility: M0^(e0(n)) (KRVZ Prop. 1, lambda != 0: e0 = n(n-1)(n-2)/6) and
      Phi_l^hom(M0, 2^L)^(e_l(n)) (KRVZ Prop. 4) divide N_n
  E4  margins for rho = 2^10/3^9 and other maps; thresholds 3/2 (elementary), 2, 2.237193
  E5  realizability: three blocks B, B', B'' of equal length and weight with R_B'' = 2 R_B - R_B'
      (signed bilateral pentagonal swap words) under 3x+1
  E6  remark: the unsigned bilateral pentagonal sum via the two-family Hankel (naive threshold 21/16)

Usage: python3 collatz_procgen_20260923_theta_euler.py [--quick]
"""
import itertools
import math
import random
import sys
import time
from fractions import Fraction

import flint
import gmpy2
from gmpy2 import mpz
import sympy

QUICK = "--quick" in sys.argv


def hdr(s):
    print()
    print("=" * 78)
    print(s)
    print("=" * 78)
    sys.stdout.flush()


def v2(x):
    x = mpz(x)
    return int(gmpy2.bit_scan1(abs(x))) if x else None


def log2abs(x):
    x = abs(mpz(x))
    b = int(x.bit_length())
    if b <= 1000:
        return math.log2(int(x))
    sh = b - 60
    return math.log2(int(x >> sh)) + sh


def e0_lam(n):          # KRVZ (3.14), lambda != 0
    return n * (n - 1) * (n - 2) // 6


def el_krvz(l, n):      # KRVZ (4.1)
    return sum((i + l) // (3 * l) + i // (3 * l) for i in range(n))


def Dn(n):
    return n * (n - 1) * (4 * n + 1) // 6


_cyc = {}


def cyc_hom(l, X, Y):
    if l not in _cyc:
        x = sympy.symbols("x")
        _cyc[l] = [int(c) for c in sympy.Poly(sympy.cyclotomic_poly(l, x), x).all_coeffs()]
    c = _cyc[l]
    d = len(c) - 1
    return sum(ci * mpz(X) ** (d - i) * mpz(Y) ** i for i, ci in enumerate(c))


MAPS = [
    ("3x+1 (A=9,L=10): rho = 2^10/3^9", 10, 3 ** 9),
    ("3x+1 A=L (rho = 2/3)", 1, 3),
    ("5x+1 1^6 0^4", 10, 5 ** 6),
    ("5x+1 A=L (rho = 2/5)", 1, 5),
]


# ------------------------------------------------------------------ 2-adic helpers
def euler_mod(N, L, M0):
    mod = mpz(1) << N
    inv = gmpy2.invert(mpz(M0), mod)
    rho = ((mpz(1) << L) * inv) % mod if L < N else mpz(0)
    # product
    prod = mpz(1)
    n = 1
    while L * n < N:
        prod = (prod * (1 - gmpy2.powmod(rho, n, mod))) % mod
        n += 1
    # pentagonal bilateral sum
    pent = mpz(0)
    k = 0
    while True:
        e1 = k * (3 * k - 1) // 2
        if L * e1 >= N:
            break
        pent += (-1) ** k * gmpy2.powmod(rho, e1, mod)
        if k > 0:
            e2 = k * (3 * k + 1) // 2
            if L * e2 < N:
                pent += (-1) ** k * gmpy2.powmod(rho, e2, mod)
        k += 1
    pent %= mod
    # Euler q-series sum (-1)^k rho^(k(k+1)/2)/(rho;rho)_k
    ser = mpz(0)
    poch = mpz(1)
    k = 0
    while L * (k * (k + 1) // 2) < N:
        if k > 0:
            poch = (poch * (1 - gmpy2.powmod(rho, k, mod))) % mod
        ser += (-1) ** k * gmpy2.powmod(rho, k * (k + 1) // 2, mod) * gmpy2.invert(poch, mod)
        k += 1
    ser %= mod
    return prod, pent, ser, rho, mod


def tails_euler(N, L, M0, count):
    """KRVZ tails v_n = (-1)^n sum_m (-1)^m rho^(mn + m(m+1)/2) / prod_(j=n+1..n+m)(1 - rho^j) mod 2^N."""
    mod = mpz(1) << N
    inv = gmpy2.invert(mpz(M0), mod)
    rho = ((mpz(1) << L) * inv) % mod
    out = []
    for n in range(count):
        t = mpz(0)
        m = 1
        den = mpz(1)
        while L * (m * n + m * (m + 1) // 2) < N:
            den = (den * (1 - gmpy2.powmod(rho, n + m, mod))) % mod
            t += (-1) ** m * gmpy2.powmod(rho, m * n + m * (m + 1) // 2, mod) * gmpy2.invert(den, mod)
            m += 1
        out.append(((-1) ** n * t) % mod)
    return out


def Nn_integer(n, a, b, L, M0):
    """N_n = 2^(L D_n) b^n V_n(q = M0/2^L, alpha = -1, mu = a/b) via exact polynomials in q."""
    q = sympy.symbols("q")
    mu_num, mu_den = a, b
    # b * v_m as polynomial in q with integer coefficients
    vs = []
    prodq = sympy.Integer(1)
    for m in range(2 * n - 1):
        if m > 0:
            prodq = sympy.expand(prodq * (q ** m - 1))
        # b v_m = a prod_{j<=m}(q^j-1) - b sum_k (-1)^k prod_{j=k+1}^m (q^j - 1)
        tail = sympy.Integer(0)
        for k in range(m + 1):
            p = sympy.Integer(1)
            for j in range(k + 1, m + 1):
                p = p * (q ** j - 1)
            tail += (-1) ** k * p
        vs.append(sympy.Poly(sympy.expand(mu_num * prodq - mu_den * tail), q))
    X, Y = mpz(M0), mpz(2) ** L
    D = Dn(n)
    # evaluate each entry homogeneously later: compute det of the polynomial matrix numerically at q = X/Y
    # using exact rationals, then scale by 2^(L D).
    vals = [Fraction(int(sum(int(c) * X ** (len(P.all_coeffs()) - 1 - i) * Y ** i
                                for i, c in enumerate(P.all_coeffs()))),
                     int(Y ** (len(P.all_coeffs()) - 1))) for P in vs]
    M = flint.fmpq_mat([[flint.fmpq(vals[i + j].numerator, vals[i + j].denominator) for j in range(n)] for i in range(n)])
    d = M.det()
    d = Fraction(int(d.p), int(d.q))
    Nn = d * Fraction(int(Y)) ** D
    assert Nn.denominator == 1, "N_n not an integer"
    return mpz(Nn.numerator)


# ------------------------------------------------------------------ sections
def section_E0():
    hdr("E0  (rho;rho)_inf: product = Euler q-series = pentagonal bilateral sum (mod 2^N)")
    for name, L, M0 in MAPS:
        N = 4000
        prod, pent, ser, rho, mod = euler_mod(N, L, M0)
        print(f"  {name:34s} product==pentagonal: {prod == pent};  product==q-series: {prod == ser}")


def section_E1():
    hdr("E1  v2(V_n) of the Hankel determinant of the 2-adic tails: prediction L n(n+1)(2n+1)/6")
    nmax = 7 if QUICK else 12
    for name, L, M0 in MAPS:
        ok = True
        got_list = []
        for n in range(1, nmax + 1):
            pred = L * n * (n + 1) * (2 * n + 1) // 6
            N = pred + 20 * L + 64
            tl = tails_euler(N, L, M0, 2 * n)
            H = flint.fmpz_mat([[int(tl[i + j]) for j in range(n)] for i in range(n)])
            got = v2(int(H.det()) % (1 << N))
            got_list.append(got)
            ok = ok and (got == pred)
        print(f"  {name:34s} n=1..{nmax}: {got_list[:6]}...  all == prediction: {ok}")


def section_E2_E3():
    hdr("E2/E3  integer N_n: identity, v2 = L n^3, archimedean bound; divisibility (KRVZ Props. 1, 4)")
    random.seed(11)
    nmax = 5 if QUICK else 7
    for name, L, M0 in MAPS[:3]:
        okint = True
        okdiv = True
        worst = -1e9
        for trial in range(3):
            a = random.randint(-10 ** 5, 10 ** 5)
            b = random.choice([1, 3, 5, 7, 11]) * random.choice([1, -1])
            for n in range(2, nmax + 1):
                Nn = Nn_integer(n, a, b, L, M0)
                Dv = mpz(M0) ** e0_lam(n)
                for l in range(1, n):
                    if 2 * l < n:
                        Dv *= cyc_hom(l, M0, 2 ** L) ** el_krvz(l, n)
                if Nn % Dv != 0:
                    okdiv = False
                C = abs(a) + 2 * abs(b)
                bound = Dn(n) * math.log2(M0) + n * (n - 1) + n * math.log2(C) + math.lgamma(n + 1) / math.log(2)
                if Nn != 0:
                    worst = max(worst, log2abs(Nn) - bound)
        print(f"  {name:34s} N_n integer: {okint};  M0^e0 * prod Phi_l^hom divides N_n: {okdiv};"
              f"  max log2|N_n| - bound = {worst:.1f} (<= 0)")


def margin(L, M0, n, level):
    m = L * n ** 3 - Dn(n) * math.log2(M0) - n * (n - 1) - math.lgamma(n + 1) / math.log(2)
    if level >= 2:
        m += e0_lam(n) * math.log2(M0)
    if level >= 3:
        for l in range(1, n):
            if 2 * l < n:
                m += el_krvz(l, n) * log2abs(cyc_hom(l, M0, 2 ** L))
    return m


def section_E4():
    hdr("E4  margins E_k(n) = v2(N_n) + log2(odd divisor) - log2(archimedean bound w/o the C^n term)")
    C_cyc = 5 / 54 - 0.676627737606436 / (math.pi ** 2 * math.sqrt(3))
    c1, c2, c3 = 1.5, 2.0, 1 / (0.5 - C_cyc)
    print(f"  thresholds: E1 (elementary) mu_bar < 3/2;  E2 (+KRVZ Prop. 1) < 2;  E3 (+KRVZ Prop. 4) < {c3:.6f}")
    print(f"  real-case dictionary: KRVZ Thm 2 (lambda != 0) gamma > 1.80828115  <->  mu_bar < {1.80828115/0.80828115:.6f}")
    ns = [5, 10, 20, 40] if QUICK else [5, 10, 20, 40, 80, 160]
    for name, L, M0 in MAPS:
        mu = math.log2(M0) / L
        lead = [L - (2 / 3) * math.log2(M0), L - 0.5 * math.log2(M0), L - (0.5 - C_cyc) * math.log2(M0)]
        print(f"  {name:34s} mu_bar={mu:.5f} leading n^3 coeff: E1 {lead[0]:+.4f}  E2 {lead[1]:+.4f}  E3 {lead[2]:+.4f}")
        for level in (1, 2, 3):
            print(f"      E{level}(n), n={ns}: " + " ".join(f"{margin(L, M0, n, level):.4g}" for n in ns))
    L, M0 = 10, 3 ** 9
    print("  rho = 2^10/3^9, elementary level: least n with E1(n) > n log2(3 * 2^h) (excludes all a/b of height <= 2^h):")
    for h in (100, 10 ** 4, 10 ** 6):
        n = 2
        while not margin(L, M0, n, 1) > n * (math.log2(3) + h):
            n += 1 if n < 200 else max(1, n // 100)
        print(f"      h = {h}: n = {n}")


def section_E5():
    hdr("E5  signed bilateral pentagonal swap words: blocks with R_B'' = 2 R_B - R_B' (3x+1)")
    def R(z):
        r = 0
        for i, ch in enumerate(z):
            if ch == "1":
                r = 3 * r + 2 ** i
        return r
    found = []
    for L in range(3, 15 if not QUICK else 12):
        for A in range(1, L + 1):
            if A / L <= math.log(2) / math.log(3):
                continue
            blocks = ["".join("1" if i in S else "0" for i in range(L)) for S in itertools.combinations(range(L), A)]
            if len(blocks) > 3000:
                continue
            rv = {}
            for z in blocks:
                rv.setdefault(R(z), []).append(z)
            vals = sorted(rv)
            sv = set(vals)
            hit = None
            for i, x in enumerate(vals):
                for y in vals[i + 1:]:
                    if 2 * y - x in sv and 2 * y - x != y:
                        hit = (rv[x][0], rv[y][0], rv[2 * y - x][0])
                        break
                if hit:
                    break
            if hit:
                found.append((L, A, A * math.log2(3) / L, hit))
    for L, A, mu, (b1, b2, b3) in found[:12]:
        print(f"  L={L:2d} A={A:2d} mu_bar={mu:.4f}: B'={b1}  B={b2}  B''={b3}  (R in arithmetic progression)")
    if not any(L == 10 and A == 9 for L, A, _, _ in found):
        print("  (L, A) = (10, 9): no arithmetic progression of R-values among the 10 blocks with nine 1s, so no")
        print("  signed-pentagonal swap word has exactly rho = 2^10/3^9; the NUMBER (rho;rho)_inf is still proved")
        print("  irrational by Theorem E1, and every listed (L, A) gives signed pentagonal swap words covered by E1/E2.")
    return found


def section_E6():
    hdr("E6  unsigned bilateral pentagonal sum_(k in Z) rho^(k(3k-1)/2): two-family Hankel (naive)")
    L, M0 = 10, 3 ** 9
    g = []
    e = []
    for m in range(1, 40):
        g += [m * (3 * m - 1) // 2, m * (3 * m + 1) // 2]
        e += [3 * m, 3 * m + 1]
    ok = True
    nmax = 6 if QUICK else 9
    for n in range(1, nmax + 1):
        pred = L * (sum(g[:n]) + 2 * sum((n - 1 - i) * e[i] for i in range(n)))
        N = pred + 20 * L + 64
        mod = mpz(1) << N
        inv = gmpy2.invert(mpz(M0), mod)
        rho = ((mpz(1) << L) * inv) % mod
        tl = []
        for s in range(2 * n):
            t = mpz(0)
            m = 1
            while L * (m * (3 * m - 1) // 2 + 3 * m * s) < N:
                t += (-1) ** m * (gmpy2.powmod(rho, m * (3 * m - 1) // 2 + 3 * m * s, mod)
                                  + gmpy2.powmod(rho, m * (3 * m + 1) // 2 + (3 * m + 1) * s, mod))
                m += 1
            tl.append(t % mod)
        H = flint.fmpz_mat([[int(tl[i + j]) for j in range(n)] for i in range(n)])
        got = v2(int(H.det()) % (1 << N))
        ok = ok and (got == pred)
    print(f"  signed two-family tails (= Euler's function): v2(T_n) = L[sum g_i + 2 sum (n-i) e_i] for n <= {nmax}: {ok}")
    print("  leading terms (5/8) L n^3 against a clearing 2 n^3 log2 M0: naive threshold mu_bar < 21/16 = 1.3125;")
    print("  the q-exponential representation (E1) is strictly better (3/2), and the unsigned sum has no")
    print("  single q-exponential representation (it is a product of three q-Pochhammer symbols).")


if __name__ == "__main__":
    t0 = time.time()
    print("collatz_procgen_20260923_theta_euler.py" + (" --quick" if QUICK else ""))
    section_E0()
    section_E1()
    section_E2_E3()
    section_E4()
    section_E5()
    section_E6()
    print(f"\n[euler done in {time.time()-t0:.1f}s]")
