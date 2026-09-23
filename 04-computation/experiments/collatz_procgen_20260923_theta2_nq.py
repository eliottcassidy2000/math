#!/usr/bin/env python3
"""
collatz_procgen_20260923_theta2_nq.py

HYP-9132 (partial): a 2-adic transcription of the NON-QUADRATICITY argument of
Krattenthaler-Rochev-Vaananen-Zudilin [KRVZ, Theorem 1 and section 5].

Theorem NQ (note section 1). Let rho = 2^L/M0 with M0 > 2^L odd, and mu_bar = log2(M0)/L. If
    mu_bar < c_NQ = 126 pi^2 / (79 pi^2 + 72 sqrt3 Im Li_2(e^(2 pi i/3))) = 1.43918...,
then Theta(rho) = sum_(k>=0) rho^(k^2) in Z_2 is neither rational nor quadratic over Q.
rho = 2^10/3^9 has mu_bar = 1.42647, so theta_3(2^10/3^9) = 2 Theta - 1 has degree >= 3 (or is transcendental).

Mechanism: if Theta = alpha/w with alpha in O_K, [K:Q] = 2, then K embeds in Q_2, so 2 splits in K.
N_n = det(G_(i+j) 2^(L(i-j)^2)) (G_m = alpha M0^(m^2) - w sum_(k<=m) 2^(L k^2) M0^(m^2-k^2)) lies in O_K.
Under the embedding iota_1 with iota_1(alpha)/w = Theta, N_n = 2^(L K_n) w^n T_n, so
v_2(iota_1 N_n) = L K_n + n v_2(w) + L n(n+1)(2n-1)/2. The other 2-adic embedding gives v_2 >= 0.
The forced odd divisor G_n (KRVZ Props. 2 and 4, transferred in the round-2 note, Lemma D) divides N_n,
so G_n^2 divides the nonzero integer Norm(N_n), and |Norm(N_n)| <= (n! C^n M0^(K_n))^2.

  Q0  thresholds (theta and Euler; irrationality d = 1 and non-quadraticity d = 2) and the KRVZ dictionary
  Q1  random alpha in Z[sqrt 17], Z[sqrt -7]: G_n^2 | Norm(N_n), archimedean bound (n <= 5)
  Q2  Theta-approximants alpha/w in Q(sqrt 17), Q(sqrt -7) (LLL, 2-adic precision 2^4000):
      v_2(iota_1 N_n) equals the Hankel valuation exactly; v_2(Norm) >= it; G_n^2 | Norm; Norm != 0
  Q3  the exact margin A_NQ(n) at rho = 2^10/3^9; the n from which it stays positive; heights
  Q4  Euler's function (lambda = 1): the same transcription needs mu_bar < 1.1186 (not reached)

Usage: python3 collatz_procgen_20260923_theta2_nq.py [--quick]
"""
import math
import random
import sys
import time

import flint
import gmpy2
from gmpy2 import mpz
import mpmath
import sympy

QUICK = "--quick" in sys.argv
I_LI2 = float(mpmath.polylog(2, mpmath.exp(2j * mpmath.pi / 3)).imag)
C_CYC = 5 / 54 - I_LI2 / (math.pi ** 2 * math.sqrt(3))


def hdr(s):
    print()
    print("=" * 78)
    print(s)
    print("=" * 78)
    sys.stdout.flush()


def v2(x):
    x = mpz(x)
    return 10 ** 9 if x == 0 else int(gmpy2.bit_scan1(abs(x)))


# ---------------------------------------------------------------- KRVZ exponents (round-2 Lemma D)
def e0_krvz(n):
    return n * (n - 2) * (5 * n - 2) // 24 if n % 2 == 0 else n * (n - 1) * (5 * n - 7) // 24


def eprime(n):
    return max(0, 2 * e0_krvz(n) - n * (n - 1))


def _S(M, m):
    """sum_(j=0)^(M-1) floor(j/m)."""
    q, r = divmod(M, m)
    return m * q * (q - 1) // 2 + r * q


def el_krvz(l, n):
    """e_l(n) = sum_(i<n) (floor((i+l)/3l) + floor(i/3l)), closed form."""
    m = 3 * l
    return _S(n + l, m) - _S(l, m) + _S(n, m)


def el_direct(l, n):
    return sum((i + l) // (3 * l) + i // (3 * l) for i in range(n))


def Kn(n):
    return 4 * sum(i * i for i in range(n))


def v2T(L, n):
    return L * n * (n + 1) * (2 * n - 1) // 2


_cyc = {}


def cyc_hom(l, X, Y):
    """Phi_l^hom(X, Y) = X^phi(l) Phi_l(Y/X) (up to sign for l = 1)."""
    if l not in _cyc:
        _cyc[l] = [int(c) for c in flint.fmpz_poly.cyclotomic(l).coeffs()]  # low -> high
    c = _cyc[l]
    d = len(c) - 1
    return sum(ci * mpz(Y) ** i * mpz(X) ** (d - i) for i, ci in enumerate(c))


def forced_odd(n, L, M0):
    G = mpz(M0) ** eprime(n)
    for l in range(1, n):
        if 2 * l < n:
            G *= cyc_hom(l, mpz(M0) ** 2, mpz(2) ** (2 * L)) ** el_krvz(l, n)
    return G


# ---------------------------------------------------------------- Q0
def section_Q0():
    hdr("Q0  thresholds: irrationality (d = 1) and non-quadraticity (d = 2), 2-adic transcription of KRVZ")
    print(f"  C_cyc = 5/54 - Im Li2(e^(2 pi i/3))/(pi^2 sqrt3) = {C_CYC:.8f}   (Im Li2 = {I_LI2:.12f})")
    print("  per n^3: 2-adic gain (7/3) L; clearing 4/3 log2 M0 per archimedean conjugate; forced divisor kappa log2 M0")
    for name, kappa in (("no divisor", 0.0), ("+ q-order (Prop. 2)", 5 / 12), ("+ cyclotomic (Prop. 4)", 5 / 12 + 2 * C_CYC)):
        c1 = (7 / 3) / (4 / 3 - kappa)
        c2 = (7 / 3) / (2 * (4 / 3 - kappa))
        c3 = (7 / 3) / (3 * (4 / 3 - kappa))
        print(f"  theta, {name:24s}: d=1 mu_bar < {c1:.5f};  d=2 mu_bar < {c2:.5f};  d=3 mu_bar < {c3:.5f}")
    cNQ = 126 * math.pi ** 2 / (79 * math.pi ** 2 + 72 * math.sqrt(3) * I_LI2)
    gK = 126 * math.pi ** 2 / (47 * math.pi ** 2 - 72 * math.sqrt(3) * I_LI2)
    print(f"  closed form c_NQ = 126 pi^2/(79 pi^2 + 72 sqrt3 Im Li2) = {cNQ:.6f}")
    print(f"  dictionary: KRVZ Thm 1 (lambda = 0) gamma > {gK:.8f} (they print 3.27694460); gamma/(gamma-1) = {gK/(gK-1):.6f}")
    print("  d = 3 would need mu_bar < 0.9595 < 1: as in KRVZ (5.11), only d = 1, 2 are reachable.")
    for name, kappa in (("no divisor", 0.0), ("+ q-order (Prop. 1)", 1 / 6), ("+ cyclotomic (Prop. 4)", 1 / 6 + C_CYC)):
        c1 = 1 / (2 / 3 - kappa)
        c2 = 1 / (2 * (2 / 3 - kappa))
        print(f"  Euler, {name:24s}: d=1 mu_bar < {c1:.5f};  d=2 mu_bar < {c2:.5f}")
    gE = 27 * math.pi ** 2 / (5 * math.pi ** 2 - 18 * math.sqrt(3) * I_LI2)
    print(f"  dictionary: KRVZ Thm 1 (lambda != 0) gamma > {gE:.8f} (they print 9.43194241); gamma/(gamma-1) = {gE/(gE-1):.6f}")
    for lab, L, M0 in (("3x+1 Y3 rho = 2^10/3^9", 10, 3 ** 9), ("5x+1 HYP-9131 rho = 2^33/5^23", 33, 5 ** 23),
                       ("3x+1 2^8/3^7", 8, 3 ** 7), ("5x+1 2^5/5^3", 5, 5 ** 3), ("7x+1 2^2/7", 2, 7)):
        mu = math.log2(M0) / L
        print(f"  {lab:32s}: mu_bar = {mu:.5f} -> theta not quadratic: {mu < cNQ}")
    return cNQ


# ---------------------------------------------------------------- exact N_n in Z[sqrt d]
def _bareiss(M):
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
                q, r = divmod(A[i][j] * A[k][k] - A[i][k] * A[k][j], prev)
                assert r == 0
                A[i][j] = q
        prev = A[k][k]
    return sign * A[n - 1][n - 1]


def Nn_quadratic(n, u, v, w, d, L, M0):
    """N_n(alpha, w) for alpha = u + v sqrt(d): returns (x, y) with N_n = x + y sqrt(d)."""
    X, Y = mpz(M0), mpz(2) ** L
    alpha = flint.fmpz_poly([int(u), int(v)])        # u + v s, s^2 = d

    def G(m):
        return alpha * int(X ** (m * m)) - int(w) * int(sum(Y ** (k * k) * X ** (m * m - k * k) for k in range(m + 1)))
    Gs = [G(m) for m in range(2 * n - 1)]
    M = [[Gs[i + j] * int(Y ** ((i - j) ** 2)) for j in range(n)] for i in range(n)]
    det = _bareiss(M)
    cs = [int(c) for c in det.coeffs()]
    x = sum(c * d ** (k // 2) for k, c in enumerate(cs) if k % 2 == 0)
    y = sum(c * d ** (k // 2) for k, c in enumerate(cs) if k % 2 == 1)
    return mpz(x), mpz(y)


def sqrt2adic(d, N):
    """s with s^2 = d mod 2^N (d = 1 mod 8)."""
    assert d % 8 == 1
    s = mpz(1)
    for k in range(3, N):
        if (s * s - d) % (mpz(1) << (k + 1)) != 0:
            s += mpz(1) << (k - 1)
    s %= mpz(1) << N
    assert (s * s - d) % (mpz(1) << N) == 0
    return s


def theta_mod(L, M0, P):
    mod = mpz(1) << P
    inv = gmpy2.invert(mpz(M0), mod)
    t = mpz(0)
    k = 0
    while L * k * k < P:
        t += (mpz(1) << (L * k * k)) * gmpy2.powmod(inv, k * k, mod)
        k += 1
    return t % mod


def arch_bound_log2(n, Cm, M0):
    lc = float(gmpy2.log2(gmpy2.mpfr(Cm)))
    return math.lgamma(n + 1) / math.log(2) + n * lc + Kn(n) * math.log2(M0)


# ---------------------------------------------------------------- Q1
def section_Q1():
    hdr("Q1  random alpha = u + v sqrt(d): forced divisor and archimedean bound in O_K (2 splits: d = 17, -7)")
    random.seed(9132)
    L, M0 = 10, 3 ** 9
    nmax = 4 if QUICK else 5
    for d in (17, -7):
        ok_div, worst, cnt = True, -1e18, 0
        for trial in range(3):
            u, v, w = random.randint(-50, 50), random.randint(1, 20), random.choice([1, 2, 3, 5, 7])
            for n in range(2, nmax + 1):
                x, y = Nn_quadratic(n, u, v, w, d, L, M0)
                norm = x * x - d * y * y
                G = forced_odd(n, L, M0)
                ok_div &= (norm % (G * G) == 0)
                Cm = abs(u) + abs(v) * (gmpy2.isqrt(abs(d)) + 1) + (2 * n - 1) * w
                if norm != 0:
                    worst = max(worst, float(gmpy2.log2(gmpy2.mpfr(abs(norm)))) - 2 * arch_bound_log2(n, Cm, M0))
                cnt += 1
        print(f"  d = {d:3d}: {cnt} cases (n <= {nmax}); G_n^2 | Norm(N_n): {ok_div}; "
              f"max[log2|Norm| - 2 log2(n! C^n M0^K_n)] = {worst:.2f} (must be <= 0)")


# ---------------------------------------------------------------- Q2
def lll_approximant(Th, s, P):
    """short (u, v, w) with u + v s - w Th = 0 mod 2^P and v, w != 0."""
    m = int(mpz(1) << P)
    B = flint.fmpz_mat([[m, 0, 0], [int((-s) % m), 1, 0], [int(Th % m), 0, 1]])
    R = B.lll()
    best = None
    for i in range(3):
        u, v, w = (int(R[i, 0]), int(R[i, 1]), int(R[i, 2]))
        if w != 0 and v != 0:
            size = max(abs(u), abs(v), abs(w))
            if best is None or size < best[3]:
                best = (u, v, w, size)
    return best


def section_Q2():
    hdr("Q2  Theta-approximants in Q(sqrt 17) and Q(sqrt -7): one embedding carries the full Hankel valuation")
    L, M0 = 10, 3 ** 9
    P = 4000
    Th = theta_mod(L, M0, P)
    nmax = 4 if QUICK else 5
    print(f"  rho = 2^10/3^9, Theta mod 2^{P}; alpha = u + v sqrt(d) with iota_1(alpha)/w = Theta mod 2^(P - v2(w))")
    print("   d |  n | v2(iota_1 N_n) | pred. L K_n + v2(T_n) + n v2(w) | v2(iota_2 N_n) | v2(Norm) | G_n^2 | Norm"
          " | arch. bound")
    for d in (17, -7):
        s = sqrt2adic(d, P)
        u, v, w, size = lll_approximant(Th, s, P)
        vw = v2(w)
        err = v2(((u + v * s) - w * Th) % (mpz(1) << P))
        for n in range(2, nmax + 1):
            x, y = Nn_quadratic(n, u, v, w, d, L, M0)
            mod = mpz(1) << P
            i1 = (x + y * s) % mod
            i2 = (x - y * s) % mod
            norm = x * x - d * y * y
            pred = L * Kn(n) + v2T(L, n) + n * vw
            G = forced_odd(n, L, M0)
            Cm = abs(u) + abs(v) * (gmpy2.isqrt(abs(d)) + 1) + (2 * n - 1) * abs(w)
            arch = float(gmpy2.log2(gmpy2.mpfr(abs(norm)))) <= 2 * arch_bound_log2(n, Cm, M0) + 1e-9
            print(f"  {d:3d} | {n:2d} | {v2(i1):14d} | {pred:31d} | {v2(i2):14d} | {v2(norm):8d} | "
                  f"{str(norm % (G * G) == 0):5s} | {'!=0' if norm != 0 else '=0':4s} | {arch}")
        print(f"      (approximant: log2 max(|u|,|v|,|w|) = {float(gmpy2.log2(gmpy2.mpfr(size))):.1f}, v2(w) = {vw}, "
              f"v2(u + v s - w Theta) >= {min(err, P)})")
    print("  v2(iota_1 N_n) equals the prediction: the embedding that sees Theta carries the whole Hankel")
    print("  valuation, the other embedding contributes >= 0. This is the only 2-adic input of Theorem NQ.")


# ---------------------------------------------------------------- Q3
_logs = {}


def log2_cyc_hom_float(l, lX, t):
    """log2|Phi_l^hom(X^2, Y^2)| = 2 phi(l) log2 X + sum_(e|l) mu(l/e) log2|1 - t^e|, t = Y^2/X^2."""
    if l not in _logs:
        phi = int(sympy.totient(l))
        corr = 0.0
        for e in sympy.divisors(l):
            te = t ** e
            if te > 1e-300:
                corr += int(sympy.mobius(l // e)) * math.log2(abs(1 - te))
        _logs[l] = 2 * phi * lX + corr
    return _logs[l]


def A_NQ(n, L, M0):
    lX = math.log2(M0)
    t = (2.0 ** (2 * L)) / float(M0) ** 2
    m = L * Kn(n) + v2T(L, n)
    logG = eprime(n) * lX
    for l in range(1, n):
        if 2 * l < n:
            logG += el_krvz(l, n) * log2_cyc_hom_float(l, lX, t)
    return m + 2 * logG - 2 * (Kn(n) * lX + math.lgamma(n + 1) / math.log(2))


def section_Q3():
    hdr("Q3  exact margin A_NQ(n) = L K_n + v2(T_n) + 2 log2 G_n - 2 log2(n! M0^K_n) at rho = 2^10/3^9")
    L, M0 = 10, 3 ** 9
    ok = all(el_krvz(l, n) == el_direct(l, n) for n in range(1, 120) for l in range(1, n))
    print(f"  closed form of e_l(n) equals the defining sum for n < 120: {ok}")
    worst = 0.0
    for n in range(2, 26):
        Gi = forced_odd(n, L, M0)
        lf = eprime(n) * math.log2(M0) + sum(el_krvz(l, n) * log2_cyc_hom_float(l, math.log2(M0), 2.0 ** 20 / 3.0 ** 18)
                                               for l in range(1, n) if 2 * l < n)
        lg = float(gmpy2.log2(gmpy2.mpfr(abs(Gi)))) if abs(Gi) > 1 else 0.0
        worst = max(worst, abs(lg - lf))
    print(f"  float log2 G_n agrees with the exact integer to {worst:.2e} bits (n < 26)")
    lead = (7 / 3) * L - 2 * (4 / 3 - 5 / 12 - 2 * C_CYC) * math.log2(M0)
    print(f"  leading coefficient of A_NQ(n)/n^3: {lead:+.5f} bits (positive iff mu_bar < c_NQ)")
    ns = [10, 50, 100, 200, 400, 800, 1600] + ([] if QUICK else [3200, 6400])
    vals = [A_NQ(n, L, M0) for n in ns]
    print("  A_NQ(n):      " + ", ".join(f"n={n}: {v:.4g}" for n, v in zip(ns, vals)))
    print("  A_NQ(n)/n^3:  " + ", ".join(f"{v / n ** 3:.4f}" for n, v in zip(ns, vals)))
    top = 2500 if QUICK else 8000
    lastneg = None
    for n in range(2, top + 1):
        if A_NQ(n, L, M0) <= 0:
            lastneg = n
    print(f"  A_NQ(n) <= 0 exactly for n <= {lastneg} among n <= {top}; positive for {lastneg + 1} <= n <= {top}")
    for h in (10, 100, 1000, 10000):
        n = lastneg + 1
        while A_NQ(n, L, M0) <= 2 * n * (math.log2(2 * n + 1) + h):
            n += 1 if n < 2000 else 50
        print(f"  quadratic equations of height <= 2^{h}: excluded at n = {n}")
    print("  (a root of A x^2 + B x + C with max|A|,|B|,|C| <= H has alpha = A Theta in O_K, w = |A| and")
    print("   |sigma(alpha)| <= 2H, so C_sigma <= (2n+1)H; exclusion needs A_NQ(n) > 2n log2((2n+1)H)).")


# ---------------------------------------------------------------- Q4
def section_Q4():
    hdr("Q4  Euler's function (rho;rho)_inf (lambda = 1): the same transcription at rho = 2^10/3^9")
    L, M0 = 10, 3 ** 9
    lead = L - 2 * (2 / 3 - 1 / 6 - C_CYC) * math.log2(M0)
    print(f"  leading n^3 coefficient of the norm margin: {lead:+.4f} bits: NOT reached")
    print(f"  (non-quadraticity needs mu_bar < 1/(1 - 2 C_cyc) = {1 / (1 - 2 * C_CYC):.5f}; irrationality holds for"
          f" mu_bar < {1 / (0.5 - C_CYC):.5f}, round-2 Theorem E)")


if __name__ == "__main__":
    t0 = time.time()
    print("collatz_procgen_20260923_theta2_nq.py" + (" --quick" if QUICK else ""))
    section_Q0()
    section_Q1()
    section_Q2()
    section_Q3()
    section_Q4()
    print(f"\n[nq done in {time.time() - t0:.1f}s]")
