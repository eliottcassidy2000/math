#!/usr/bin/env python3
"""
collatz_procgen_20260923_theta_hankel.py

2-adic transcription of Bezivin's Hankel-determinant method (Bezivin 1998; Choulet 2001;
Krattenthaler-Rochev-Vaananen-Zudilin [KRVZ], Acta Arith. 136 (2009), arXiv 0812.2921) for the
2-adic theta value  Theta(rho) = sum_{k>=0} rho^(k^2),  rho = 2^L/M0  (M0 odd, mu_bar = log2(M0)/L).

Objects (all exact integers):
  tails      t_n = sum_{m>=1} rho^(m^2 + 2mn)            (= rho^(-n^2)(Theta - S_n), S_n partial sum)
  Hankel     T_n = det(t_(i+j))_(0<=i,j<n)
  integers   N_n = det( G_(i+j) 2^(L (i-j)^2) ),  G_m = a M0^(m^2) - b sum_(k<=m) 2^(L k^2) M0^(m^2-k^2)
             so that  N_n = 2^(L K_n) b^n T_n  if Theta = a/b,  K_n = 4 sum_(i<n) i^2.

Theorem H (note, section 1):  v2(T_n) = L n(n+1)(2n-1)/2 exactly (ultrametric Cauchy-Binet), so
N_n != 0 and 2^(v2 N_n) <= |N_n| <= n! C^n M0^(K_n):  irrational if mu_bar < 7/4 (H1).
With the q-order of KRVZ Prop. 2:  M0^(e'(n)) | N_n, e'(n) = 2 e0(n) - n(n-1):  mu_bar < 28/11 (H2).
With the cyclotomic factors of KRVZ Prop. 4:  prod_l Phi_l^hom(M0^2, 2^(2L))^(e_l(n)) | N_n:
mu_bar < 2.87838... (H3).

Sections:
  K0 constants and thresholds
  K1 exact v2(T_n) versus the formula, several maps (incl. beyond phi, beyond 7/4, beyond 28/11)
  K2 integer identity N_n = 2^(L K_n) b^n T_n (fake rationals) and the archimedean bound
  K3 divisibility: X-order e(n) (tied) >= e'(n); cyclotomic factors divide N_n (random a, b)
  K4 margin functions A1, A2, A3 (exact) for the target maps; n needed to exclude height 2^h
  K5 coverage table for square-swap words under m x + 1
  K6 other quadratic families (triangular, pentagonal, 2k^2+k): v2(T_n) formula
  K7 cubes: the naive Hankel determinant of cube tails (why the method stops)

Usage: python3 collatz_procgen_20260923_theta_hankel.py [--quick]
"""
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
PHI = (1 + 5 ** 0.5) / 2


def hdr(s):
    print()
    print("=" * 78)
    print(s)
    print("=" * 78)
    sys.stdout.flush()


def v2(x):
    x = mpz(x)
    return int(gmpy2.bit_scan1(abs(x))) if x else None


def vp(x, p):
    x = mpz(x)
    if x == 0:
        return None
    c = 0
    while x % p == 0:
        x //= p
        c += 1
    return c


def log2abs(x):
    x = abs(mpz(x))
    b = int(x.bit_length())
    if b <= 1000:
        return math.log2(int(x))
    sh = b - 60
    return math.log2(int(x >> sh)) + sh


# ------------------------------------------------------------------ formulas
def Kn(n):
    return 4 * sum(i * i for i in range(n))


def v2T_formula(L, n):
    return L * n * (n + 1) * (2 * n - 1) // 2


def v2N_formula(L, n):
    return L * n * (2 * n - 1) * (7 * n - 1) // 6


def e0_krvz(n):
    """KRVZ (3.18): q-order of V_n for lambda = 0."""
    if n % 2 == 0:
        return n * (n - 2) * (5 * n - 2) // 24
    return n * (n - 1) * (5 * n - 7) // 24


def eprime(n):
    return max(0, 2 * e0_krvz(n) - n * (n - 1))


def el_krvz(l, n):
    """KRVZ (4.1): exponent of Phi_l(q) in V_n, valid for 1 <= l < n/2."""
    return sum((i + l) // (3 * l) + i // (3 * l) for i in range(n))


_cyc_cache = {}


def cyc_hom(l, X, Y):
    """Phi_l^hom(X, Y) = Y^phi(l) Phi_l(X/Y) (exact integer)."""
    if l not in _cyc_cache:
        x = sympy.symbols("x")
        _cyc_cache[l] = [int(c) for c in sympy.Poly(sympy.cyclotomic_poly(l, x), x).all_coeffs()]
    coeffs = _cyc_cache[l]              # highest degree first
    d = len(coeffs) - 1
    tot = mpz(0)
    for i, c in enumerate(coeffs):
        deg = d - i
        tot += c * mpz(X) ** deg * mpz(Y) ** (d - deg)
    return tot


# ------------------------------------------------------------------ 2-adic tails and T_n
def tails_mod(N, L, M0, count, P=lambda m: m * m, a=2):
    """t_s = sum_{m>=1} rho^(P(m) + a m s) mod 2^N for s < count (quadratic family P, step a)."""
    mod = mpz(1) << N
    inv = gmpy2.invert(mpz(M0), mod)
    out = []
    for s in range(count):
        t = mpz(0)
        m = 1
        while L * (P(m) + a * m * s) < N:
            e = P(m) + a * m * s
            t += (mpz(1) << (L * e)) * gmpy2.powmod(inv, e, mod)
            m += 1
        out.append(t % mod)
    return out


def hankel_v2(tl, n, N):
    H = flint.fmpz_mat([[int(tl[i + j]) for j in range(n)] for i in range(n)])
    d = int(H.det()) % (1 << N)
    return v2(d)


def Nn_integer(n, a, b, L, M0):
    Y = mpz(2) ** L
    X = mpz(M0)

    def G(m):
        return a * X ** (m * m) - b * sum(Y ** (k * k) * X ** (m * m - k * k) for k in range(m + 1))
    Gs = [G(m) for m in range(2 * n - 1)]
    M = [[int(Gs[i + j] * Y ** ((i - j) ** 2)) for j in range(n)] for i in range(n)]
    return mpz(int(flint.fmpz_mat(M).det()))


MAPS = [
    ("3x+1 Y (A=9,L=10)", 10, 3 ** 9),
    ("5x+1 HYP-9131 (A=23,L=33)", 33, 5 ** 23),
    ("5x+1 1^7 0^3", 10, 5 ** 7),
    ("5x+1 1^8 0^2", 10, 5 ** 8),
    ("5x+1 1^10 (A=L)", 1, 5),
    ("7x+1 1^10 0 (A=10,L=11)", 11, 7 ** 10),
    ("7x+1 A=L", 1, 7),
    ("9x+1 1^9 0", 10, 9 ** 9),
]


# ------------------------------------------------------------------ K0
def section_K0():
    hdr("K0  thresholds (mu_bar = log2(M0)/L)")
    C_cyc = 5 / 54 - float(sympy.im(sympy.polylog(2, sympy.exp(2 * sympy.pi * sympy.I / 3)))) / (math.pi ** 2 * math.sqrt(3))
    c3 = (7 / 3) / (11 / 12 - 2 * C_cyc)
    print(f"Tschakaloff-Zudilin (Lemma P): phi = {PHI:.6f}")
    print(f"H1 (Hankel, elementary):       7/4 = {7/4:.6f}")
    print(f"H2 (+ KRVZ Prop. 2 q-order):   28/11 = {28/11:.6f}")
    print(f"H3 (+ KRVZ Prop. 4 cyclotomic): (7/3)/(11/12 - 2 C_cyc) = {c3:.6f},  C_cyc = sum_l e_l phi(l)/n^3 = {C_cyc:.8f}")
    print("dictionary check with the real-case gamma = log|num q|/log|den q|: mu_bar = gamma/(gamma-1)")
    for name, g in (("Tschakaloff/Bundschuh", (3 + 5 ** 0.5) / 2), ("naive Hankel", 7 / 3),
                    ("Bezivin 1998", 28 / 15), ("Choulet 2001", 28 / 17), ("KRVZ 2009 Thm 2", 1.53237645)):
        print(f"   {name:22s} gamma > {g:.8f}  <->  mu_bar < {g/(g-1):.6f}")
    return c3


# ------------------------------------------------------------------ K1
def section_K1():
    hdr("K1  exact v2(T_n) of the Hankel determinant of the 2-adic tails (Cauchy-Binet prediction)")
    nmax = 7 if QUICK else 12
    for name, L, M0 in MAPS:
        mu = math.log2(M0) / L
        ok = True
        vals = []
        for n in range(1, nmax + 1):
            pred = v2T_formula(L, n)
            N = pred + 40 * L + 64
            tl = tails_mod(N, L, M0, 2 * n)
            got = hankel_v2(tl, n, N)
            vals.append(got)
            ok = ok and (got == pred)
        print(f"  {name:28s} mu_bar={mu:.4f}: v2(T_n), n=1..{nmax}: {vals[:6]}...  all == L n(n+1)(2n-1)/2: {ok}")


# ------------------------------------------------------------------ K2
def section_K2():
    hdr("K2  integer identity N_n = 2^(L K_n) b^n T_n(a/b) and the archimedean bound (random fake rationals)")
    random.seed(20260923)
    nmax = 5 if QUICK else 7
    for name, L, M0 in MAPS[:3]:
        ok_id = True
        worst = -1e9
        for trial in range(3):
            a = random.randint(-10 ** 6, 10 ** 6)
            b = random.choice([1, 3, 5, 7, 11, 13]) * random.choice([1, -1])
            mu = Fraction(a, b)
            rho = Fraction(2 ** L, M0)
            for n in range(1, nmax + 1):
                Nint = Nn_integer(n, a, b, L, M0)
                # exact rational Hankel determinant of the fake tails rho^(-m^2)(mu - S_m)
                S = [sum(rho ** (k * k) for k in range(m + 1)) for m in range(2 * n - 1)]
                tl = [rho ** (-(m * m)) * (mu - S[m]) for m in range(2 * n - 1)]
                Hq = flint.fmpq_mat([[flint.fmpq(tl[i + j].numerator, tl[i + j].denominator) for j in range(n)]
                                     for i in range(n)])
                dq = Hq.det()
                dq = Fraction(int(dq.p), int(dq.q))
                rhs = Fraction(2) ** (L * Kn(n)) * Fraction(b) ** n * dq
                ok_id = ok_id and (rhs == Fraction(int(Nint)))
                C = abs(a) + (2 * n - 1) * abs(b)
                bound = math.lgamma(n + 1) / math.log(2) + n * math.log2(C) + Kn(n) * math.log2(M0)
                if Nint != 0 and n >= 2:
                    worst = max(worst, log2abs(Nint) - bound)
        print(f"  {name:28s} identity holds: {ok_id};  max over trials, 2<=n<={nmax}, of"
              f" log2|N_n| - log2(n! C^n M0^K_n) = {worst:.2f} (must be <= 0)")


# ------------------------------------------------------------------ K3
def section_K3():
    hdr("K3  divisibility of N_n: X-order (KRVZ Prop. 2) and cyclotomic factors (KRVZ Prop. 4)")
    random.seed(7)
    nmax = 8 if QUICK else 12
    print("  X-order with X = 10007 (prime), Y = 2, min over 5 random (a,b):  e(n) versus e'(n) = 2e0(n) - n(n-1)")
    row = []
    for n in range(1, nmax + 1):
        ev = []
        for t in range(5):
            a = random.randint(-10 ** 6, 10 ** 6)
            b = random.choice([1, 3, 5, 7, 9]) * random.choice([1, -1])
            ev.append(vp(Nn_integer(n, a, b, 1, 10007), 10007))
        e = min(ev)
        row.append((n, e, eprime(n), e >= eprime(n)))
    print("   " + "  ".join(f"n={n}:{e}>={ep}{'' if ok else '!!'}" for n, e, ep, ok in row))
    print(f"   e(n)/n^3 at n={nmax}: {row[-1][1]/nmax**3:.4f}  (limit 5/12 = {5/12:.4f})")
    print("  actual maps: N_n divisible by M0^(e'(n)) * prod_(l<n/2) Phi_l^hom(M0^2, 2^(2L))^(e_l(n))  (5 random (a,b) each)")
    nmax2 = 7 if QUICK else 10
    for name, L, M0 in MAPS[:3] + [MAPS[5]]:
        allok = True
        for n in range(2, nmax2 + 1):
            D = mpz(M0) ** eprime(n)
            for l in range(1, (n + 1) // 2):
                if 2 * l < n:
                    D *= cyc_hom(l, mpz(M0) ** 2, mpz(2) ** (2 * L)) ** el_krvz(l, n)
            for t in range(5):
                a = random.randint(-10 ** 8, 10 ** 8)
                b = random.choice([1, 3, 5, 7, 11]) * random.choice([1, -1])
                Nint = Nn_integer(n, a, b, L, M0)
                if Nint % D != 0:
                    allok = False
        print(f"   {name:28s} n=2..{nmax2}: divisibility holds for all trials: {allok}")


# ------------------------------------------------------------------ K4
def margins(L, M0, n, level):
    """log2 of (2^v2(N_n) * guaranteed odd divisor) - log2(n! M0^K_n); the certificate is
    'no rational a/b with n log2(|a| + (2n-1)|b|) < margin'."""
    m = v2N_formula(L, n) - Kn(n) * math.log2(M0) - math.lgamma(n + 1) / math.log(2)
    if level >= 2:
        m += eprime(n) * math.log2(M0)
    if level >= 3:
        for l in range(1, n):
            if 2 * l < n:
                m += el_krvz(l, n) * log2abs(cyc_hom(l, mpz(M0) ** 2, mpz(2) ** (2 * L)))
    return m


def section_K4(c3):
    hdr("K4  margin functions (exact): A_k(n) = v2(N_n) + log2(odd divisor, level k) - log2(n! M0^K_n)")
    targets = [MAPS[0], MAPS[1], MAPS[3], MAPS[4], MAPS[5], MAPS[6]]
    ns = [5, 10, 20, 40] if QUICK else [5, 10, 20, 40, 80, 120]
    for name, L, M0 in targets:
        mu = math.log2(M0) / L
        lead1 = (7 / 3) * L - (4 / 3) * math.log2(M0)
        lead2 = (7 / 3) * L - (11 / 12) * math.log2(M0)
        lead3 = (7 / 3) * L - ((7 / 3) / c3) * math.log2(M0)
        print(f"  {name:28s} mu_bar={mu:.5f}: leading n^3 coeff  H1 {lead1:+.3f}  H2 {lead2:+.3f}  H3 {lead3:+.3f}")
        for level in (1, 2, 3):
            vals = [margins(L, M0, n, level) for n in ns]
            print(f"      A{level}(n) for n={ns}: " + " ".join(f"{v:.3g}" for v in vals))
    print("  Reading: a rational a/b = Theta(rho) forces n log2(|a|+(2n-1)|b|) >= A_k(n) for EVERY n;")
    print("  A_k(n) ~ c n^3 with c > 0 exactly when mu_bar < 7/4, 28/11, 2.8784 (k = 1, 2, 3).")
    # explicit n needed to exclude heights 2^h at the HYP-9131 map (level 1: elementary)
    name, L, M0 = MAPS[1]
    print(f"  HYP-9131 map, elementary level: least n with A1(n) > n log2(2n 2^h + 2^h):")
    for h in (100, 10 ** 4, 10 ** 6, 10 ** 9):
        n = 1
        while not (margins(L, M0, n, 1) > n * (math.log2(2 * n) + h + 1)):
            n += 1 if n < 100 else max(1, n // 50)
        print(f"      height 2^{h}: n = {n}")


# ------------------------------------------------------------------ K5
def section_K5(c3):
    hdr("K5  coverage: square-swap block words under m x + 1 (blocks of length L with A ones)")
    print("   m | supercritical A/L > | H1 covers A/L < | H2 covers A/L < | H3 covers A/L < | all blocks covered at level")
    for m in (3, 5, 7, 9, 11, 13, 15, 17):
        lm = math.log2(m)
        crit = 1 / lm
        c = [min(1.0, t / lm) for t in (7 / 4, 28 / 11, c3)]
        lev = next((k + 1 for k, t in enumerate((7 / 4, 28 / 11, c3)) if lm < t), None)
        print(f"  {m:2d} | {crit:.4f}             | {c[0]:.4f}          | {c[1]:.4f}          | {c[2]:.4f}          |"
              f" {('H' + str(lev)) if lev else 'none (partial)'}")


# ------------------------------------------------------------------ K6
def section_K6():
    hdr("K6  other quadratic families P(m) = a m(m-1)/2 + s m: v2(T_n) = L[sum_(m<=n) P(m) + a n(n+1)(n-1)/3]")
    nmax = 6 if QUICK else 9
    fams = [("squares", 2, 1), ("triangular", 1, 1), ("pent k(3k-1)/2", 3, 1), ("pent k(3k+1)/2", 3, 2), ("2k^2+k", 4, 3)]
    L, M0 = 10, 3 ** 9
    for fname, a, s in fams:
        P = lambda m, a=a, s=s: a * m * (m - 1) // 2 + s * m
        ok = True
        for n in range(1, nmax + 1):
            pred = L * (sum(P(m) for m in range(1, n + 1)) + a * n * (n + 1) * (n - 1) // 3)
            N = pred + 40 * L + 64
            tl = tails_mod(N, L, M0, 2 * n, P=P, a=a)
            got = hankel_v2(tl, n, N)
            ok = ok and (got == pred)
        print(f"   {fname:16s} (a={a}, s={s}) under 3x+1 Y-map: formula exact for n <= {nmax}: {ok}")


# ------------------------------------------------------------------ K7
def section_K7():
    hdr("K7  cubes: naive Hankel determinant of the cube tails t_s = sum_m rho^(m^3+3m^2 s+3m s^2)")
    L, M0 = 10, 3 ** 9
    nmax = 5 if QUICK else 7
    mu = math.log2(M0) / L
    print("   n | v2(T_n^cube) | clearing K3_n = 8 sum i^3 | margin v2(T) - (mu_bar-1) L K3_n (Y3 map)")
    for n in range(1, nmax + 1):
        K3 = 8 * sum(i ** 3 for i in range(n))
        N = 10 * (3 * n ** 3 + 12 * n ** 2 + 40) + L * 8 * n ** 3
        mod = mpz(1) << N
        inv = gmpy2.invert(mpz(M0), mod)
        tl = []
        for s in range(2 * n):
            t = mpz(0)
            m = 1
            while L * (m ** 3 + 3 * m * m * s + 3 * m * s * s) < N:
                e = m ** 3 + 3 * m * m * s + 3 * m * s * s
                t += (mpz(1) << (L * e)) * gmpy2.powmod(inv, e, mod)
                m += 1
            tl.append(t % mod)
        got = hankel_v2(tl, n, N)
        law = L * (3 * n ** 3 - 3 * n ** 2 + n)
        print(f"  {n:2d} | {got:12d} (law {law}{'' if got == law else ' MISMATCH'}) | {K3:24d} | {got - (mu - 1) * L * K3:10.1f}")
    print("   exact law (proved, note section 2.3): v2(T_n^cube) = L(3n^3 - 3n^2 + n), attained by the")
    print("   anti-diagonal permutation alone (no Cauchy-Binet cancellation), checked above.")
    print("   v2(T_n^cube) grows like 3 L n^3 while the clearing is 2 n^4 (quartic): the margin is")
    print("   negative from small n on, for every mu_bar > 1. The cube tails are not exponential sums in s")
    print("   (factor (rho^(3m))^(s^2)), so Cauchy-Binet gives no rank gain.")


if __name__ == "__main__":
    t0 = time.time()
    print("collatz_procgen_20260923_theta_hankel.py" + (" --quick" if QUICK else ""))
    c3 = section_K0()
    section_K1()
    section_K2()
    section_K3()
    section_K4(c3)
    section_K5(c3)
    section_K6()
    section_K7()
    print(f"\n[hankel done in {time.time()-t0:.1f}s]")
