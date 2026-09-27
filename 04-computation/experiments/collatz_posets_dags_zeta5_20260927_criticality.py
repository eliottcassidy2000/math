#!/usr/bin/env python3
"""collatz_posets_dags_zeta5_20260927_criticality.py -- the three Diophantine criticalities of a Collatz orbit
(session collatz-posets-zeta5-20260927, opus, 2026-09-27).

Setting (Syracuse map U(m) = (3m+1)/2^v on odd m, valuations v_l, d_l = v_1 + ... + v_l, S_(L-1) = sum_(k<L) 3^(L-1-k) 2^(d_k),
C_L = prod_(j<L) (1 + 1/(3 m_j)), so that 2^(d_L) m_L = 3^L n + S_(L-1) = 3^L n C_L). For a word of finite length K
(here: the orbit of n down to 1, truncated before the last odd step, as in the S14 shadow-error note) put
eta_l = sum_(k>=l, k<K) 2^(d_k - d_l)/3^(k-l+1) and xi = n + eta_0, so that xi 3^l/2^(d_l) = m_l + eta_l for l < K.
The rational r_L = 2^(d_L) m_L/3^L = n + S_(L-1)/3^L is the orbit's own approximant. Exact identities checked here:

 (6) Apery / Ridout / product criticality (all L < K, exact Fractions):
      |xi - r_L|_oo = eta_L 2^(d_L)/3^L,   |r_L|_2 = 2^(-d_L),   |1/r_L|_3 = 3^(-L)  (3 does not divide m_L for L >= 1),
      |xi - r_L|_oo * |r_L|_2 = eta_L/3^L,   three-place product = eta_L/9^L = eta_L (n C_L)^2 / H_L^2  with H_L = 2^(d_L) m_L,
      exponents (denominator 3^L): mu_L = -log_3|xi - r_L|/L, nu_L = d_L/(L log_2 3):  mu_L + nu_L = 1 - log_3(eta_L)/L,
      nu_L - 1 = log_2(n C_L/m_L)/(L log_2 3)  (so nu_L > 1 iff coefficient descent at L).
 (7) The Apery-form recurrence A_(L+1) = 3 A_L + 2^(v_2(A_L)), A_0 = n: oddpart(A_L) = m_L and v_2(A_L) = d_L (all odd n < 2*10^5,
      whole orbit), and 2^(d_L) m_L = 3^L n C_L exactly.
 (8) Borel-Dwork radii of F_w(z) = sum_k 2^(d_k) z^k: R_oo = 2^(-limsup d_k/k), R_2 = 2^(liminf d_k/k), so R_oo R_2 =
      2^(liminf - limsup) <= 1; running means d_k/k for the Sturmian valuation word d_k = ceil(k log_2 3), the negative-cycle
      words (1 2)^oo and (1 1 1 2 1 1 4)^oo, and random 2-adic points (mean valuation 2).
Usage: python3 collatz_posets_dags_zeta5_20260927_criticality.py
"""
import math, random
from fractions import Fraction

LOG23 = math.log2(3)


def U(m):
    m = 3 * m + 1; v = 0
    while m % 2 == 0:
        m //= 2; v += 1
    return m, v


def v2(a):
    v = 0
    while a % 2 == 0:
        a //= 2; v += 1
    return v


def orbit_word(n):
    w = []; m = n; ms = [n]
    while m != 1:
        m, v = U(m); w.append(v); ms.append(m)
    return w, ms


def part6():
    print("== (6) the three criticalities, exact, on truncated orbit words ==")
    for n in (27, 703, 6171):
        w, ms = orbit_word(n)
        K = len(w)
        d = [0]
        for v in w:
            d.append(d[-1] + v)
        etas = [sum(Fraction(2 ** (d[k] - d[l]), 3 ** (k - l + 1)) for k in range(l, K)) for l in range(K)]
        xi = n + etas[0]
        C = [Fraction(1)]
        for j in range(K):
            C.append(C[-1] * (1 + Fraction(1, 3 * ms[j])))
        S = [0]  # S[L] = S_(L-1) = sum_(k<L) 3^(L-1-k) 2^(d_k)
        for L in range(1, K + 1):
            S.append(3 * S[-1] + 2 ** d[L - 1])
        ok = True; rows = []
        for L in range(1, K):
            m = ms[L]
            r = Fraction(2 ** d[L] * m, 3 ** L)
            ok &= (r == n + Fraction(S[L], 3 ** L))
            ok &= (xi * 3 ** L / 2 ** d[L] == m + etas[L])
            err = xi - r
            ok &= (err == etas[L] * Fraction(2 ** d[L], 3 ** L)) and err > 0
            ok &= (m % 2 == 1) and (m % 3 != 0)            # |r_L|_2 = 2^(-d_L), |1/r_L|_3 = 3^(-L)
            ok &= (err * Fraction(1, 2 ** d[L]) == etas[L] / 3 ** L)
            three = err * Fraction(1, 2 ** d[L]) * Fraction(1, 3 ** L)
            ok &= (three == etas[L] / 9 ** L)
            H = max(2 ** d[L] * m, 3 ** L)
            ok &= (H == 2 ** d[L] * m) and (Fraction(2 ** d[L] * m) == 3 ** L * n * C[L])
            ok &= (three * H * H == etas[L] * (n * C[L]) ** 2)
            mu = -math.log(float(err), 3) / L
            nu = d[L] / (L * LOG23)
            ok &= abs(mu + nu - (1 - math.log(float(etas[L]), 3) / L)) < 1e-9
            ok &= abs((nu - 1) - math.log2(float(n * C[L] / m)) / (L * LOG23)) < 1e-9
            if L in (1, 5, 10, 20, 30, 40, K - 1):
                rows.append((L, m, round(nu, 4), round(mu, 4), round(float(etas[L]), 3), round(float(n * C[L]), 3), round(float(three * H * H), 3)))
        print(" n = %d (K = %d odd steps to 1): all exact identities hold for L = 1..%d: %s" % (n, K, K - 1, ok))
        print("   L, m_L, nu_L (2-adic exponent), mu_L (real exponent), eta_L, n C_L, three-place product * H^2 = eta_L (n C_L)^2:")
        for row in rows:
            print("   ", row)
    print(" reading: nu_L > 1 exactly when m_L < n C_L (coefficient descent); mu_L + nu_L = 1 - log_3(eta_L)/L <= 1; the Ridout")
    print("  product sits at exponent exactly 2 with constant eta_L (n C_L)^2 >= 1 (Ridout needs 2 + eps and an algebraic target)")


def part7():
    print("== (7) the Apery-form recurrence A_(L+1) = 3 A_L + 2^(v_2(A_L)) ==")
    ok = True; tot = 0
    for n in range(1, 2 * 10 ** 5, 2):
        A = n; m = n; L = 0; d = 0; C = Fraction(1)
        while m != 1:
            C *= (1 + Fraction(1, 3 * m))
            m, v = U(m); d += v; L += 1
            A = 3 * A + 2 ** v2(A)
            ok &= (v2(A) == d) and (A >> d == m)
            if n < 2000:
                ok &= (Fraction(A) == 3 ** L * n * C)
            tot += 1
    print(" oddpart(A_L) = m_L and v_2(A_L) = d_L along every orbit of odd n < 2*10^5 (%d steps), and A_L = 3^L n C_L exactly for n < 2000: %s" % (tot, ok))
    print(" so d_L = log_2(3^L n C_L) - log_2 m_L: the 2-adic exponent d_L/(L log_2 3) exceeds 1 by exactly log_2(n C_L/m_L)/(L log_2 3)")


def running_means(d):
    return [round(d[k] / k, 4) for k in (10, 100, 1000, 10000) if k < len(d)]


def part8():
    print("== (8) Borel-Dwork radii of F_w(z) = sum 2^(d_k) z^k: R_oo R_2 = 2^(liminf d_k/k - limsup d_k/k) <= 1 ==")
    N = 20001
    d_st = [math.ceil(k * LOG23) for k in range(N)]
    v_st = [d_st[k] - d_st[k - 1] for k in range(1, N)]
    print(" Sturmian valuation word d_k = ceil(k log_2 3): letters %s, running means d_k/k at k = 10, 100, 1000, 10000: %s -> log_2 3 = %.4f; R_oo = 1/3, R_2 = 3, product 1" % (sorted(set(v_st)), running_means(d_st), LOG23))
    for name, per in (("(1 2)^oo (cycle of -5)", (1, 2)), ("(1 1 1 2 1 1 4)^oo (cycle of -17)", (1, 1, 1, 2, 1, 1, 4))):
        A = sum(per); p = len(per)
        print(" %s: mean valuation %d/%d = %.4f < log_2 3 (in E_inf), R_oo = 2^(-%.4f), R_2 = 2^(%.4f), product 1; F_w rational (periodic exponents)" % (name, A, p, A / p, A / p, A / p))
    random.seed(1)
    means = []
    for _ in range(20):
        x = random.getrandbits(4000) | 1
        d = [0]; m = x
        for _ in range(1500):
            m, v = U(m); d.append(d[-1] + v)
        means.append(d[-1] / 1500)
    print(" random 2-adic points (4000-bit odd seeds, 1500 steps): mean valuation d_k/k in [%.3f, %.3f] (Terras: 2), all > log_2 3: coefficient descent, R_oo < 1/3" % (min(means), max(means)))
    print(" verdict: every Collatz series has R_oo R_2 <= 1 -- the Borel-Dwork criterion (R_oo R_p > 1 forces rationality) is never applicable; the boundary is populated by irrational series (Sturmian words: Theorem S)")


def main():
    part6(); part7(); part8()


if __name__ == '__main__':
    main()
