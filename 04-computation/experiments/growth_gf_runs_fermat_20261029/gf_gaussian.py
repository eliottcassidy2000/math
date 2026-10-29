"""A000568 (and A000088) as Gaussian averages of a cycle-index exponential; the meromorphic normalized series.

Checks, all in exact rational arithmetic:
 (1) Davis/Burnside a_n for n <= NMAX, against the first OEIS terms.
 (2) The Gaussian identity: with X a centred Gaussian process on odd j with E[X_p X_q] = ln2 gcd(p,q) and
     W_j = 2^(-1/2) e^(X_j),  a_n = E[Z_n(W)],  n Z_n = sum_{p odd} W_p Z_{n-p}.
     Z_n is computed as a polynomial in the W_p by the recursion, and E[prod W_p^m_p] = 2^((m^T G m - |m|)/2)
     (G the gcd matrix) is applied monomialwise; the realisation X_p = sqrt(ln2) sum_{d|p} sqrt(phi(d)) xi_d
     has covariance ln2 * sum_{d | gcd} phi(d) = ln2 gcd (Smith), checked symbolically.
 (3) The same functional over all j (W_j = 2^(-[j odd]/2) e^(X_j), covariance ln2 gcd on all j) gives the
     graph counts A000088; with covariance 2 ln2 gcd and W_j = e^(X_j)/2 it gives the digraph counts A000273.
 (4) The Cameron-Martin shift reproduces the inherited cycle-insertion PDE (checked coefficientwise).
 (5) f_n = n! a_n / 2^C(n,2) = sum_k 4^(-kn) P_k(n) with explicit polynomials P_k of degree 3k
     (THM-514's core kernel), so F(y) = sum f_n y^n = sum_cores c_mu y^|mu| / (1 - y 4^-k(mu))^(|mu|+1) is
     meromorphic on C with poles at 4^k of order exactly 3k+1; remainders after K poles are O(4^-(K+1)n n^3(K+1)).
"""
import sys
from fractions import Fraction
from math import comb, factorial, gcd
from functools import lru_cache
import itertools

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 40
OEIS_A000568 = [1, 1, 1, 2, 4, 12, 56, 456, 6880, 191536, 9733056, 903753248, 154108311168,
                48542114686912, 28401423719122304, 31021002160355166848]
OEIS_A000088 = [1, 1, 2, 4, 11, 34, 156, 1044, 12346, 274668, 12005168, 1018997864, 165091172592]
OEIS_A000273 = [1, 1, 3, 16, 218, 9608, 1540944, 882033440, 1793359192848]
OEIS_A001174 = [1, 1, 2, 7, 42, 582, 21480, 2142288, 575016219, 415939243032]


def phi(n):
    r, m, p = n, n, 2
    while p * p <= m:
        if m % p == 0:
            while m % p == 0:
                m //= p
            r -= r // p
        p += 1
    if m > 1:
        r -= r // m
    return r


def partitions(n, parts):
    """Multiplicity dicts of partitions of n into the allowed parts (descending)."""
    parts = sorted(parts, reverse=True)

    def rec(rem, i):
        if rem == 0:
            yield {}
            return
        if i == len(parts):
            return
        p = parts[i]
        for k in range(rem // p, -1, -1):
            for rest in rec(rem - k * p, i + 1):
                d = dict(rest)
                if k:
                    d[p] = k
                yield d
    yield from rec(n, 0)


def zee(m):
    z = 1
    for j, k in m.items():
        z *= j ** k * factorial(k)
    return z


def quad(m):
    """m^T G m with G the gcd matrix."""
    items = list(m.items())
    return sum(gcd(a, b) * ka * kb for a, ka in items for b, kb in items)


def davis(n):
    tot = Fraction(0)
    for m in partitions(n, range(1, n + 1, 2)):
        e = (quad(m) - sum(m.values())) // 2
        tot += Fraction(2 ** e, zee(m))
    assert tot.denominator == 1
    return int(tot)


def cycle_index_poly(n, parts):
    """Z_0..Z_n as dicts {monomial(tuple of (p,mult)): coefficient} via n Z_n = sum_p W_p Z_{n-p}."""
    Z = [{(): Fraction(1)}]
    for N in range(1, n + 1):
        acc = {}
        for p in parts:
            if p > N:
                break
            for mono, c in Z[N - p].items():
                d = dict(mono)
                d[p] = d.get(p, 0) + 1
                key = tuple(sorted(d.items()))
                acc[key] = acc.get(key, Fraction(0)) + c / N
        Z.append(acc)
    return Z


def expect(mono_dict, q2_exp):
    """E over the Gaussian field: monomial -> 2^(exponent) with exponent a function of the multiplicity dict."""
    tot = Fraction(0)
    for mono, c in mono_dict.items():
        m = dict(mono)
        e2 = q2_exp(m)            # exponent of 2, may be a half-integer times 2 -> keep as Fraction
        assert e2.denominator == 1, (m, e2)
        tot += c * (Fraction(2) ** int(e2))
    return tot


def main():
    out = []
    # (1) Davis
    a = [davis(n) for n in range(NMAX + 1)]
    assert a[:len(OEIS_A000568)] == OEIS_A000568, 'Davis mismatch'
    out.append('(1) Davis a_n n<=%d; first %d agree with OEIS A000568' % (NMAX, len(OEIS_A000568)))
    # Smith: covariance of the realisation
    for p in range(1, 60):
        for q in range(1, 60):
            assert sum(phi(d) for d in range(1, min(p, q) + 1) if p % d == 0 and q % d == 0) == gcd(p, q)
    out.append('(2a) sum_{d|p,d|q} phi(d) = gcd(p,q) for p,q < 60 (covariance of X_p = sqrt(ln2) sum sqrt(phi d) xi_d)')
    # (2) Gaussian identity for tournaments
    NG = min(NMAX, 30)
    Z = cycle_index_poly(NG, list(range(1, NG + 1, 2)))
    tour_exp = lambda m: Fraction(quad(m) - sum(m.values()), 2)       # W_p = 2^-1/2 e^X_p, Var-form ln2 gcd
    for n in range(NG + 1):
        assert expect(Z[n], tour_exp) == a[n]
    out.append('(2) a_n = E[Z_n(W)] with the odd cycle-index recursion, n<=%d (exact)' % NG)
    # (3) graphs and digraphs from the same functional on all j
    NH = 12
    Zall = cycle_index_poly(NH, list(range(1, NH + 1)))
    graph_exp = lambda m: Fraction(quad(m) - sum(k for j, k in m.items() if j % 2), 2)
    digraph_exp = lambda m: Fraction(quad(m) - sum(m.values()), 1)
    g = [expect(Zall[n], graph_exp) for n in range(NH + 1)]
    dg = [expect(Zall[n], digraph_exp) for n in range(min(NH, len(OEIS_A000273) - 1) + 1)]
    assert g == OEIS_A000088[:NH + 1], g
    assert dg == OEIS_A000273[:len(dg)], dg
    out.append('(3) graphs A000088 n<=%d and digraphs A000273 n<=%d from the same Gaussian functional' % (NH, len(dg) - 1))
    # (3b) oriented graphs: 3 states per non-reversing pair orbit, 1 per reversing (antipodal) orbit of an even cycle;
    #      base-3 field: covariance ln3 gcd, W_j = 3^(-1/2) e^X_j (j odd), 3^(-1) e^X_j (j even)
    def orient_exp(m):
        return Fraction(quad(m) - sum(k for j, k in m.items() if j % 2), 2) - sum(k for j, k in m.items() if j % 2 == 0)
    NO = min(NH, len(OEIS_A001174) - 1)
    og = []
    for n in range(NO + 1):
        tot = Fraction(0)
        for mono, c in Zall[n].items():
            e3 = orient_exp(dict(mono))
            assert e3.denominator == 1
            tot += c * Fraction(3) ** int(e3)
        og.append(tot)
    assert og == OEIS_A001174[:NO + 1], og
    out.append('(3b) oriented graphs A001174 n<=%d from the base-3 field (covariance ln3 gcd)' % NO)
    # (4) Cameron-Martin = cycle-insertion PDE, coefficientwise on F(y;2) = sum_m 2^t(m) prod y^m/(j^m m!):
    #     d/dy_p F = 2^((p-1)/2)/p * F(2^gcd(p,j) y_j)
    NP = 21
    def coeff(m):
        return Fraction(2 ** ((quad(m) - sum(m.values())) // 2), zee(m))
    checked = 0
    for n in range(1, NP + 1):
        for m in partitions(n, range(1, n + 1, 2)):
            for p in list(m):
                # coefficient of y^(m - e_p) in dF/dy_p: m_p * coeff(m)
                mm = dict(m); mm[p] -= 1
                if mm[p] == 0:
                    del mm[p]
                lhs = m[p] * coeff(m)
                rhs = Fraction(2 ** ((p - 1) // 2), p) * coeff(mm) * Fraction(2) ** sum(gcd(p, j) * k for j, k in mm.items())
                assert lhs == rhs
                checked += 1
    out.append('(4) Cameron-Martin shift = inherited PDE: %d coefficient identities, degree <= %d' % (checked, NP))
    # (5) meromorphic normalised series
    # core kernel K[m,t] (THM-514): sum over partitions of m into odd parts >= 3 with t parts of (m!/z) 2^e
    NF = NMAX
    K = {}
    for msum in range(0, NF + 1):
        for mu in partitions(msum, range(3, msum + 1, 2)):
            t = sum(mu.values())
            e = (quad(mu) - t) // 2
            K[(msum, t)] = K.get((msum, t), 0) + Fraction(factorial(msum) * 2 ** e, zee(mu))
    for key, v in K.items():
        assert v.denominator == 1
    def P(k, n):
        if k == 0:
            return Fraction(1)
        tot = Fraction(0)
        for msum in range(2 * k + 1, 3 * k + 1):
            t = msum - 2 * k
            if (msum, t) in K and msum <= n:
                tot += K[(msum, t)] * Fraction(2) ** (-comb(msum, 2) + 2 * k * msum) * comb(n, msum)
        return tot
    for n in range(NF + 1):
        f = Fraction(a[n] * factorial(n), 2 ** comb(n, 2))
        s = sum(P(k, n) * Fraction(1, 4 ** (k * n)) for k in range(0, n // 2 + 1))
        assert f == s, n
    out.append('(5) f_n = n! a_n / 2^C(n,2) = sum_k 4^(-kn) P_k(n) exactly for n <= %d' % NF)
    # explicit low P_k as integer combinations of binomials C(n,m)
    for k in range(1, 5):
        terms = []
        for msum in range(2 * k + 1, 3 * k + 1):
            t = msum - 2 * k
            if (msum, t) in K:
                c = K[(msum, t)] * Fraction(2) ** (-comb(msum, 2) + 2 * k * msum)
                terms.append('%s*C(n,%d)' % (c, msum))
        out.append('    P_%d(n) = %s' % (k, ' + '.join(terms)))
    # remainder decay: f_n - sum_{k<=K} 4^-kn P_k(n) versus 4^-(K+1)n P_{K+1}(n)
    for Kc in (0, 1, 2):
        rows = []
        for n in (20, 30, 40):
            if n > NF:
                continue
            f = Fraction(a[n] * factorial(n), 2 ** comb(n, 2))
            rem = f - sum(P(k, n) * Fraction(1, 4 ** (k * n)) for k in range(Kc + 1))
            lead = P(Kc + 1, n) * Fraction(1, 4 ** ((Kc + 1) * n))
            rows.append('n=%d rem/lead=%.12f' % (n, float(rem / lead)))
        out.append('    after poles k<=%d: %s' % (Kc, '; '.join(rows)))
    # (6) graphs: g_n n!/2^C(n,2) = sum_d 2^(-dn) Q_d(n), cores = fixed-point-free permutations, d = |mu| - l(mu)
    NGr = min(NMAX, 30)
    def graph_count(n):
        tot = Fraction(0)
        for m in partitions(n, range(1, n + 1)):
            e = (quad(m) - sum(k for j, k in m.items() if j % 2)) // 2
            tot += Fraction(2 ** e, zee(m))
        assert tot.denominator == 1
        return int(tot)
    gl = [graph_count(n) for n in range(NGr + 1)]
    assert gl[:len(OEIS_A000088)] == OEIS_A000088
    Kg = {}
    for msum in range(0, NGr + 1):
        for mu in partitions(msum, range(2, msum + 1)):
            t = sum(mu.values())
            e = (quad(mu) - sum(k for j, k in mu.items() if j % 2)) // 2
            Kg[(msum, t)] = Kg.get((msum, t), 0) + Fraction(factorial(msum) * 2 ** e, zee(mu))
    for n in range(NGr + 1):
        f = Fraction(gl[n] * factorial(n), 2 ** comb(n, 2))
        s2 = Fraction(0)
        for (msum, t), v in Kg.items():
            if msum <= n:
                d = msum - t
                s2 += v * Fraction(2) ** (-comb(msum, 2) - d * (n - msum)) * comb(n, msum)
        assert f == s2, n
    q1 = [str(Kg[(m, m - 1)] * Fraction(2) ** (-comb(m, 2) + 1 * m)) + '*C(n,%d)' % m for m in (2,) if (m, m - 1) in Kg]
    out.append('(6) graphs: g_n n!/2^C(n,2) = sum_d 2^(-dn) Q_d(n) exactly for n <= %d (poles at 2^d of order 2d+1); Q_1(n) = %s'
               % (NGr, ' + '.join(q1)))
    print('\n'.join(out))


if __name__ == '__main__':
    main()
