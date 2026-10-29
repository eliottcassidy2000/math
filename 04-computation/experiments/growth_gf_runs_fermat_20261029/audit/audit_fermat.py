"""Independent audit of THM-4847.

Strong counts are computed from the exponential reciprocal  S(x) = 1 - 1/T(x),  T(x) = sum_m b^C(m,2) x^m/m!
(the sequence construction of the linear condensation), by FLINT power-series inversion modulo p -- not by
Wright's binomial convolution used in fermat_camion.py.  Exact small values are cross-checked by Wright's recurrence.
"""
import sys
from math import comb, factorial
from fractions import Fraction
import flint
import sympy


def strong_exact(mmax, b):
    """Exact strong counts by Wright's recurrence (big integers), for cross-checks."""
    T = [b ** comb(m, 2) for m in range(mmax + 1)]
    s = [0] * (mmax + 1)
    for m in range(1, mmax + 1):
        s[m] = T[m] - sum(comb(m, k) * s[k] * T[m - k] for k in range(1, m))
    return s


def strong_mod_p(r, p, b):
    """s_b(m) mod p for m <= r < p via 1/T(x) mod p (FLINT inverse_series_trunc)."""
    assert r < p
    fact = [1] * (r + 1)
    for i in range(1, r + 1):
        fact[i] = fact[i - 1] * i % p
    coeffs = [pow(b, comb(m, 2), p) * pow(fact[m], p - 2, p) % p for m in range(r + 1)]
    Tp = flint.nmod_poly(coeffs, p)
    U = Tp.inverse_series_trunc(r + 1)
    u = [int(x) for x in U.coeffs()] + [0] * (r + 1)
    return [0] + [(-u[m]) * fact[m] % p for m in range(1, r + 1)]


def legendre(a, p):
    return sympy.jacobi_symbol(a, p)


def main():
    out = []
    # exact small strong counts
    s2 = strong_exact(40, 2)
    s3 = strong_exact(40, 3)
    out.append('s_2(1..8) = %s (A054946); s_3(1..6) = %s' % (s2[1:9], s3[1:7]))
    for p in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        assert strong_mod_p(min(40, p - 1), p, 2)[1:] == [x % p for x in s2[1:min(40, p - 1) + 1]]
    out.append('series inversion mod p agrees with Wright exact values for p <= 37')
    # s(p) = (2/p), s_3(p) = (3/p) at primes
    rows = []
    for p in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        rows.append('p=%d: s2(p)=%d (2/p)=%d; s3(p)=%d (3/p)=%d' % (p, s2[p] % p if p <= 40 else -9, legendre(2, p),
                                                                   s3[p] % p if p <= 40 else -9, legendre(3, p) if p != 3 else 0))
    out.append('strong counts at primes: ' + '; '.join(rows))
    out.append('  => at the Fermat prime 5: s_2(5) = %d = %d mod 5, i.e. -1, NOT +1; at 17: s_2(17) mod 17 = %d'
               % (s2[5], s2[5] % 5, s2[17] % 17))
    # quasi-periodicity laws (THM-4746 (1),(8)) used for the composite F_5
    for b, sym in ((2, 2), (3, 3)):
        sx = strong_exact(3 * 37, b)
        for p in (5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
            if b == 3 and p == 3:
                continue
            e = legendre(sym, p)
            for m in range(1, 2 * p + 2):
                assert (sx[m + p - 1] - e * sx[m]) % p == 0, (b, p, m)
    out.append('quasi-periodicity s_b(m+p-1) = (b/p) s_b(m) mod p holds for b = 2, 3, primes 5..37, m <= 2p+1')
    # Fermat number tests
    F = lambda n: 2 ** (2 ** n) + 1
    blind = []
    pep = []
    for n in range(1, 15):
        N = F(n)
        c = pow(2, (N - 1) // 2, N)
        T = pow(c, N, N)                          # 2^C(N,2) = (2^((N-1)/2))^N
        blind.append((n, c == 1, T == 1, T == c))
        t3 = pow(pow(3, (N - 1) // 2, N), N, N)   # 3^C(N,2)
        pep.append((n, t3 == N - 1))
    out.append('(C1) n: [2^((N-1)/2)==1, 2^C(N,2)==1, T==c] for F_1..F_14: %s' % blind)
    out.append('(C2) 3^C(N,2) == -1 mod F_n for n=1..14: %s (F_1..F_4 prime, F_5..F_14 composite)' % pep)
    # F_5 at its factors
    N = F(5)
    for p in (641, 6700417):
        assert N % p == 0
        r = N % (p - 1) or (p - 1)
        q = (N - r) // (p - 1)
        sv2 = strong_mod_p(r, p, 2)
        sv3 = strong_mod_p(r, p, 3)
        D = (sv2[r] - pow(2, (r - 1) // 2, p)) % p
        e3 = legendre(3, p)
        s3N = pow(e3 % p, q, p) * sv3[r] % p
        # strong tournament residue s(N) mod p and the circulant count c(N) mod p
        e2 = legendre(2, p)
        sN = pow(e2 % p, q, p) * sv2[r] % p
        cN = pow(2, (N - 1) // 2, p)
        out.append('F_5 at p=%d: r=%d (cofactor mod p-1 = %d), q=%d (parity %d); D(r) mod p = %d; s(N) mod p = %d vs '
                   'c(N) mod p = %d; (3/p) = %d; s_3(N) mod p = %d (Pepin-type pass iff %d)'
                   % (p, r, (N // p) % (p - 1), q, q % 2, D, sN, cN, e3, s3N, p - 1))
    # F_6 at its two factors (extension: needs r = 65537 and r = 274177)
    if '--f6' in sys.argv:
        N = F(6)
        for p in (274177, 67280421310721):
            assert N % p == 0
            r = N % (p - 1) or (p - 1)
            q = (N - r) // (p - 1)
            sv2 = strong_mod_p(r, p, 2)
            D = (sv2[r] - pow(2, (r - 1) // 2, p)) % p
            sv3 = strong_mod_p(r, p, 3)
            e3 = legendre(3, p)
            s3N = pow(e3 % p, q, p) * sv3[r] % p
            out.append('F_6 at p=%d: r=%d, q parity %d; D(r) mod p = %d (strong test passes iff 0); (3/p) = %d; '
                       's_3(N) mod p = %d (Pepin-type pass iff %d)' % (p, r, q % 2, D, e3, s3N, p - 1))
    # (B) function-field factorisations
    rows = []
    for q in (3, 5, 7, 11, 13, 17, 19, 23):
        for n in range(1, 6):
            f = flint.nmod_poly([1] + [0] * (2 ** n - 1) + [1], q)
            fac = f.factor()[1]
            degs = sorted(int(g.degree()) for g, e in fac for _ in range(e))
            m = 2 ** (n + 1)
            d = sympy.n_order(q, m)
            assert all(x == d for x in degs) and len(degs) == 2 ** n // d
            rows.append((q, n, len(degs)))
    red = [(q, n, k) for q, n, k in rows if k == 1]
    out.append('(B) X^(2^n)+1 over F_q (q odd prime <= 23, n <= 5): irreducible only at %s; all factor degrees = ord_(2^(n+1))(q)'
               % red)
    # composition with a nonconstant f keeps both factors nonconstant: spot check
    q = 3
    for n in (2, 3):
        Tg = flint.nmod_poly([2, 1, 0, 1], q)          # f = T^3 + T + 2
        g = Tg ** (2 ** n) + 1
        out.append('    F_3, f = T^3+T+2, n=%d: factor degrees of f^(2^n)+1 = %s'
                   % (n, sorted(int(h.degree()) for h, e in g.factor()[1] for _ in range(e))))
    # constructible polygons / phi(n) power of two
    odd_phi2 = [n for n in range(3, 2 ** 17 + 2, 2) if (lambda t: t & (t - 1) == 0)(sympy.totient(n))]
    out.append('odd n in [3, 2^17] with phi(n) a power of 2: %d values, all dividing 2^32-1: %s'
               % (len(odd_phi2), all((2 ** 32 - 1) % n == 0 for n in odd_phi2)))
    print('\n'.join(out))


if __name__ == '__main__':
    main()
