"""Independent audit of THM-4845 (A000568 Gaussian average; meromorphic normalized series).

Algorithms deliberately different from the producer's gf_gaussian.py:
 (1) direct Burnside over ALL permutations of S_n (n <= 8): explicit orbits of <sigma> on ordered and unordered
     pairs, reversing orbits detected explicitly; gives tournaments, graphs, digraphs, oriented graphs.
 (2) the Gaussian moments are evaluated through the explicit divisor realisation
     X_j = sqrt(c) sum_{d|j} sqrt(phi(d)) xi_d (variance of sum m_j X_j = c sum_d phi(d) M_d^2), not via the gcd
     matrix; plus a genuinely numerical check by tensor Gauss-Hermite quadrature in (xi_1, xi_3, xi_5) for n <= 6.
 (3) the PDE (4) as the exponent identity e(m) = e(m-e_p) + (p-1)/2 + sum_j gcd(p,j) (m-e_p)_j (all q at once).
 (4) core kernel K[m,t] by brute force over labelled permutations of [m] (m <= 9) with orbit counting, compared
     with the partition formula; P_1..P_4; exact expansion; remainder ratios.
 (5) graph expansion with fixed-point-free cores (brute force m <= 8), Q_1 and Q_2.
"""
import itertools
import math
from fractions import Fraction
from math import comb, factorial, gcd
import numpy as np


def phi(n):
    return sum(1 for k in range(1, n + 1) if gcd(k, n) == 1)


def cycles_of(perm):
    n = len(perm)
    seen = [False] * n
    cyc = []
    for i in range(n):
        if not seen[i]:
            c = []
            j = i
            while not seen[j]:
                seen[j] = True
                c.append(j)
                j = perm[j]
            cyc.append(c)
    return cyc


def pair_orbits(perm):
    """Return (#unordered orbits, #reversing unordered orbits, #ordered orbits) of <perm> on pairs of distinct points."""
    n = len(perm)
    seen = set()
    n_unord = n_rev = 0
    for a in range(n):
        for b in range(a + 1, n):
            if (a, b) in seen:
                continue
            # orbit of the unordered pair
            orb = set()
            x, y = a, b
            reversing = False
            while True:
                key = (min(x, y), max(x, y))
                if key in orb:
                    break
                orb.add(key)
                x, y = perm[x], perm[y]
            # reversing iff ordered pair (a,b) maps to (b,a) for some power
            x, y = a, b
            while True:
                x, y = perm[x], perm[y]
                if (x, y) == (b, a):
                    reversing = True
                    break
                if (x, y) == (a, b):
                    break
            seen |= orb
            n_unord += 1
            n_rev += reversing
    # ordered orbits: each non-reversing unordered orbit splits into 2 ordered orbits, reversing gives 1
    n_ord = 2 * (n_unord - n_rev) + n_rev
    return n_unord, n_rev, n_ord


def burnside_all(n):
    tour = graph = dig = ori = Fraction(0)
    for perm in itertools.permutations(range(n)):
        u, r, o = pair_orbits(perm)
        graph += 2 ** u
        dig += 2 ** o
        ori += 3 ** (u - r)
        if r == 0:
            tour += 2 ** u
    f = factorial(n)
    return tuple(int(x / f) if (x / f).denominator == 1 else x / f for x in (tour, graph, dig, ori))


OEIS = {
    'A000568': [1, 1, 1, 2, 4, 12, 56, 456, 6880, 191536, 9733056, 903753248, 154108311168],
    'A000088': [1, 1, 2, 4, 11, 34, 156, 1044, 12346, 274668, 12005168],
    'A000273': [1, 1, 3, 16, 218, 9608, 1540944, 882033440, 1793359192848],
    'A001174': [1, 1, 2, 7, 42, 582, 21480, 2142288, 575016219, 415939243032],
}


def partitions(n, allowed):
    allowed = sorted(allowed)

    def rec(rem, i):
        if rem == 0:
            yield {}
            return
        if i < 0:
            return
        p = allowed[i]
        for k in range(rem // p, -1, -1):
            for rest in rec(rem - k * p, i - 1):
                d = dict(rest)
                if k:
                    d[p] = k
                yield d
    yield from rec(n, len(allowed) - 1)


def zee(m):
    z = 1
    for j, k in m.items():
        z *= j ** k * factorial(k)
    return z


def divisor_quadratic(m):
    """sum_d phi(d) M_d^2 with M_d = sum_{d | j} m_j (explicit divisor realisation)."""
    J = max(m) if m else 0
    tot = 0
    for d in range(1, J + 1):
        Md = sum(k for j, k in m.items() if j % d == 0)
        if Md:
            tot += phi(d) * Md * Md
    return tot


def gauss_moment(m, family):
    """E[prod W_j^{m_j}] as an exact rational, via the divisor realisation."""
    Q = divisor_quadratic(m)          # Var(sum m_j X_j) = c * Q
    odd = sum(k for j, k in m.items() if j % 2)
    even = sum(k for j, k in m.items() if j % 2 == 0)
    if family in ('tour', 'graph'):    # c = ln2, W = 2^{-[odd]/2} e^X
        e2 = Fraction(Q, 2) - Fraction(odd, 2)
        assert e2.denominator == 1
        return Fraction(2) ** int(e2)
    if family == 'digraph':            # c = 2 ln2, W = e^X / 2
        return Fraction(2) ** (Q - odd - even)
    if family == 'oriented':           # c = ln3, W = 3^{-1/2} e^X (odd), 3^{-1} e^X (even)
        e3 = Fraction(Q, 2) - Fraction(odd, 2) - even
        assert e3.denominator == 1
        return Fraction(3) ** int(e3)
    raise ValueError


def gauss_count(n, family):
    parts = range(1, n + 1, 2) if family == 'tour' else range(1, n + 1)
    tot = Fraction(0)
    for m in partitions(n, list(parts)):
        tot += gauss_moment(m, family) / zee(m)
    return tot


def davis_orbit_formula(n):
    """Davis via orbit counts of the cycle type (not the gcd quadratic form)."""
    tot = Fraction(0)
    for m in partitions(n, list(range(1, n + 1, 2))):
        items = sorted(m.items())
        e = sum(k * (p - 1) // 2 + p * comb(k, 2) for p, k in items)
        e += sum(gcd(p, q) * kp * kq for (p, kp), (q, kq) in itertools.combinations(items, 2))
        tot += Fraction(2 ** e, zee(m))
    assert tot.denominator == 1
    return int(tot)


def gauss_hermite_check(nmax=6, N=80):
    """Numerical E[Z_n(W)] by tensor Gauss-Hermite in xi_1, xi_3, xi_5 (all that Z_n, n <= 6, needs)."""
    x, w = np.polynomial.hermite_e.hermegauss(N)
    w = w / math.sqrt(2 * math.pi)
    ln2 = math.log(2)
    res = []
    for n in range(nmax + 1):
        val = 0.0
        terms = list(partitions(n, list(range(1, n + 1, 2))))
        # X_1 = s xi1 ; X_3 = s(xi1 + sqrt2 xi3) ; X_5 = s(xi1 + 2 xi5)
        s = math.sqrt(ln2)
        for m in terms:
            c = 1.0 / zee(m)
            m1, m3, m5 = m.get(1, 0), m.get(3, 0), m.get(5, 0)
            a1 = s * (m1 + m3 + m5)      # coefficient of xi1
            a3 = s * math.sqrt(2) * m3    # xi3
            a5 = s * 2 * m5               # xi5
            pref = 2 ** (-(m1 + m3 + m5) / 2)
            E1 = float(np.sum(w * np.exp(a1 * x)))
            E3 = float(np.sum(w * np.exp(a3 * x)))
            E5 = float(np.sum(w * np.exp(a5 * x)))
            val += c * pref * E1 * E3 * E5
        res.append(val)
    return res


def e_tour(m):
    items = sorted(m.items())
    e = sum(k * (p - 1) // 2 + p * comb(k, 2) for p, k in items)
    e += sum(gcd(p, q) * kp * kq for (p, kp), (q, kq) in itertools.combinations(items, 2))
    return e


def pde_exponent_check(D=25):
    cnt = 0
    for n in range(1, D + 1):
        for m in partitions(n, list(range(1, n + 1, 2))):
            for p in list(m):
                mm = dict(m)
                mm[p] -= 1
                if not mm[p]:
                    del mm[p]
                rhs = e_tour(mm) + (p - 1) // 2 + sum(gcd(p, j) * k for j, k in mm.items())
                assert e_tour(m) == rhs, (m, p)
                # coefficient identity with the multiplicity factor: m_p / z_m = (1/p) / z_{m - e_p}
                assert Fraction(m[p], zee(m)) == Fraction(1, p * zee(mm))
                cnt += 1
    return cnt


def brute_core_kernel(mmax, kind):
    """K[m,t] = sum over labelled permutations of [m] that are cores of 2^{#pair orbits}; kind 'tour': odd cycles >=3;
    'graph': fixed-point-free (any lengths >= 2), exponent = #unordered pair orbits."""
    K = {}
    for m in range(0, mmax + 1):
        for perm in itertools.permutations(range(m)):
            cyc = cycles_of(perm)
            lens = [len(c) for c in cyc]
            if kind == 'tour' and any(l % 2 == 0 or l < 3 for l in lens):
                continue
            if kind == 'graph' and any(l < 2 for l in lens):
                continue
            u, r, o = pair_orbits(perm) if m >= 2 else (0, 0, 0)
            if kind == 'tour':
                assert r == 0
            t = len(cyc)
            K[(m, t)] = K.get((m, t), 0) + 2 ** u
    return K


def partition_core_kernel(mmax):
    K = {}
    for m in range(0, mmax + 1):
        for mu in partitions(m, list(range(3, m + 1, 2))):
            t = sum(mu.values())
            K[(m, t)] = K.get((m, t), 0) + Fraction(factorial(m) * 2 ** e_tour(mu), zee(mu))
    return {k: int(v) for k, v in K.items()}


def main():
    out = []
    # (1) direct Burnside
    for n in range(0, 9):
        t, g, d, o = burnside_all(n)
        assert t == OEIS['A000568'][n], (n, t)
        assert g == OEIS['A000088'][n], (n, g)
        assert d == OEIS['A000273'][n], (n, d)
        assert o == OEIS['A001174'][n], (n, o)
        # (2) Gaussian functional through the divisor realisation
        assert gauss_count(n, 'tour') == t
        assert gauss_count(n, 'graph') == g
        assert gauss_count(n, 'digraph') == d
        assert gauss_count(n, 'oriented') == o
    out.append('(1) direct Burnside over S_n (explicit pair orbits) = OEIS A000568, A000088, A000273, A001174 for n <= 8')
    out.append('(2) Gaussian functional via the divisor realisation X_j = sqrt(c) sum_{d|j} sqrt(phi d) xi_d '
               'reproduces all four families for n <= 8')
    for n in range(9, 13):
        assert gauss_count(n, 'tour') == OEIS['A000568'][n]
        if n < len(OEIS['A000088']):
            assert gauss_count(n, 'graph') == OEIS['A000088'][n]
        if n < len(OEIS['A001174']):
            assert gauss_count(n, 'oriented') == OEIS['A001174'][n]
    out.append('    and for n = 9..12 against the OEIS terms held here (tournaments to 12, graphs to 10, oriented to 9)')
    a = [davis_orbit_formula(n) for n in range(0, 61)]
    for n in range(0, 31):
        assert gauss_count(n, 'tour') == a[n]
    out.append('(2b) tournaments: Gaussian functional = Davis orbit formula for n <= 30')
    gh = gauss_hermite_check(6)
    out.append('(2c) Gauss-Hermite numerical E[Z_n(W)], n=0..6: ' + ', '.join('%.9f' % v for v in gh))
    # (3) PDE
    c = pde_exponent_check(25)
    out.append('(3) PDE (4) exponent identity (all q): %d checks, degree <= 25' % c)
    # (4) core kernel brute force
    Kb = brute_core_kernel(9, 'tour')
    Kp = partition_core_kernel(60)
    for key, v in Kb.items():
        assert Kp.get(key, 0) == v, (key, v, Kp.get(key))
    for key, v in Kp.items():
        if key[0] <= 9:
            assert Kb.get(key, 0) == v
    out.append('(4a) K[m,t] brute force over labelled permutations (m <= 9) = partition formula: %s'
               % sorted((k, v) for k, v in Kb.items()))
    # P_k coefficients
    Pk = {}
    for k in range(1, 5):
        coeffs = {}
        for m in range(2 * k + 1, 3 * k + 1):
            t = m - 2 * k
            if (m, t) in Kp:
                coeffs[m] = Fraction(Kp[(m, t)]) * Fraction(2) ** (2 * k * m - comb(m, 2))
        Pk[k] = coeffs
        out.append('    P_%d coefficients of C(n,m): %s' % (k, {m: str(v) for m, v in coeffs.items()}))
    for k in range(1, 8):
        lead = Fraction(factorial(3 * k), 3 ** k * factorial(k)) * 2 ** (3 * k * k + k)
        assert Fraction(Kp[(3 * k, k)]) * Fraction(2) ** (2 * k * 3 * k - comb(3 * k, 2)) == lead
    out.append('    leading coefficient ((3k)!/(3^k k!)) 2^(3k^2+k) at C(n,3k) confirmed for k = 1..7')

    def P(k, n):
        if k == 0:
            return Fraction(1)
        s = Fraction(0)
        for m in range(2 * k + 1, min(3 * k, n) + 1):
            t = m - 2 * k
            if (m, t) in Kp:
                s += Fraction(Kp[(m, t)]) * Fraction(2) ** (2 * k * m - comb(m, 2)) * comb(n, m)
        return s
    for n in range(0, 61):
        f = Fraction(a[n] * factorial(n), 2 ** comb(n, 2))
        s = sum(P(k, n) / Fraction(4) ** (k * n) for k in range(0, n // 2 + 1))
        assert f == s, n
    out.append('(4b) f_n = n! a_n / 2^C(n,2) = sum_k 4^(-kn) P_k(n) exactly for n <= 60 (orbit-formula a_n)')
    rows = []
    for Kc in (0, 1, 2, 3):
        rr = []
        for n in (12, 20, 30, 40, 60):
            f = Fraction(a[n] * factorial(n), 2 ** comb(n, 2))
            rem = f - sum(P(k, n) / Fraction(4) ** (k * n) for k in range(Kc + 1))
            lead = P(Kc + 1, n) / Fraction(4) ** ((Kc + 1) * n)
            rr.append('n=%d:%.10f' % (n, float(rem / lead)))
        rows.append('K=%d ' % Kc + ' '.join(rr))
    out.append('(4c) remainder / next term: ' + ' | '.join(rows))
    # meromorphic F: compare (5) at y = 0.3 with the power series, and evaluate near the pole y = 4
    def F_cores(y, mmax=60):
        tot = 0
        for (m, t), K in Kp.items():
            if m > mmax:
                continue
            k = (m - t) // 2
            kappa = Fraction(K) / Fraction(2) ** comb(m, 2)
            tot += float(kappa) * y ** m / (1 - y / 4 ** k) ** (m + 1)
        return tot
    y = 0.3
    ps = sum(float(Fraction(a[n] * factorial(n), 2 ** comb(n, 2))) * y ** n for n in range(0, 61))
    out.append('(4d) F(0.3): power series %.15f, core sum (5) %.15f' % (ps, F_cores(y)))
    for eps in (1e-2, 1e-3, 1e-4):
        yy = 4 * (1 - eps)
        out.append('     (1 - y/4)^4 F(y) at y = 4(1-%g): %.9f (pole coefficient of the single 3-cycle core: 32)'
                   % (eps, (1 - yy / 4) ** 4 * F_cores(yy)))
    # (5) graphs
    gq = []
    for n in range(0, 31):
        tot = Fraction(0)
        for m in partitions(n, list(range(1, n + 1))):
            items = sorted(m.items())
            e = sum(k * (j // 2) + j * comb(k, 2) for j, k in items)
            e += sum(gcd(i, j) * ki * kj for (i, ki), (j, kj) in itertools.combinations(items, 2))
            tot += Fraction(2 ** e, zee(m))
        assert tot.denominator == 1
        gq.append(int(tot))
    assert gq[:len(OEIS['A000088'])] == OEIS['A000088']
    Kg = brute_core_kernel(8, 'graph')
    # partition kernel for graphs (fixed-point-free) up to 30
    Kgp = {}
    for m in range(0, 31):
        for mu in partitions(m, list(range(2, m + 1))):
            items = sorted(mu.items())
            e = sum(k * (j // 2) + j * comb(k, 2) for j, k in items)
            e += sum(gcd(i, j) * ki * kj for (i, ki), (j, kj) in itertools.combinations(items, 2))
            t = sum(mu.values())
            Kgp[(m, t)] = Kgp.get((m, t), 0) + Fraction(factorial(m) * 2 ** e, zee(mu))
    for key, v in Kg.items():
        assert Kgp[key] == v
    for n in range(0, 31):
        f = Fraction(gq[n] * factorial(n), 2 ** comb(n, 2))
        s = Fraction(0)
        for (m, t), v in Kgp.items():
            if m <= n:
                s += v * Fraction(2) ** (-comb(m, 2) - (m - t) * (n - m)) * comb(n, m)
        assert f == s
    Q1 = {m: str(Kgp[(m, m - 1)] * Fraction(2) ** (m - comb(m, 2))) for m in range(2, 3)}
    Q2 = {m: str(Kgp[(m, m - 2)] * Fraction(2) ** (2 * m - comb(m, 2))) for m in range(3, 5) if (m, m - 2) in Kgp}
    out.append('(5) graphs: n! g_n/2^C(n,2) = sum_d 2^(-dn) Q_d(n) exact for n <= 30; brute-force fpf kernel m <= 8 agrees;'
               ' Q_1 coefficients %s, Q_2 coefficients %s' % (Q1, Q2))
    # convergence bound for tournament cores: e(mu) <= C(m,2)/3
    worst = max(Fraction(e_tour(mu), comb(m, 2)) for m in range(3, 41) for mu in partitions(m, list(range(3, m + 1, 2))))
    out.append('(6) max e(mu)/C(m,2) over tournament cores m <= 40: %s (bound 1/3)' % worst)
    print('\n'.join(out))


if __name__ == '__main__':
    main()
