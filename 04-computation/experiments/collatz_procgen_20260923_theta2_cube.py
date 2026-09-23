#!/usr/bin/env python3
"""
collatz_procgen_20260923_theta2_cube.py

HYP-9127 (cubes), round 3: the natural determinant families for X = sum_(k>=0) rho^(k^3), rho = 2^10/3^9,
and a no-go theorem for them (note: collatz_procgen_20260923_theta_round2.md, section 2).

  C1  Parabola lemma. G(u, v) = sum_k u^k v^(k^2) rho^(k^3) satisfies G(u,v) = 1 + u v rho G(u v^2 rho^3, v rho^3).
      The X-tails are the points P(s) = (3s^2, 3s) of the (a, b)-plane of G(rho^a, rho^b); this parabola has
      no nondegenerate parallelogram and meets every line in <= 2 points.
  C2  Leading layers at the two places (2 and xi = M0/2^L) for sum patterns s_ij = x_i + y_j and the
      normalizations psi(s) = s^3 + c s^2 (c = -1..4) and raw tails (psi = 0): the determinant valuation
      equals the minimal matching (cancellation-free) exactly when that matching is unique; the only
      degenerate cases are (c = 0, place xi) and (c = 3, place 2), where the leading layer has rank one.
  C3  The two degenerate cases (Hankel pattern), proved by row-difference peeling:
        c = 0 (round-2 K7):  ord_xi det F = e3(n) = (n-1)(3n^2 - 15n + 19)   (n >= 3)
        c = 3:               v2 det T = L (n-1)(3n^2 - 3n + 1)               (n >= 2)
      checked exactly, with uniqueness of the peeled optimal assignment.
  C4  Certificate margins of the Hankel family (c = 0, 1, 2, 3): exact v2, exact forced divisor
      (xi-order and cyclotomic factors, gcd over (a, b)), permutation bound. Reach per n (largest log2 C
      excluded) against the single-tail reach L B(mu_bar); cyclotomic multiplicities m_d(n).
  C5  v-direction: Z_n = G(1, rho^n) is an exponential sum in n; v2 det(Z_(i+j)) = L n(n-1)^2(5n-4)/12
      exactly (quartic). The line {(0, n)} meets the X-orbit only at X = Z_0.
  C6  Weighted and residue-class tails: mod-2 collapse of sum k^j x^(k^3); F_even(x) = F(x^8),
      F_(0 mod 3)(x) = F(x^27); two-sequence block-Hankel determinants are cancellation-free;
      Hermite-Pade dimension count mu_bar < (r+2)/(r+1).

Usage: python3 collatz_procgen_20260923_theta2_cube.py [--quick]
"""
import itertools
import math
import random
import sys
import time

import flint
import gmpy2
from gmpy2 import mpz
import numpy as np

QUICK = "--quick" in sys.argv
L, M0 = 10, 3 ** 9
MU = math.log2(M0) / L
LOG2M0 = math.log2(M0)
PRIME = 1048573            # < 2^20, for truncated xi-adic series


def hdr(s):
    print()
    print("=" * 78)
    print(s)
    print("=" * 78)
    sys.stdout.flush()


def v2(x):
    x = mpz(x)
    return 10 ** 9 if x == 0 else int(gmpy2.bit_scan1(abs(x)))


def B_single():
    """single-tail reach: max_s L[(s+1)^3 - mu_bar s^3] (bits)."""
    vals = [(L * ((s + 1) ** 3 - MU * s ** 3), s) for s in range(0, 40)]
    return max(vals)


PSI = {
    "raw": lambda s: 0,
    "c=-1": lambda s: s ** 3 - s * s,
    "c=0": lambda s: s ** 3,
    "c=1": lambda s: s ** 3 + s * s,
    "c=2": lambda s: s ** 3 + 2 * s * s,
    "c=3": lambda s: s ** 3 + 3 * s * s,
    "c=4": lambda s: s ** 3 + 4 * s * s,
}


# ---------------------------------------------------------------- C1
def section_C1():
    hdr("C1  parabola lemma: the X-orbit in the (a,b)-plane of G(rho^a, rho^b) = sum_k rho^(k^3 + b k^2 + a k)")
    N = 3000
    mod = mpz(1) << N
    inv = gmpy2.invert(mpz(M0), mod)

    def rp(e):
        return (mpz(1) << (L * e)) * gmpy2.powmod(inv, e, mod) % mod if L * e < N else mpz(0)

    def G(a, b, shift=0):
        t, k = mpz(0), 0
        while L * (k ** 3 + b * k * k + a * k) < N:
            t += rp(k ** 3 + b * k * k + a * k)
            k += 1
        return t % mod
    ok_fe, ok_tail = True, True
    for a in range(0, 4):
        for b in range(0, 4):
            lhs = G(a, b)
            rhs = (1 + rp(a + b + 1) * G(a + 2 * b + 3, b + 3)) % mod
            ok_fe &= (lhs == rhs)
    X = G(0, 0)
    for s in range(0, 6):
        tail = sum(rp(k ** 3) for k in range(s, 40) if L * k ** 3 < N) % mod
        ok_tail &= ((rp(s ** 3) * G(3 * s * s, 3 * s)) % mod == tail)
    print(f"  functional equation G(rho^a, rho^b) = 1 + rho^(a+b+1) G(rho^(a+2b+3), rho^(b+3)) (mod 2^{N}): {ok_fe}")
    print(f"  tails: sum_(k>=s) rho^(k^3) = rho^(s^3) G(rho^(3s^2), rho^(3s)), i.e. the point P(s) = (3s^2, 3s): {ok_tail}")
    # parallelograms and collinearity on the parabola
    S = 60
    pts = [(3 * s * s, 3 * s) for s in range(S)]
    para = 0
    for s1, s2, s3 in itertools.product(range(S), repeat=3):
        s4 = s1 + s2 - s3
        if 0 <= s4 < S and {s1, s2} != {s3, s4}:
            if pts[s1][0] + pts[s2][0] == pts[s3][0] + pts[s4][0]:
                para += 1
    col = 0
    for s1, s2, s3 in itertools.combinations(range(S), 3):
        (a1, b1), (a2, b2), (a3, b3) = pts[s1], pts[s2], pts[s3]
        if (a2 - a1) * (b3 - b1) - (a3 - a1) * (b2 - b1) == 0:
            col += 1
    print(f"  s < {S}: nondegenerate parallelograms P(s1)+P(s2) = P(s3)+P(s4): {para};  collinear triples: {col}")
    print("  (proof: s1+s2 = s3+s4 and s1^2+s2^2 = s3^2+s4^2 force {s1,s2} = {s3,s4}; equivalently")
    print("   (x1+y1)^2 + (x2+y2)^2 - (x1+y2)^2 - (x2+y1)^2 = 2(x1-x2)(y1-y2) != 0 for a sum pattern)")
    print("  v-direction: the points (0, n) lie on the line a = 0, which meets the parabola only at (0, 0) = X.")


# ---------------------------------------------------------------- 2-adic tools
def tails_psi(psi, svals, N, K):
    """2^(L K) rho^(-psi(s)) R_(s+1) mod 2^N, R_(s+1) = sum_(k>s) rho^(k^3), for s in svals."""
    mod = mpz(1) << N
    inv = gmpy2.invert(mpz(M0), mod)
    out = {}
    for s in svals:
        t, m = mpz(0), 1
        while True:
            e = (s + m) ** 3 - psi(s)
            if L * (e + K) >= N:
                break
            assert e + K >= 0
            term = mpz(1) << (L * (e + K))
            term *= gmpy2.powmod(inv, e, mod) if e >= 0 else gmpy2.powmod(mpz(M0), -e, mod)
            t += term
            m += 1
        out[s] = t % mod
    return out


def v2_det_mod(M, N):
    """2-adic valuation of det(M), entries known mod 2^N; minimal-valuation pivoting. None if not determined."""
    mod = mpz(1) << N
    A = [[mpz(x) % mod for x in row] for row in M]
    n = len(A)
    total = 0
    for k in range(n):
        best = None
        for i in range(k, n):
            for j in range(k, n):
                if A[i][j] != 0:
                    vv = v2(A[i][j])
                    if best is None or vv < best[0]:
                        best = (vv, i, j)
        if best is None:
            return None
        vv, bi, bj = best
        A[k], A[bi] = A[bi], A[k]
        for row in A:
            row[k], row[bj] = row[bj], row[k]
        total += vv
        inv = gmpy2.invert(A[k][k] >> vv, mod)
        for i in range(k + 1, n):
            if A[i][k] != 0:
                f = ((A[i][k] >> vv) * inv) % mod
                for j in range(k, n):
                    A[i][j] = (A[i][j] - f * A[k][j]) % mod
    return total


# ---------------------------------------------------------------- xi-adic tools (truncated series mod PRIME)
def ser_mul(a, b, T):
    return np.convolve(a, b)[:T] % PRIME


def ser_inv(u, T):
    w = np.array([pow(int(u[0]), PRIME - 2, PRIME)], dtype=np.int64)
    k = 1
    while k < T:
        k2 = min(2 * k, T)
        uw = np.convolve(u[:k2], w)[:k2] % PRIME
        corr = (-uw) % PRIME
        corr[0] = (corr[0] + 2) % PRIME
        w = np.convolve(w, corr)[:k2] % PRIME
        k = k2
    return w


def ser_order(a):
    nz = np.nonzero(a)[0]
    return int(nz[0]) if len(nz) else len(a)


def xi_order_det(psi, xs, ys, a, b, T, shift):
    """ord_xi det(F_ij), F_ij = a xi^(psi(s)) - b sum_(k<=s) xi^(psi(s)-k^3), s = x_i + y_j, via series mod xi^T
    (every entry multiplied by xi^shift; the shift is removed at the end). None if T is too small."""
    n = len(xs)
    A = []
    for i in range(n):
        row = []
        for j in range(n):
            s = xs[i] + ys[j]
            c = np.zeros(T, dtype=np.int64)
            e = psi(s) + shift
            if e < T:
                c[e] = (c[e] + a) % PRIME
            for k in range(s + 1):
                e = psi(s) - k ** 3 + shift
                assert e >= 0
                if e < T:
                    c[e] = (c[e] - b) % PRIME
            row.append(c)
        A.append(row)
    total = 0
    for k in range(n):
        best = None
        for i in range(k, n):
            for j in range(k, n):
                o = ser_order(A[i][j])
                if o < T and (best is None or o < best[0]):
                    best = (o, i, j)
        if best is None:
            return None
        v, bi, bj = best
        A[k], A[bi] = A[bi], A[k]
        for row in A:
            row[k], row[bj] = row[bj], row[k]
        total += v
        uinv = ser_inv(A[k][k][v:], T - v)
        for i in range(k + 1, n):
            if ser_order(A[i][k]) < T:
                f = ser_mul(A[i][k][v:], uinv, T - v)
                for j in range(k, n):
                    prod = ser_mul(f, A[k][j][v:], T - v)
                    A[i][j][v:] = (A[i][j][v:] - prod) % PRIME
    return total - n * shift


def match_stats(cost):
    n = len(cost)
    best, count, arg = None, 0, None
    for p in itertools.permutations(range(n)):
        v = sum(cost[i][p[i]] for i in range(n))
        if best is None or v < best:
            best, count, arg = v, 1, p
        elif v == best:
            count += 1
    return best, count, arg


# ---------------------------------------------------------------- C2
def section_C2():
    hdr("C2  leading layers at the places 2 and xi: determinant valuation versus minimal matching")
    random.seed(9127)
    ns = [3, 4, 5] if QUICK else [3, 4, 5, 6]
    patterns = []
    for n in ns:
        patterns.append(("hankel", list(range(n)), list(range(n))))
        for _ in range(3):
            xs = sorted(random.sample(range(0, n + 3), n))
            ys = sorted(random.sample(range(0, n + 3), n))
            patterns.append(("random", xs, ys))
    print("  cost at place 2 (units L): c2(s) = (s+1)^3 - psi(s) (level m = 1);  at xi: cM(s) = psi(s) - s^3 (k = s)")
    print("  precision is set just above the minimal matching m0, so 'det valuation = m0' is decided exactly and")
    print("  'cancellation' means det valuation > m0")
    print("  psi   | place | #patterns | min matching unique | det valuation = min matching | cancellation seen")
    for name, psi in PSI.items():
        for place in ("2", "xi"):
            if place == "xi" and name == "raw":
                continue
            uniq_all, eq_all, canc, tot, bad = 0, 0, 0, 0, 0
            for (kind, xs, ys) in patterns:
                n = len(xs)
                svals = sorted({x + y for x in xs for y in ys})
                if place == "2":
                    cost = [[(xs[i] + ys[j] + 1) ** 3 - psi(xs[i] + ys[j]) for j in range(n)] for i in range(n)]
                else:
                    cost = [[psi(xs[i] + ys[j]) - (xs[i] + ys[j]) ** 3 for j in range(n)] for i in range(n)]
                mn, cnt, _ = match_stats(cost)
                if place == "2":
                    K = max(0, max(psi(s) - (s + 1) ** 3 for s in svals))
                    N = L * (mn + n * K + 5)
                    tl = tails_psi(psi, svals, N, K)
                    val = v2_det_mod([[tl[xs[i] + ys[j]] for j in range(n)] for i in range(n)], N)
                    val = None if val is None else val - n * L * K
                    ref = L * mn
                else:
                    shift = max(0, max(s ** 3 - psi(s) for s in svals))
                    T = mn + n * shift + 2
                    val = xi_order_det(psi, xs, ys, 12345, 678, T, shift)
                    ref = mn
                tot += 1
                uniq_all += (cnt == 1)
                eq_all += (val == ref)
                canc += (val is None or val > ref)
                bad += (val is not None and val < ref)
            print(f"  {name:5s} | {place:5s} | {tot:9d} | {uniq_all:19d} | {eq_all:28d} | {canc}" + (f"  (!! {bad} below m0)" if bad else ""))
    print("  => unique minimal matching <=> no cancellation; ties occur only where the leading cost is affine in s:")
    print("     c = 3 at 2 (c2 = 3s + 1) and c = 0 at xi (cM = 0); there the leading layer has rank one.")


# ---------------------------------------------------------------- exact Hankel determinants over Z[xi]
def F_entry(s, a, b, psi):
    P = psi(s)
    c = [0] * (P + 1)
    c[P] += a
    for k in range(s + 1):
        c[P - k ** 3] -= b
    return flint.fmpz_poly(c)


def bareiss(M):
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


PTS = [(0, 1), (1, 0), (1, 1), (2, 1), (1, 2), (3, 1), (1, 3), (5, 2), (2, 5), (7, 3), (3, 7), (11, 4), (4, 11),
       (13, 6), (6, 13), (17, 5), (5, 17), (19, 7)]
_CYC_BY_DEG = {}
for _d in range(1, 4000):
    _CYC_BY_DEG.setdefault(flint.fmpz_poly.cyclotomic(_d).degree(), []).append(_d)
XI = flint.fmpz_poly([0, 1])


def cyc_index(f):
    for d in _CYC_BY_DEG.get(f.degree(), []):
        if flint.fmpz_poly.cyclotomic(d) == f or flint.fmpz_poly.cyclotomic(d) == -f:
            return d
    return None


def forced_divisor(n, psi):
    """gcd over (a, b) of det F(xi; a, b): returns (xi-order, {d: m_d}, other-degree)."""
    g = None
    for (a, b) in PTS[: n + 2]:
        d = bareiss([[F_entry(i + j, a, b, psi) for j in range(n)] for i in range(n)])
        g = d if g is None else g.gcd(d)
    if g == 0:
        return None
    fac = g.factor()
    xo, cyc, other = 0, {}, 0
    for f, e in fac[1]:
        if f == XI:
            xo += e
        else:
            d = cyc_index(f)
            if d is None:
                other += f.degree() * e
            else:
                cyc[d] = cyc.get(d, 0) + e
    return xo, cyc, other


def e3(n):
    return (n - 1) * (3 * n * n - 15 * n + 19)


def V2_c3(n):
    return L * (n - 1) * (3 * n * n - 3 * n + 1)


# ---------------------------------------------------------------- C3
_FD = {}


def fd_cached(n, key):
    if (n, key) not in _FD:
        psi = PSI[key] if key != "c=3'" else (lambda s: (s + 1) ** 3)
        _FD[(n, key)] = forced_divisor(n, psi)
    return _FD[(n, key)]


def section_C3():
    hdr("C3  the two rank-one cases, Hankel pattern, proved by row-difference peeling")
    nmax = 8 if QUICK else 12
    print("  (i) c = 0 (K7): rows i >= 1 replaced by row_i - row_(i-1); xi-order of det F = e3(n) = (n-1)(3n^2-15n+19)")
    print("   n | xi-order of forced divisor | e3(n) | peeled assignment unique")
    for n in range(2, nmax + 1):
        xo = fd_cached(n, "c=0")[0]
        # peeled cost matrix: row 0 -> 0; row 1, col 0 -> 0; rows i >= 1: E'(i+j-1), E'(0)=0, E'(s)=3s^2-3s+1
        if n <= 8:
            Ep = lambda s: 0 if s == 0 else 3 * s * s - 3 * s + 1
            cost = [[0] * n] + [[Ep(i + j - 1) for j in range(n)] for i in range(1, n)]
            mn, cnt, _ = match_stats(cost)
            u = f"{cnt == 1} (min {mn})"
        else:
            u = "(not enumerated)"
        print(f"  {n:2d} | {xo:26d} | {e3(n) if n >= 3 else 0:5d} | {u}")
    print("  (ii) c = 3, psi(s) = (s+1)^3: t'_s = 1 + rho^(3s^2+9s+7) + ..., leading 2-adic layer = all-ones matrix;")
    print("       after row differences v2 det = L (n-1)(3n^2-3n+1)")
    print("   n | v2 det(t'_(i+j)) | L(n-1)(3n^2-3n+1) | min matching (units L, all tied)")
    psi = lambda s: (s + 1) ** 3
    for n in range(2, (7 if QUICK else 9) + 1):
        svals = list(range(2 * n - 1))
        N = L * ((n - 1) * (3 * n * n - 3 * n + 1) + 400)
        tl = tails_psi(psi, svals, N, 0)
        val = v2_det_mod([[tl[i + j] for j in range(n)] for i in range(n)], N)
        cost = [[(i + j + 1) ** 3 - psi(i + j) for j in range(n)] for i in range(n)]
        mn = sum(cost[i][n - 1 - i] for i in range(n))
        print(f"  {n:2d} | {val:16d} | {V2_c3(n):17d} | {mn} (x{math.factorial(n)} permutations)")


# ---------------------------------------------------------------- C4
def log2_cyc_hom(d):
    """log2 |Phi_d^hom(M0, 2^L)| = phi(d) log2 M0 + log2|Phi_d(rho)|."""
    c = [int(x) for x in flint.fmpz_poly.cyclotomic(d).coeffs()]
    r = 2.0 ** L / M0
    val = sum(ci * r ** i for i, ci in enumerate(c))
    return (len(c) - 1) * LOG2M0 + math.log2(abs(val))


def section_C4():
    hdr("C4  certificate margins of the Hankel family psi_c(s) = s^3 + c s^2 (c = 3 taken as (s+1)^3)")
    Bs, s_star = B_single()
    print(f"  single-tail reach: max_s L[(s+1)^3 - mu_bar s^3] = {Bs:.1f} bits (s = {s_star}); a single tail excludes a/b if")
    print("  log2(|a| + (s+1)|b|) < L[(s+1)^3 - mu_bar s^3]")
    print("  margin_n = V2 + (xi-order) log2 M0 + sum_d m_d log2|Phi_d^hom| - (mu_bar-1) L max_pi sum psi(s) - log2 n!;")
    print("  a rational with |a| + (2n-1)|b| <= C is excluded at size n if margin_n > n log2 C: reach_n = margin_n / n")
    nmax = 7 if QUICK else 10
    table = {}
    for key, lab in (("c=0", "c=0 (K7)"), ("c=1", "c=1"), ("c=2", "c=2"), ("c=3'", "c=3")):
        psi = PSI[key] if key != "c=3'" else (lambda s: (s + 1) ** 3)
        print(f"  {lab}:")
        print("     n |      V2 | xi-ord | cyc deg | other |  quartic loss | margin_n (bits) | reach_n (bits)")
        best = (-1e18, None)
        for n in range(2, nmax + 1):
            if key == "c=0":
                V2 = L * n * (3 * n * n - 3 * n + 1)
            elif key == "c=3'":
                V2 = V2_c3(n)
            else:
                c = int(key[2:])
                V2 = L * n * ((3 - c) * (n - 1) ** 2 + 3 * (n - 1) + 1)
            fd = fd_cached(n, key)
            xo, cyc, other = fd
            cyc_bits = sum(m * log2_cyc_hom(d) for d, m in cyc.items())
            cyc_deg = sum(m * flint.fmpz_poly.cyclotomic(d).degree() for d, m in cyc.items())
            loss = (MU - 1) * L * sum(psi(2 * i) for i in range(n))
            margin = V2 + xo * LOG2M0 + cyc_bits - loss - math.lgamma(n + 1) / math.log(2)
            reach = margin / n
            best = max(best, (reach, n))
            print(f"    {n:2d} | {V2:7d} | {xo:6d} | {cyc_deg:7d} | {other:5d} | {loss:13.1f} | {margin:15.1f} | {reach:13.1f}")
        table[lab] = best
        print(f"     best reach {best[0]:.1f} bits at n = {best[1]}  (single tail: {Bs:.1f})")
    return table


def section_C4b():
    hdr("C4b cyclotomic multiplicities m_d(n) of the forced divisor, c = 0 (K7), versus KRVZ's e_1(n)")
    nmax = 9 if QUICK else 12
    print("   n | xi-order | cyc deg | cyc deg / n^3 | needed for a certificate ~ 2(1-1/mu_bar) n^4 | m_1 | e_1(n) | m_d")
    for n in range(3, nmax + 1):
        xo, cyc, other = fd_cached(n, "c=0")
        deg = sum(m * flint.fmpz_poly.cyclotomic(d).degree() for d, m in cyc.items())
        e1 = (n - 1) ** 2 // 3 if n % 3 == 1 else n * (n - 2) // 3
        need = 2 * (1 - 1 / MU) * n ** 4
        md = ", ".join(f"{d}:{m}" for d, m in sorted(cyc.items()))
        print(f"  {n:2d} | {xo:8d} | {deg:7d} | {deg / n ** 3:13.3f} | {need:44.0f} | {cyc.get(1, 0):3d} | {e1:6d} | {md}")
    print("  The Phi_1 multiplicity equals KRVZ's e_1(n) in range; all other m_d are small. The cyclotomic degree")
    print("  grows like ~0.23 n^3 in range, against the ~0.6 n^4 a certificate would need.")


# ---------------------------------------------------------------- C5
def section_C5():
    hdr("C5  v-direction: Z_n = G(1, rho^n) = sum_k rho^(k^3) (rho^(k^2))^n, Hankel determinants in n")
    nmax = 6 if QUICK else 8
    N = L * (nmax * (nmax - 1) ** 2 * (5 * nmax - 4) // 12 + 200)
    mod = mpz(1) << N
    inv = gmpy2.invert(mpz(M0), mod)

    def Z(n):
        t, k = mpz(0), 0
        while L * (k ** 3 + n * k * k) < N:
            e = k ** 3 + n * k * k
            t += (mpz(1) << (L * e)) * gmpy2.powmod(inv, e, mod)
            k += 1
        return t % mod
    Zs = [Z(n) for n in range(2 * nmax - 1)]
    print("   n | v2 det(Z_(i+j)) | L n(n-1)^2(5n-4)/12")
    for n in range(1, nmax + 1):
        val = v2_det_mod([[Zs[i + j] for j in range(n)] for i in range(n)], N)
        law = L * n * (n - 1) ** 2 * (5 * n - 4) // 12
        print(f"  {n:2d} | {val:15d} | {law}")
    print("  Proof: Cauchy-Binet over subsets K = {k_1 < ... < k_n} of the nodes rho^(k^2); the term has valuation")
    print("  L[sum k^3 + 2 sum_(l<l') k_l^2], uniquely minimal at K = {0, ..., n-1}. The valuation is quartic,")
    print("  but the entries Z_n (n >= 1) are values of G off the X-orbit: they are not in Q + Q X, so under")
    print("  'X rational' the determinant is not an integer and gives no contradiction.")


# ---------------------------------------------------------------- C6
def section_C6():
    hdr("C6  weighted and residue-class tails (Hermite-Pade companions)")
    # mod-2 collapse
    W = 4000
    cubes = [k ** 3 for k in range(0, 20) if k ** 3 < W]
    ok = True
    for j in range(1, 6):
        for k in range(len(cubes)):
            ok &= ((k ** j) % 2 == k % 2)
    print(f"  sum_k k^j x^(k^3) = sum_(k odd) x^(k^3) mod 2 for j = 1..5: {ok}  (over F_2 at most two companions survive)")
    print("  F_even(x) = sum_(k even) x^(k^3) = F(x^8); F_(0 mod 3)(x) = F(x^27): residue classes bring in the new")
    print("  numbers X(rho^8), X(rho^27) (the same problem at another rho, same mu_bar).")
    # block Hankel with weights k^0 and k^1 (normalized tails)
    random.seed(7)
    print("  two-sequence block-Hankel (rows: weight k^0 and weight k^1 tails, columns shifted):")
    print("   m | n = 2m | min-matching (bits) unique | v2 det")
    for m in range(1, (3 if QUICK else 4) + 1):
        n = 2 * m
        svals = list(range(0, 2 * n))
        N = L * (n * ((2 * n + 3) ** 3) + 100)
        mod = mpz(1) << N
        inv = gmpy2.invert(mpz(M0), mod)

        def tw(s, j):
            t, mm = mpz(0), 1
            while L * (mm ** 3 + 3 * mm * mm * s + 3 * mm * s * s) < N:
                e = mm ** 3 + 3 * mm * mm * s + 3 * mm * s * s
                t += (s + mm) ** j * (mpz(1) << (L * e)) * gmpy2.powmod(inv, e, mod)
                mm += 1
            return t % mod
        rows = [(0, i) for i in range(m)] + [(1, i) for i in range(m)]
        Mx = [[tw(r + c, j) for c in range(n)] for (j, r) in rows]
        cost = [[L * (1 + 3 * (r + c) + 3 * (r + c) ** 2) + j * v2(r + c + 1) for c in range(n)] for (j, r) in rows]
        mn, cnt, _ = match_stats(cost)
        val = v2_det_mod(Mx, N)
        print(f"  {m:2d} | {n:6d} | {mn} ({'unique' if cnt == 1 else f'{cnt} ties'}) | {val}")
    print("  dimension count (Siegel's lemma, type-I forms in r+1 companions with the common cube support):")
    for r in range(0, 4):
        thr = (r + 2) / (r + 1)
        print(f"    r = {r}: needs mu_bar < {thr:.4f}  -> {'within' if MU < thr else 'OUT OF'} reach at mu_bar = {MU:.4f}")


if __name__ == "__main__":
    t0 = time.time()
    print("collatz_procgen_20260923_theta2_cube.py" + (" --quick" if QUICK else ""))
    print(f"rho = 2^{L}/3^9, mu_bar = {MU:.5f}")
    section_C1()
    section_C2()
    section_C3()
    section_C4()
    section_C4b()
    section_C5()
    section_C6()
    print(f"\n[cube done in {time.time() - t0:.1f}s]")
