#!/usr/bin/env python3
"""AUDIT (independent referee, 2026-10-14) of H(P_27), H(P_31) claimed in
05-knowledge/results/paley_h31_exact_rate_and_apery_wheel_20261014.md.

Genuinely different algorithm from the author's colex Held-Karp DP from the arc 0 -> 1:
  KARP INCLUSION-EXCLUSION (Design B of paley_H31_compute_design_kpc1.md):
     H = sum_{S subset V} (-1)^{n-|S|} W_{n-1}(S),
  W_{n-1}(S) = number of walks with n-1 arcs inside the induced subtournament on S (vector iteration).
  * prime q: Z_q rotation reduction, one canonical representative per orbit (all orbits of size q except
    the empty and the full set, q prime);
  * q = 27: translation reduction by the identity sum_{S != 0} f(S) = q sum_{S containing 0} f(S)/|S| for
    translation-invariant f (weights inverted mod p);
  * arithmetic modulo three ~58-bit primes (NOT the author's 31-bit primes), compile-time constants, CRT on two,
    the third as a check;
  * GF(27) built from a DIFFERENT irreducible polynomial (y^3 + 2y^2 + 1) than the author's (y^3 + 2y + 1).
Validation: q = 7, 11 by brute-force DFS and by an exact int64 Held-Karp from all start vertices;
            q = 19, 23 by the same exact Held-Karp (independent of the repo constants).
Usage: python audit_H_karp_ie.py [qlist]   (default 7,11,19,23,27,31)
"""
import sys, time, math
from fractions import Fraction
import numpy as np
import numba as nb

P1 = 288230376151711717
P2 = 288230376151711687
P3 = 288230376151711681
MODS = (P1, P2, P3)


def is_prime(n):
    if n < 2:
        return False
    i = 2
    while i * i <= n:
        if n % i == 0:
            return False
        i += 1
    return True


# ---------------------------------------------------------------- tournaments
def paley_prime(q):
    sq = {(x * x) % q for x in range(1, q)}
    arc = [[((y - x) % q) in sq for y in range(q)] for x in range(q)]
    return arc


def gf27_tables():
    """GF(27) = F_3[y]/(y^3 + 2y^2 + 1): y^3 = -2y^2 - 1 = y^2 + 2.  element a0 + 3 a1 + 9 a2."""
    def digits(a):
        return [a % 3, (a // 3) % 3, (a // 9) % 3]

    def num(d):
        return d[0] + 3 * d[1] + 9 * d[2]

    # irreducibility: no root in F_3
    for t in range(3):
        assert (t ** 3 + 2 * t * t + 1) % 3 != 0

    def mul(a, b):
        A, B = digits(a), digits(b)
        pr = [0] * 5
        for i in range(3):
            for j in range(3):
                pr[i + j] += A[i] * B[j]
        for k in (4, 3):          # y^k = y^(k-3) * y^3 = y^(k-3) (y^2 + 2) = y^(k-1) + 2 y^(k-3)
            c = pr[k]
            pr[k] = 0
            pr[k - 1] += c
            pr[k - 3] += 2 * c
        return num([x % 3 for x in pr[:3]])

    def sub(a, b):
        A, B = digits(a), digits(b)
        return num([(A[i] - B[i]) % 3 for i in range(3)])
    # field sanity: multiplicative group cyclic of order 26 -> some element has order 26
    def order(g):
        x, k = g, 1
        while x != 1:
            x = mul(x, g)
            k += 1
        return k
    assert max(order(g) for g in range(1, 27)) == 26
    sq = {mul(x, x) for x in range(1, 27)}
    assert len(sq) == 13 and 1 in sq
    minus1 = sub(0, 1)
    assert minus1 not in sq                      # -1 nonsquare: q = 3 mod 4 tournament
    arc = [[sub(y, x) in sq for y in range(27)] for x in range(27)]
    return arc


def check_tournament(arc):
    n = len(arc)
    for x in range(n):
        assert not arc[x][x]
        for y in range(x + 1, n):
            assert arc[x][y] != arc[y][x]


def inmasks(arc):
    n = len(arc)
    m = np.zeros(n, dtype=np.int64)
    for v in range(n):
        for u in range(n):
            if arc[u][v]:
                m[v] |= np.int64(1) << np.int64(u)
    return m


# ---------------------------------------------------------------- Karp kernel
@nb.njit(cache=True)
def walks_mod(S, n, inm, pos, verts, ptr, idx, w1, w2, w3, n1, n2, n3):
    """walk counts with n-1 arcs inside S, modulo P1, P2, P3"""
    k = 0
    x = S
    b = 0
    while x:
        if x & 1:
            verts[k] = b
            pos[b] = k
            k += 1
        x >>= 1
        b += 1
    e = 0
    for i in range(k):
        ptr[i] = e
        m = inm[verts[i]] & S
        bb = 0
        while m:
            if m & 1:
                idx[e] = pos[bb]
                e += 1
            m >>= 1
            bb += 1
    ptr[k] = e
    for i in range(k):
        w1[i] = 1
        w2[i] = 1
        w3[i] = 1
    for step in range(n - 1):
        for i in range(k):
            s1 = 0
            s2 = 0
            s3 = 0
            for t in range(ptr[i], ptr[i + 1]):
                j = idx[t]
                s1 += w1[j]
                s2 += w2[j]
                s3 += w3[j]
            n1[i] = s1 % P1
            n2[i] = s2 % P2
            n3[i] = s3 % P3
        for i in range(k):
            w1[i] = n1[i]
            w2[i] = n2[i]
            w3[i] = n3[i]
    t1 = 0
    t2 = 0
    t3 = 0
    for i in range(k):
        t1 += w1[i]
        t2 += w2[i]
        t3 += w3[i]
    return t1 % P1, t2 % P2, t3 % P3, k


@nb.njit(parallel=True, cache=True)
def karp_rotation(n, inm, nchunks):
    """prime n: sum over canonical rotation representatives S (not empty, not full) of n (-1)^(n-|S|) W(S)."""
    full = (np.int64(1) << np.int64(n)) - 1
    lo_all = np.int64(1)
    hi_all = full                      # S in [1, full): excludes empty and full
    total = hi_all - lo_all
    chunk = (total + nchunks - 1) // nchunks
    part = np.zeros((nchunks, 4), dtype=np.int64)
    for c in nb.prange(nchunks):
        lo = lo_all + c * chunk
        hi = min(hi_all, lo + chunk)
        pos = np.zeros(n, dtype=np.int64)
        verts = np.zeros(n, dtype=np.int64)
        ptr = np.zeros(n + 1, dtype=np.int64)
        idx = np.zeros(n * n, dtype=np.int64)
        w1 = np.zeros(n, dtype=np.int64)
        w2 = np.zeros(n, dtype=np.int64)
        w3 = np.zeros(n, dtype=np.int64)
        n1 = np.zeros(n, dtype=np.int64)
        n2 = np.zeros(n, dtype=np.int64)
        n3 = np.zeros(n, dtype=np.int64)
        a1 = np.int64(0)
        a2 = np.int64(0)
        a3 = np.int64(0)
        cnt = np.int64(0)
        S = lo
        while S < hi:
            # canonical = minimal over all rotations
            x = S
            canon = True
            for r in range(1, n):
                x = ((x << 1) | (x >> (n - 1))) & full
                if x < S:
                    canon = False
                    break
            if canon:
                r1, r2, r3, k = walks_mod(S, n, inm, pos, verts, ptr, idx, w1, w2, w3, n1, n2, n3)
                cnt += 1
                v1 = (n * r1) % P1
                v2 = (n * r2) % P2
                v3 = (n * r3) % P3
                if (n - k) % 2 == 0:
                    a1 = (a1 + v1) % P1
                    a2 = (a2 + v2) % P2
                    a3 = (a3 + v3) % P3
                else:
                    a1 = (a1 + P1 - v1) % P1
                    a2 = (a2 + P2 - v2) % P2
                    a3 = (a3 + P3 - v3) % P3
            S += 1
        part[c, 0] = a1
        part[c, 1] = a2
        part[c, 2] = a3
        part[c, 3] = cnt
    return part


@nb.njit(parallel=True, cache=True)
def karp_contains0(n, inm, winv1, winv2, winv3, nchunks):
    """sum over S containing vertex 0 of (-1)^(n-|S|) W(S) * n/|S| (weights given mod P_i)."""
    total = np.int64(1) << np.int64(n - 1)
    chunk = (total + nchunks - 1) // nchunks
    part = np.zeros((nchunks, 4), dtype=np.int64)
    for c in nb.prange(nchunks):
        lo = c * chunk
        hi = min(total, lo + chunk)
        pos = np.zeros(n, dtype=np.int64)
        verts = np.zeros(n, dtype=np.int64)
        ptr = np.zeros(n + 1, dtype=np.int64)
        idx = np.zeros(n * n, dtype=np.int64)
        w1 = np.zeros(n, dtype=np.int64)
        w2 = np.zeros(n, dtype=np.int64)
        w3 = np.zeros(n, dtype=np.int64)
        n1 = np.zeros(n, dtype=np.int64)
        n2 = np.zeros(n, dtype=np.int64)
        n3 = np.zeros(n, dtype=np.int64)
        a1 = np.int64(0)
        a2 = np.int64(0)
        a3 = np.int64(0)
        cnt = np.int64(0)
        T = lo
        while T < hi:
            S = (T << 1) | 1
            r1, r2, r3, k = walks_mod(S, n, inm, pos, verts, ptr, idx, w1, w2, w3, n1, n2, n3)
            cnt += 1
            # r < 2^58, winv < 2^58: product overflows int64 -> split multiplication
            v1 = mulmod1(r1, winv1[k])
            v2 = mulmod2(r2, winv2[k])
            v3 = mulmod3(r3, winv3[k])
            if (n - k) % 2 == 0:
                a1 = (a1 + v1) % P1
                a2 = (a2 + v2) % P2
                a3 = (a3 + v3) % P3
            else:
                a1 = (a1 + P1 - v1) % P1
                a2 = (a2 + P2 - v2) % P2
                a3 = (a3 + P3 - v3) % P3
            T += 1
        part[c, 0] = a1
        part[c, 1] = a2
        part[c, 2] = a3
        part[c, 3] = cnt
    return part


@nb.njit(cache=True)
def _mulmod(a, b, p):
    # a, b < p < 2^58; schoolbook on 29-bit halves of b
    r = np.int64(0)
    lo = b & ((np.int64(1) << 29) - 1)
    hi = b >> 29
    # a * hi * 2^29 + a * lo ; a*lo < 2^87 too large -> split a as well
    al = a & ((np.int64(1) << 29) - 1)
    ah = a >> 29
    # a*b = (ah 2^29 + al)(hi 2^29 + lo) = ah hi 2^58 + (ah lo + al hi) 2^29 + al lo
    t = (ah * hi) % p                       # < 2^58
    for _ in range(29):
        t = (t * 2) % p
    mid = ((ah * lo) % p + (al * hi) % p) % p
    t = (t + mid) % p
    for _ in range(29):
        t = (t * 2) % p
    t = (t + (al * lo) % p) % p
    return t


@nb.njit(cache=True)
def mulmod1(a, b):
    return _mulmod(a, b, P1)


@nb.njit(cache=True)
def mulmod2(a, b):
    return _mulmod(a, b, P2)


@nb.njit(cache=True)
def mulmod3(a, b):
    return _mulmod(a, b, P3)


# ---------------------------------------------------------------- independent exact methods for small q
@nb.njit(cache=True)
def held_karp_exact(n, outm):
    """exact int64 count of Hamiltonian paths (any start): dp[S][v] = paths covering S ending at v."""
    N = 1 << n
    dp = np.zeros((N, n), dtype=np.int64)
    for v in range(n):
        dp[1 << v, v] = 1
    for S in range(1, N):
        for v in range(n):
            c = dp[S, v]
            if c == 0:
                continue
            m = outm[v] & ~S
            w = 0
            while m:
                if m & 1:
                    dp[S | (1 << w), w] += c
                m >>= 1
                w += 1
    tot = 0
    for v in range(n):
        tot += dp[N - 1, v]
    return tot


def outmasks(arc):
    n = len(arc)
    m = np.zeros(n, dtype=np.int64)
    for u in range(n):
        for v in range(n):
            if arc[u][v]:
                m[u] |= np.int64(1) << np.int64(v)
    return m


def dfs_count(arc):
    n = len(arc)
    out = [[v for v in range(n) if arc[u][v]] for u in range(n)]
    sys.setrecursionlimit(10000)
    total = 0
    used = [False] * n

    def rec(u, depth):
        nonlocal total
        if depth == n:
            total += 1
            return
        for v in out[u]:
            if not used[v]:
                used[v] = True
                rec(v, depth + 1)
                used[v] = False
    for s in range(n):
        used[s] = True
        rec(s, 1)
        used[s] = False
    return total


def crt2(r1, r2, m1, m2):
    t = ((r2 - r1) * pow(m1, -1, m2)) % m2
    return r1 + m1 * t, m1 * m2


def run(q, nchunks=2048):
    if q == 27:
        arc = gf27_tables()
    else:
        assert is_prime(q) and q % 4 == 3
        arc = paley_prime(q)
    check_tournament(arc)
    n = q
    inm = inmasks(arc)
    t0 = time.time()
    if q == 27:
        # translation-invariant f; H = n * sum_{S containing 0} (-1)^(n-|S|) W(S)/|S|  (W(empty) = 0)
        winv = []
        for p in MODS:
            arr = np.zeros(n + 1, dtype=np.int64)
            for k in range(1, n + 1):
                arr[k] = (n * pow(k, -1, p)) % p
            winv.append(arr)
        part = karp_contains0(n, inm, winv[0], winv[1], winv[2], nchunks)
        res = [int(part[:, j].astype(object).sum()) % MODS[j] for j in range(3)]
        nsub = int(part[:, 3].sum())
        method = 'Karp IE over the 2^26 subsets containing 0, weight 27/|S|'
    else:
        part = karp_rotation(n, inm, nchunks)
        res = [int(part[:, j].astype(object).sum()) % MODS[j] for j in range(3)]
        nsub = int(part[:, 3].sum())
        # full set (orbit of size 1): W(V) = n ((n-1)/2)^(n-1) for a regular tournament; also direct
        full = (1 << n) - 1
        W_full_closed = n * ((n - 1) // 2) ** (n - 1)
        res = [(res[j] + W_full_closed) % MODS[j] for j in range(3)]
        method = f'Karp IE, Z_{q} rotation representatives ({nsub} canonical classes, expected (2^q-2)/q = {(2 ** q - 2) // q})'
        assert nsub == (2 ** q - 2) // q
    dt = time.time() - t0
    H12, M12 = crt2(res[0], res[1], MODS[0], MODS[1])
    ok3 = (H12 % MODS[2]) == res[2]
    return H12, M12, ok3, res, method, dt


def main():
    qlist = [int(x) for x in sys.argv[1].split(',')] if len(sys.argv) > 1 else [7, 11, 19, 23, 27, 31]
    for p in MODS:
        assert is_prime(p) if p < 10 ** 12 else True
    from sympy import isprime
    print('moduli', MODS, 'prime:', [isprime(p) for p in MODS])
    claimed = {7: 189, 11: 95095, 19: 1172695746915, 23: 15760206976379349,
               27: 418425945688463432031, 31: 19876729650849384218015865}
    for q in qlist:
        H, M, ok3, res, method, dt = run(q)
        print(f'\nq = {q}: {method}; {dt:.1f} s')
        print(f'   residues {res}; CRT(P1,P2) = {H}; third modulus consistent: {ok3}; H < M/2^20: {H < M >> 20}')
        print(f'   claimed H(P_{q}) = {claimed[q]}; MATCH: {H == claimed[q]}')
        if q <= 23:
            arc = paley_prime(q)
            hk = int(held_karp_exact(q, outmasks(arc)))
            print(f'   independent exact Held-Karp (int64, all starts): {hk}; matches Karp IE: {hk == H}')
            if q <= 11:
                d = dfs_count(arc)
                print(f'   brute-force DFS: {d}; matches: {d == H}')
        aff = q * (q - 1) // 2
        N1, rem = divmod(H, aff)
        R = Fraction(H * 2 ** (q - 1), math.factorial(q))
        print(f'   H/|Aff| = H/{aff} = {N1} (remainder {rem}); N_1 odd: {N1 % 2 == 1}; N_1 mod 8 = {N1 % 8}; N_1 mod 256 = {N1 % 256}')
        print(f'   H odd: {H % 2 == 1};  H mod {2 * aff} = {H % (2 * aff)}')
        print(f'   R(q) = H 2^(q-1)/q! = {R.numerator / R.denominator:.12f}  (exact fraction denominator {R.denominator})')
        if q == 31:
            lo, hi = Fraction(258949, 100000), Fraction(260249, 100000)
            print(f'   HYP-2371 interval [2.58949, 2.60249]: inside = {lo <= R <= hi};  R - 2.59599 = {float(R) - 2.59599:+.7f}')
        if q == 27:
            print(f'   H/1053 = {H // 1053} (remainder {H % 1053}); mod 8 = {(H // 1053) % 8};  H mod 8424 = {H % 8424}')


if __name__ == '__main__':
    main()
