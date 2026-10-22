"""Audit S52 (independent code) - part E: numerology_tests.py claims.

E1  H(P_q) recomputed by my own start-vertex subset DP (numba, uint64; H = q * h_0 by vertex-transitivity) for
    q = 3, 7, 11, 19, 23 and, with 14.5 GB, q = 27 (h_0 = H/27 < 2^64, so uint64 is exact).  q = 31 is taken from the
    repo (paley_H31_exact_20261014.out, audited by the S34 audit with Karp inclusion-exclusion, MISTAKE-609).
E2  |Aut(P_q)| by exhaustive backtracking (q = 3..31, including F_27 = F_3[x]/(x^3 - x - 1)).
E3  H/|Aut| (cube pattern) and H/|Aff| mod 8 (HYP-9273 (B)).
E4  Kervaire-Milnor |bP_4m| (Milnor 2011 Thm 5 form), MISTAKE-563 values, |Theta_23| = 48 |bP_24|; Adams im J;
    Eisenstein coefficients; 691.
E5  SL2(Z)^ab, Mp2(Z)^ab by my own Smith normal form; eta multiplier e^(2 pi i/24) numerically; j coefficients, 744;
    Heegner j-values; e^(pi sqrt 163).
E6  Projective-plane sizes: the strong inventory below 2100 is regenerated with the S50 AUDIT tools (alib, grow; the
    same tools the session used - not independent; the count 723,544 was independently reproduced by S51's audit c4),
    then MY analysis of the session's "within a factor 2 of the local median" statement: base rate over all odd d,
    and a regression with the HYP-9359 arithmetic covariates.
"""
import os, sys, math, time, json, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
os.environ['NUMBA_CACHE_DIR'] = os.path.join(HERE, '_numba_cache')     # keep numba caches inside audit_S52
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
from collections import Counter
import numpy as np
import numba as nb
import mpmath as mp

FAIL = []
T0 = time.time()


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


# ---------------------------------------------------------------- finite fields and Paley tournaments
def gf27():
    # elements a0 + 3 a1 + 9 a2 <-> a0 + a1 x + a2 x^2, x^3 = x + 1
    def dec(e):
        return [e % 3, (e // 3) % 3, e // 9]

    def enc(c):
        return c[0] + 3 * c[1] + 9 * c[2]

    def mul(a, b):
        A, B = dec(a), dec(b)
        prod = [0] * 5
        for i in range(3):
            for j in range(3):
                prod[i + j] += A[i] * B[j]
        # reduce x^4 = x^2 + x, x^3 = x + 1
        prod[2] += prod[4]; prod[1] += prod[4]; prod[4] = 0
        prod[1] += prod[3]; prod[0] += prod[3]; prod[3] = 0
        return enc([prod[0] % 3, prod[1] % 3, prod[2] % 3])

    def sub(a, b):
        A, B = dec(a), dec(b)
        return enc([(x - y) % 3 for x, y in zip(A, B)])
    return mul, sub


def paley(q):
    if q == 27:
        mul, sub = gf27()
        sq = {mul(y, y) for y in range(1, 27)}
        adj = [[1 if a != b and sub(b, a) in sq else 0 for b in range(27)] for a in range(27)]
    else:
        sq = {y * y % q for y in range(1, q)}
        adj = [[1 if (b - a) % q in sq else 0 for b in range(q)] for a in range(q)]
    assert all(adj[a][b] + adj[b][a] == 1 for a in range(q) for b in range(q) if a != b)
    return adj


DEBRUIJN = 0x077CB531
TAB = np.zeros(32, dtype=np.int64)
for _i in range(32):
    TAB[((1 << _i) * DEBRUIJN & 0xFFFFFFFF) >> 27] = _i


@nb.njit(cache=True)
def h0_dp(q, out_others, tab):
    """Number of Hamiltonian paths starting at vertex 0.  out_others[v] = out-neighbours of v among 1..q-1 as bits
    0..q-2 (bit w-1 for vertex w).  dp[mask, v], mask over vertices 1..q-1, uint64."""
    M = 1 << (q - 1)
    dp = np.zeros((M, q), dtype=np.uint64)
    dp[0, 0] = 1
    for mask in range(M):
        for v in range(q):
            c = dp[mask, v]
            if c == 0:
                continue
            free = out_others[v] & ~mask
            while free:
                low = free & (-free)
                w = tab[((low * 0x077CB531) & 0xFFFFFFFF) >> 27]
                dp[mask | low, w + 1] += c
                free ^= low
    tot = np.uint64(0)
    for v in range(q):
        tot += dp[M - 1, v]
    return tot


def H_paley(q):
    adj = paley(q)
    out = np.zeros(q, dtype=np.int64)
    for v in range(q):
        m = 0
        for w in range(1, q):
            if adj[v][w]:
                m |= 1 << (w - 1)
        out[v] = m
    h0 = int(h0_dp(q, out, TAB))
    return q * h0


def aut_order(adj):
    """|Aut| of a vertex-transitive tournament = n * |Stab(0)|; Stab(0) by backtracking."""
    n = len(adj)
    order = [0] + [v for v in range(1, n)]
    count = 0
    img = [-1] * n
    used = [False] * n

    def rec(k):
        nonlocal count
        if k == n:
            count += 1
            return
        v = order[k]
        for c in range(n):
            if used[c]:
                continue
            ok = True
            for j in range(k):
                u = order[j]
                if adj[u][v] != adj[img[u]][c]:
                    ok = False
                    break
            if ok:
                img[v] = c
                used[c] = True
                rec(k + 1)
                used[c] = False
                img[v] = -1

    img[0] = 0
    used[0] = True
    rec(1)
    return n * count


def snf_2x2(M):
    """invariant factors of an integer 2x2 matrix: d1 = gcd of entries, d1 d2 = |det|."""
    g = 0
    for row in M:
        for x in row:
            g = math.gcd(g, abs(x))
    det = abs(M[0][0] * M[1][1] - M[0][1] * M[1][0])
    return [g, det // g if g else 0]


def bern(m):
    B = [Fr(1)]
    for k in range(1, m + 1):
        B.append(-sum(math.comb(k + 1, j) * B[j] for j in range(k)) / (k + 1))
    return B


def main():
    print('E1  Paley H by a start-vertex subset DP (own code)')
    known = {3: 3, 7: 189, 11: 95095, 19: 1172695746915, 23: 15760206976379349,
             27: 418425945688463432031, 31: 19876729650849384218015865}     # repo values (session numerology_tests.py)
    mine = {}
    for q in (3, 7, 11, 19, 23, 27):
        t = time.time()
        mine[q] = H_paley(q)
        check(mine[q] == known[q], f'H(P_{q}) = {mine[q]} (own DP, {time.time() - t:.0f}s) = repo value')
    print('     H(P_31) not recomputed (2^30 x 31 states); repo value from paley_H31_exact_20261014.out, independently reproduced by the S34 audit (Karp inclusion-exclusion)')
    print('E2  |Aut(P_q)| by backtracking')
    aut = {}
    for q in (3, 7, 11, 19, 23, 27, 31):
        t = time.time()
        aut[q] = aut_order(paley(q))
        f = 3 if q == 27 else 1
        check(aut[q] == q * (q - 1) // 2 * f, f'|Aut(P_{q})| = {aut[q]} = q(q-1)/2 * {f} ({time.time() - t:.1f}s)')
    print('E3  cube pattern and HYP-9273')
    rows = []
    for q in (3, 7, 11, 19, 23, 27, 31):
        Hq = known[q]
        aff = q * (q - 1) // 2
        assert Hq % aut[q] == 0 and Hq % aff == 0
        t, ta = Hq // aut[q], Hq // aff
        c = (t - 1) // 8 if (t - 1) % 8 == 0 else None
        cube = False
        if c is not None:
            lo, hi = 0, 1
            while hi ** 3 < c:
                hi *= 2
            while lo < hi:
                mid = (lo + hi) // 2
                if mid ** 3 < c:
                    lo = mid + 1
                else:
                    hi = mid
            cube = lo ** 3 == c
        rows.append((q, t, ta % 8, t % 8, cube))
        print(f'     q = {q:2d}: H/|Aut| = {t}, (H/|Aut| - 1)/8 a cube: {cube}; H/|Aff| mod 8 = {ta % 8}; H/|Aut| mod 8 = {t % 8}')
    check([r[4] for r in rows] == [True, True, True, False, False, False, False], 'H/|Aut| - 1 = 8 m^3 exactly at q = 3, 7, 11 (m = 0, 1, 6)')
    check(all(r[2] == 1 for r in rows) and rows[5][3] == 3,
          'H/|Aff| = 1 mod 8 at all seven q; at q = 27 H/|Aut| = 3 mod 8 (as already recorded in HYP-9273 (B), S33/S34)')
    print('E4  Kervaire-Milnor, Adams, Eisenstein')
    B = bern(40)
    bp = {}
    for m in range(2, 16):
        Bm = abs(B[2 * m])
        bp[m] = 2 ** (2 * m - 2) * (2 ** (2 * m - 1) - 1) * (4 * Bm / m).numerator
    check([bp[m] for m in (2, 3, 4, 5, 6, 7)] == [28, 992, 8128, 261632, 1448424448, 67100672],
          f'|bP_4m| = 2^(2m-2)(2^(2m-1)-1) num(4|B_2m|/m), m = 2..7: {[bp[m] for m in range(2, 8)]} (MISTAKE-563 corrected values)')
    # the alternative a_m-form (a_m = 1 for m even, 2 for m odd) with num(|B_2m|/(4m)) gives the same numbers;
    # the doubled legacy form a_m * num(4|B_2m|/m) does not (MISTAKE-563 (1))
    alt = {m: (1 if m % 2 == 0 else 2) * 2 ** (2 * m - 2) * (2 ** (2 * m - 1) - 1) * (abs(B[2 * m]) / (4 * m)).numerator for m in range(2, 16)}
    legacy = {m: (1 if m % 2 == 0 else 2) * bp[m] for m in range(2, 16)}
    check(alt == bp and legacy[3] == 1984, 'equivalent form a_m 2^(2m-2)(2^(2m-1)-1) num(|B_2m|/4m) agrees for m <= 15; the legacy double-count gives 1984 at m = 3')
    check(48 * bp[6] == 69524373504, '|Theta_23| = |coker J_23| * |bP_24| = 48 * 1448424448 = 69524373504 (MISTAKE-563 reference)')
    primes_m = [m for m in range(2, 16) if all((2 ** (2 * m - 1) - 1) % d for d in range(2, math.isqrt(2 ** (2 * m - 1) - 1) + 1))]
    check(primes_m == [2, 3, 4, 7, 9, 10], f'2^(2m-1) - 1 prime for m in {primes_m} (perfect Euclid factors)')
    imj = [(B[2 * s] / (4 * s)).denominator for s in range(1, 7)]
    eis = [-Fr(4 * s) / B[2 * s] for s in range(1, 7)]
    check(imj == [24, 240, 504, 480, 264, 65520] and [abs(e.numerator) for e in eis] == imj,
          f'Adams |im J| in stem 4s-1 = denom(B_2s/4s) = {imj}; Eisenstein coefficients -4s/B_2s = {[str(e) for e in eis]}: numerators equal the im J orders TAUTOLOGICALLY (one is the reciprocal of the other)')
    check(abs((B[12] / 24).numerator) == 691 and bp[6] % 691 == 0, '691 = |num(B_12/24)| divides |bP_24|')
    print('E5  abelianisations, eta, j, Heegner')
    check(snf_2x2([[1, 3], [4, 0]]) == [1, 12] and snf_2x2([[1, 3], [8, 0]]) == [1, 24],
          'SL2(Z) = <S,T | S^4, S^2 = (ST)^3>, Mp2(Z) = <S,T | S^8, S^2 = (ST)^3>: abelianisations Z/12, Z/24')
    mp.mp.dps = 40
    tau0 = mp.mpc('0.123', '0.917')
    eta = lambda t: mp.exp(2j * mp.pi * t / 24) * mp.qp(mp.exp(2j * mp.pi * t))
    mult = eta(tau0 + 1) / eta(tau0)
    check(abs(mult - mp.exp(2j * mp.pi / 24)) < mp.mpf(10) ** -30, 'eta(tau + 1) = e^(2 pi i/24) eta(tau): the multiplier has order 24')
    N = 6
    s3 = [0] + [sum(d ** 3 for d in range(1, n + 1) if n % d == 0) for n in range(1, N + 2)]
    E4 = [1] + [240 * s3[n] for n in range(1, N + 2)]
    E43 = [0] * (N + 2)
    for i in range(N + 2):
        for j in range(N + 2 - i):
            for k in range(N + 2 - i - j):
                E43[i + j + k] += E4[i] * E4[j] * E4[k]
    inv = [1] + [0] * (N + 1)                  # prod (1 - q^n)^(-24)
    for n in range(1, N + 2):
        for _ in range(24):
            for i in range(n, N + 2):
                inv[i] += inv[i - n]
    qj = [sum(E43[i] * inv[k - i] for i in range(k + 1)) for k in range(N + 2)]
    check(qj[:4] == [1, 744, 196884, 21493760] and 744 == 720 + 24 == 3 * 248, f'q j(q) = {qj[:4]}; 744 = 720 + 24 = 3 * 248')
    heeg = {3: 0, 7: -15, 11: -32, 19: -96, 43: -960, 67: -5280, 163: -640320}
    ok = all(abs(mp.kleinj((1 + mp.sqrt(-d)) / 2) * 1728 - c ** 3) < mp.mpf(10) ** -20 * max(1, abs(c) ** 3) for d, c in heeg.items())
    near = mp.exp(mp.pi * mp.sqrt(163)) - 640320 ** 3 - 744
    check(ok and abs(near + mp.mpf('7.4993e-13')) < mp.mpf('1e-16'), f'Heegner j-values are the listed cubes; e^(pi sqrt 163) - 640320^3 - 744 = {mp.nstr(near, 6)}')

    print('E6  projective-plane sizes vs the strong spectrum below 2100')
    sys.path.insert(0, os.path.join(HERE, '..', '..', 'h_spectrum_multiplicities_20261021', 'audit_S50'))
    import alib as L
    from b3_ear_census_625 import grow
    from multiprocessing import Pool
    LIM = 2100
    inv_ = {3: {L.canon(L.chain((2,))): 3}}
    with Pool(12) as pool:
        n = 4
        while True:
            par = [c for c, h in inv_[n - 1].items() if h <= LIM]
            got, _ = grow(pool, n, par, cap=LIM) if par else ({}, None)
            got = {c: h for c, h in got.items() if h <= LIM}
            if not got:
                break
            inv_[n] = got
            n += 1
    S = Counter(h for dd in inv_.values() for h in dd.values())
    json.dump({str(k): v for k, v in sorted(S.items())}, open(os.path.join(HERE, 'e_strong_spectrum_S.json'), 'w'))
    check(sum(S.values()) == 723544, f'inventory: {sum(S.values())} strong classes with H <= 2100 [{time.time() - T0:.0f}s]')
    odd = [d for d in range(9, 1901, 2)]

    def ratio(d):
        win = sorted(S[x] for x in range(int(d * 0.9) | 1, int(d * 1.1) + 1, 2) if x != d)
        med = win[len(win) // 2] if win else 0
        return (S[d] / med if med else float('nan')), med, len(win)
    pp = [q for q in range(2, 46) if len({p for p in range(2, q + 1) if q % p == 0 and all(p % r for r in range(2, p))}) == 1]
    pg = [q * q + q + 1 for q in pp if q * q + q + 1 < 1901]
    base = [ratio(d)[0] for d in odd if d not in (21,)]
    within = sum(1 for r in base if 0.5 <= r <= 2) / len(base)
    small = [d for d in odd if ratio(d)[2] <= 3]
    print(f'     base rate: {within:.3f} of all odd d in [9, 1899] (except 21) lie within a factor 2 of their local median; windows with <= 3 values: d <= {max(small)}')
    pgr = [(d, S[d], ratio(d)[1], round(ratio(d)[0], 3)) for d in pg if d not in (7, 21)]
    print(f'     PG sizes (d, S, median, ratio): {pgr}')
    check(within > 0.9, f'the factor-2 band contains {100 * within:.1f}% of ALL odd d: the PG statement is consistent but has little power')
    # regression with the HYP-9359 covariates on d in [101, 2099], PG sizes excluded from the fit
    D = np.array([d for d in range(101, 2100, 2) if S[d] > 0 and d not in pg])
    y = np.log(np.array([S[d] for d in D], dtype=float))
    def design(ds):
        ds = np.asarray(ds, dtype=float)
        lx = np.log(ds)
        return np.column_stack([np.ones_like(lx), lx, lx ** 2, lx ** 3, (ds % 3 == 0), (ds % 9 == 0), (ds % 5 == 0)])
    X = design(D)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ beta
    sd = res.std()
    pgd = [d for d in pg if d >= 101]
    z = [(d, round(float((np.log(S[d]) - design([d]) @ beta)[0] / sd), 2)) for d in pgd]
    print(f'     regression log S(d) ~ cubic(log d) + [3|d] + [9|d] + [5|d] on odd d in [101, 2099]: coefficients {np.round(beta[4:], 3).tolist()}, residual s.d. {sd:.3f}')
    print(f'     standardised residuals of the PG sizes: {z}')
    check(all(abs(t) < 3 for _, t in z), 'no projective-plane size >= 101 is an outlier (|z| < 3) once the 3|d, 9|d, 5|d support is controlled')
    print()
    print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f'({time.time() - T0:.0f}s)')


if __name__ == '__main__':
    main()
