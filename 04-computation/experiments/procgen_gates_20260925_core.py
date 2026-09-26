#!/usr/bin/env python3
"""procgen_gates_20260925_core.py -- shared library for the cycle-gate equidistribution lane.

Session collatz-procgen-20260922, lane "gates" (2026-09-25).

Conventions (THM-4471 section 4; mod-192 note section 2.4):
  * a parity word w of length p with a ones at positions s_0 < ... < s_{a-1};
  * carry  c_w = sum_i q^(a-1-i) 2^(s_i)   (q = 3 for Collatz, q = 5 for the DRIFT control);
  * gate   G = 2^p - q^a;  the unique rational periodic point with word w is x_w = c_w / G
    for the map x -> x/2, (q x + 1)/2.  For G < 0 the point -x_w > 0 is the periodic point of
    x -> x/2, (q x - 1)/2 with the same word (SHEET transport x -> -x).
  * check: 2^p T^p(x) = q^a x + c_w for T(x) = x/2, (qx+1)/2 (one line: y_{t+1} = 3 y_t + 2^t
    at an odd step, y_t = 2^t T^t(x)).
  * exponential sum  S(h) = sum_w e(h c_w / M),  M = |G|,  e(t) = exp(2 pi i t);
    counting identity  #{w : M | c_w} = (1/M) sum_{h mod M} S(h).

Everything here is exact integer arithmetic except the float DP (checked against mpmath and an
exact finite-field DP) and float summaries.
"""
import itertools
import math
import resource
import sys

import numpy as np

TWO_PI = 2.0 * math.pi


# ----------------------------------------------------------------------------------------------
# basic objects
# ----------------------------------------------------------------------------------------------
def gate(p, a, q=3):
    return (1 << p) - q ** a


def carry(pos, a, q=3):
    """c_w for the word with ones at the increasing positions pos (len(pos) == a)."""
    return sum(q ** (a - 1 - i) * (1 << s) for i, s in enumerate(pos))


def carry_word(w, q=3):
    """c_w for a 0/1 tuple w."""
    pos = [t for t, b in enumerate(w) if b]
    return carry(pos, len(pos), q)


def cmin_cmax(p, a, q=3):
    """min and max of c_w over words of length p with a >= 1 ones (ones first / ones last)."""
    lo = (q ** a - 2 ** a) // (q - 2)
    return lo, (1 << (p - a)) * lo


def words(p, a):
    return itertools.combinations(range(p), a)


def rotate_carry(c, first_letter, G, q=3):
    """c_{Rw} from c_w, R = move the first letter to the end.  Exact:
    first letter 0: c/2 ;  first letter 1: (q c + G)/2  (G = 2^p - q^a, signed)."""
    if first_letter == 0:
        assert c % 2 == 0
        return c // 2
    v = q * c + G
    assert v % 2 == 0
    return v // 2


def lyndon_count(p, a):
    """number of primitive necklaces (Lyndon words) of length p with a ones."""
    g = math.gcd(p, a)
    tot = 0
    for d in range(1, g + 1):
        if g % d == 0:
            tot += mobius(d) * math.comb(p // d, a // d)
    assert tot % p == 0
    return tot // p


def mobius(n):
    if n == 1:
        return 1
    res, m, f = 1, n, 2
    while f * f <= m:
        if m % f == 0:
            m //= f
            if m % f == 0:
                return 0
            res = -res
        f += 1
    if m > 1:
        res = -res
    return res


def factor_small(n, bound=10 ** 6):
    """trial division up to bound; returns (list of (prime, exp), cofactor)."""
    n = abs(n)
    out = []
    f = 2
    while f <= bound and f * f <= n:
        if n % f == 0:
            e = 0
            while n % f == 0:
                n //= f
                e += 1
            out.append((f, e))
        f += 1 if f == 2 else 2
    return out, n


# ----------------------------------------------------------------------------------------------
# exponential sums by dynamic programming over positions
# ----------------------------------------------------------------------------------------------
def phase_table(p, a, M, hs, q=3):
    """integer phases (h * q^(a-1-i) * 2^t) mod M, shape (len(hs), a, p), as Python ints -> float
    fraction array theta[h, i, t] in [0,1)."""
    base = [[(pow(q, a - 1 - i, M) * pow(2, t, M)) % M for t in range(p)] for i in range(a)]
    H = len(hs)
    th = np.empty((H, a, p), dtype=np.float64)
    if M < (1 << 31):
        hv = np.array([int(h) % M for h in hs], dtype=np.int64)
        bv = np.array(base, dtype=np.int64)  # (a, p)
        prod = (hv[:, None, None] * bv[None, :, :]) % M
        th[:] = prod / M
    else:
        for k, h in enumerate(hs):
            h = int(h) % M
            for i in range(a):
                for t in range(p):
                    th[k, i, t] = ((h * base[i][t]) % M) / M
    return th


def S_dp(p, a, hs, q=3, normalize=False):
    """S(h) = sum_w e(h c_w / M) for all h in hs (vectorised float DP, O(p a len(hs))).
    If normalize, returns S(h)/C(p,a) computed with rescaling (safe for large p)."""
    M = abs(gate(p, a, q))
    hs = list(hs)
    if a == 0:
        return np.ones(len(hs), dtype=complex)
    th = phase_table(p, a, M, hs, q)
    E = np.exp(1j * TWO_PI * th)  # (H, a, p)
    H = len(hs)
    v = np.zeros((a + 1, H), dtype=complex)
    v[0, :] = 1.0
    logscale = 0.0
    for t in range(p):
        lo = max(0, a - (p - t))  # states that can still finish
        for i in range(min(t, a - 1), lo - 1, -1):
            v[i + 1] += v[i] * E[:, i, t]
        if lo > 0:
            v[:lo] = 0.0
        if normalize:
            s = np.abs(v).max()
            if s > 0:
                v /= s
                logscale += math.log(s)
    if normalize:
        lc = math.lgamma(p + 1) - math.lgamma(a + 1) - math.lgamma(p - a + 1)
        return v[a] * math.exp(logscale - lc)
    return v[a].copy()


def S_dp_mpmath(p, a, h, q=3, dps=40):
    """same DP in mpmath at dps digits (single h) -- precision certificate for S_dp."""
    import mpmath
    mpmath.mp.dps = dps
    M = abs(gate(p, a, q))
    v = [mpmath.mpc(0)] * (a + 1)
    v[0] = mpmath.mpc(1)
    for t in range(p):
        for i in range(min(t, a - 1), -1, -1):
            ph = (h * pow(q, a - 1 - i, M) * pow(2, t, M)) % M
            v[i + 1] += v[i] * mpmath.expjpi(2 * mpmath.mpf(ph) / M)
    return complex(v[a])


def S_brute(p, a, h, q=3):
    M = abs(gate(p, a, q))
    s = 0j
    for pos in words(p, a):
        s += complex(math.cos(TWO_PI * ((h * carry(pos, a, q)) % M) / M),
                     math.sin(TWO_PI * ((h * carry(pos, a, q)) % M) / M))
    return s


# exact arithmetic: the ring map Z[zeta_M] -> F_l for a prime l = 1 mod M
def is_prime(n):
    if n < 2:
        return False
    small = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    for sp in small:
        if n % sp == 0:
            return n == sp
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for b in small:
        x = pow(b, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def prime_1_mod(M, start_k=1):
    k = start_k
    while True:
        ell = k * M + 1
        if is_prime(ell):
            return ell
        k += 1


def root_of_unity(M, ell):
    """a primitive M-th root of unity mod the prime ell (M | ell - 1)."""
    fac, cof = factor_small(M, 10 ** 7)
    primes = [f for f, _ in fac] + ([cof] if cof > 1 else [])
    g = 2
    while True:
        z = pow(g, (ell - 1) // M, ell)
        if all(pow(z, M // r, ell) != 1 for r in primes):
            return z
        g += 1


def S_dp_exact_all(p, a, q=3):
    """Exact image of (S(h))_{h mod M} in F_ell (ell prime, ell = 1 mod M), all h at once.
    Returns (ell, array of S(h) mod ell).  Since sum_h S(h) = M * N exactly in Z[zeta_M] and the
    reduction map is a ring homomorphism, N = (sum_h S(h)) / M mod ell exactly."""
    M = abs(gate(p, a, q))
    ell = prime_1_mod(M, start_k=max(1, (1 << 20) // M))
    assert ell < (1 << 31)
    z = root_of_unity(M, ell)
    zp = np.empty(M, dtype=np.int64)  # zeta^k mod ell
    acc = 1
    for k in range(M):
        zp[k] = acc
        acc = acc * z % ell
    hv = np.arange(M, dtype=np.int64)
    v = np.zeros((a + 1, M), dtype=np.int64)
    v[0, :] = 1
    for t in range(p):
        for i in range(min(t, a - 1), -1, -1):
            m = (pow(q, a - 1 - i, M) * pow(2, t, M)) % M
            idx = (hv * m) % M
            v[i + 1] = (v[i + 1] + v[i] * zp[idx]) % ell
    return ell, v[a]


# ----------------------------------------------------------------------------------------------
# residue lists: prefix / suffix split (meet in the middle)
# ----------------------------------------------------------------------------------------------
def prefix_lists(p, a, m, M, q=3, want_float=False, want_mask=False, want_int=False):
    """for the prefix positions 0..m-1: dict j -> residues of sum_{i<j} q^(a-1-i) 2^(s_i) mod M
    (global exponents, i.e. already multiplied by q^(a-j)).  Optional values (float64 or exact int64,
    selected by want_float / want_int) and bitmasks."""
    L = {0: np.zeros(1, dtype=np.int64)}
    vdt = np.int64 if want_int else np.float64
    F = {0: np.zeros(1, dtype=vdt)} if (want_float or want_int) else None
    B = {0: np.zeros(1, dtype=np.int64)} if want_mask else None
    for t in range(m):
        for j in range(min(t, a - 1), -1, -1):
            if j not in L:
                continue
            if j < a - p + t:  # placing a one here still cannot reach a ones: prune (never joined)
                continue
            term = (pow(q, a - 1 - j, M) * pow(2, t, M)) % M
            new = (L[j] + term) % M
            L[j + 1] = np.concatenate([L[j + 1], new]) if (j + 1) in L else new
            if F is not None:
                tv = q ** (a - 1 - j) * (1 << t)
                if want_int:
                    assert tv < (1 << 61)
                nf = F[j] + vdt(tv)
                F[j + 1] = np.concatenate([F[j + 1], nf]) if (j + 1) in F else nf
            if want_mask:
                nb = B[j] + (1 << t)
                B[j + 1] = np.concatenate([B[j + 1], nb]) if (j + 1) in B else nb
        for j in [j for j in L if j < a - p + 1 + t]:  # states that can no longer reach a ones
            del L[j]
            if F is not None:
                del F[j]
            if want_mask:
                del B[j]
    return L, F, B


def suffix_lists(p, a, m, M, q=3, want_float=False, want_mask=False, want_int=False):
    """for the suffix positions m..p-1: dict k -> residues of sum_{i'<k} q^(k-1-i') 2^(t_i') mod M
    (exponents counted from the end, so independent of the prefix)."""
    L = {0: np.zeros(1, dtype=np.int64)}
    vdt = np.int64 if want_int else np.float64
    F = {0: np.zeros(1, dtype=vdt)} if (want_float or want_int) else None
    B = {0: np.zeros(1, dtype=np.int64)} if want_mask else None
    for t in range(p - 1, m - 1, -1):
        placed = p - 1 - t  # positions already processed after t
        for k in range(min(placed, a - 1), -1, -1):
            if k not in L:
                continue
            term = (pow(q, k, M) * pow(2, t, M)) % M
            new = (L[k] + term) % M
            L[k + 1] = np.concatenate([L[k + 1], new]) if (k + 1) in L else new
            if F is not None:
                tv = q ** k * (1 << t)
                if want_int:
                    assert tv < (1 << 61)
                nf = F[k] + vdt(tv)
                F[k + 1] = np.concatenate([F[k + 1], nf]) if (k + 1) in F else nf
            if want_mask:
                nb = B[k] + (1 << t)
                B[k + 1] = np.concatenate([B[k + 1], nb]) if (k + 1) in B else nb
    return L, F, B


def census_mitm(p, a, q=3, want_words=False):
    """N(p,a) = #{w : |G| divides c_w} by a sorted join of prefix and suffix residues (|G| < 2^62).
    Returns (N, list of integral points x_w = c_w/G) (points only if want_words)."""
    G = gate(p, a, q)
    M = abs(G)
    assert M < (1 << 62)
    if a == 0:
        return 1, [0]
    m = p // 2
    P, _, PB = prefix_lists(p, a, m, M, q, want_mask=want_words)
    Q, _, QB = suffix_lists(p, a, m, M, q, want_mask=want_words)
    N = 0
    pts = []
    for j, pr in P.items():
        k = a - j
        if k not in Q:
            continue
        order = np.argsort(pr, kind="stable")
        ps = pr[order]
        tgt = (M - Q[k]) % M
        lo = np.searchsorted(ps, tgt, "left")
        hi = np.searchsorted(ps, tgt, "right")
        cnt = hi - lo
        N += int(cnt.sum())
        if want_words and cnt.sum() > 0:
            for yi in np.flatnonzero(cnt):
                for r in range(lo[yi], hi[yi]):
                    mask = int(PB[j][order[r]]) | int(QB[k][yi])
                    pos = [t for t in range(p) if (mask >> t) & 1]
                    c = carry(pos, a, q)
                    assert c % G == 0
                    pts.append(c // G)
    return N, sorted(pts)


def census_range(p, a, q=3, chunk=1 << 20, want_points=True):
    """Independent census: every integral periodic point of the clock lies in
    [ceil(cmin/|G|), floor(cmax/|G|)] (in absolute value); iterate the map (qx+1)/2 (G>0) or the
    SHEET-transported map (qy-1)/2 (G<0) for p steps on every candidate and keep the points that
    return with exactly a odd steps.  An orbit leaving [lo, hi] is killed (all orbit points of an
    integral cycle of the clock are periodic points of the same clock)."""
    G = gate(p, a, q)
    if a == 0:
        return 1, [0]
    M = abs(G)
    cmn, cmx = cmin_cmax(p, a, q)
    lo, hi = -(-cmn // M), cmx // M
    sgn = 1 if G > 0 else -1
    if hi < lo:
        return 0, []
    assert hi * q + 1 < (1 << 62)
    N = 0
    pts = []
    x0 = lo
    while x0 <= hi:
        x1 = min(hi, x0 + chunk - 1)
        X = np.arange(x0, x1 + 1, dtype=np.int64)
        Y = X.copy()
        odd = np.zeros_like(X)
        alive = np.ones(X.shape, dtype=bool)
        for _ in range(p):
            b = Y & 1
            odd += b
            Y = np.where(b == 1, (q * Y + sgn) >> 1, Y >> 1)
            alive &= (Y >= lo) & (Y <= hi)
            Y = np.where(alive, Y, lo)
        good = alive & (Y == X) & (odd == a)
        N += int(good.sum())
        if want_points:
            pts.extend(int(v) * sgn for v in X[good])
        x0 = x1 + 1
    return N, sorted(pts)


# ----------------------------------------------------------------------------------------------
# full residue multisets, streamed in chunks
# ----------------------------------------------------------------------------------------------
def residue_chunks(p, a, q=3, max_chunk=1 << 21, want_float=True, want_exact=False):
    """yield (residues mod M as int64, values c_w) over all words of the clock, in chunks.
    values: float64 (want_float), exact int64 (want_exact; requires c_max < 2^62), or None."""
    G = gate(p, a, q)
    M = abs(G)
    assert M < (1 << 62)
    if a == 0:
        yield np.zeros(1, dtype=np.int64), (np.zeros(1, dtype=np.int64) if want_exact else np.zeros(1))
        return
    if want_exact:
        assert cmin_cmax(p, a, q)[1] < (1 << 62)
    m = p // 2
    P, PF, _ = prefix_lists(p, a, m, M, q, want_float=want_float and not want_exact, want_int=want_exact)
    Q, QF, _ = suffix_lists(p, a, m, M, q, want_float=want_float and not want_exact, want_int=want_exact)
    for j in sorted(P):
        k = a - j
        if k not in Q:
            continue
        pr, qr = P[j], Q[k]
        rows = max(1, max_chunk // max(1, len(qr)))
        for r0 in range(0, len(pr), rows):
            blk = (pr[r0:r0 + rows, None] + qr[None, :]) % M
            if PF is not None:
                fb = PF[j][r0:r0 + rows, None] + QF[k][None, :]
                yield blk.ravel(), fb.ravel()
            else:
                yield blk.ravel(), None


def S_abs_all_dp(p, a, q=3, block=4096):
    """|S(h)| for h = 0 .. M//2 by the DP in blocks (memory-lean replacement of an FFT of length M;
    numpy's pocketfft needs ~150 bytes per point for lengths with large prime factors)."""
    M = abs(gate(p, a, q))
    H = M // 2 + 1
    out = np.empty(H, dtype=np.float64)
    for h0 in range(0, H, block):
        h1 = min(H, h0 + block)
        out[h0:h1] = np.abs(S_dp(p, a, range(h0, h1), q))
    return out


def lean_malloc():
    """macOS: re-exec once with MallocLargeCache=0 so that freed large blocks are returned to the
    system (the default large-allocation cache keeps them resident and inflates RSS by up to ~2x)."""
    import os
    if sys.platform == "darwin" and os.environ.get("MallocLargeCache") != "0":
        env = dict(os.environ)
        env["MallocLargeCache"] = "0"
        os.execve(sys.executable, [sys.executable] + sys.argv, env)


def mem_mb():
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / (1024 * 1024) if sys.platform == "darwin" else r / 1024


def report_mem(tag):
    print(f"[mem] {tag}: peak {mem_mb():.0f} MB", file=sys.stderr)
