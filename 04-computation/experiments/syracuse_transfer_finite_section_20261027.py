"""Finite-section (Ulam) transfer operators of the Syracuse map, exact.

Companion to 05-knowledge/results/syracuse_transfer_finite_section_20261027.md (opus-2026-10-27-S65).

Carrier 1 (2-adic).  P_K on the 2^(K-1) odd residues a mod 2^K: P_K(a, b) = Haar measure of
{x = a mod 2^K : S(x) = b mod 2^K} / Haar measure of {x = a mod 2^K}, S(x) = (3x+1)/2^v(x), v = v_2(3x+1).
If v = v_2(3a+1) < K the row is uniform on the 2^v residues b = (3a+1)/2^v (mod 2^(K-v)); for the one residue
a* = -1/3 mod 2^K (v >= K) it is uniform on all odd residues.  The integer counting matrix is C_K = 2^(K-1) P_K:
for every other row C_K(a, b) = #{t mod 2^(K-1) : S(a + 2^K t) = b mod 2^K} (exact, since only t mod 2^v matters).

Carrier 2 (3-adic).  Tao's Syracuse chain on Z/3^J: x -> (3x+1) 2^-a mod 3^J with P(a) = 2^-a, a >= 1.  Its
stationary law is the law of Syrac(Z/3^J) (Tao 2019, Lemma 1.12, Remark 1.13).  P_J(x, y) = 2^(L-a0)/(2^L - 1),
L = 2 3^(J-1), a0 in [1, L] with 2^a0 = (3x+1)/y mod 3^J (y a unit), 0 otherwise.

Sections:
  A  build C_K for K = 8..20: row and column sums (doubly stochastic), nnz = K 2^(K-1)
  B  the class lemma, checked exactly at every level: P_K maps the uniform law on a class mod 2^m to the uniform
     law on a class of strictly smaller precision (or to the stationary law)
  C  exact consequences for K = 8..20: R_n = all classes of precision <= K - n, rank P^n = 2^(K-1-n), the Jordan
     type of P - Pi, the hitting times (h - 1 ~ Binomial(K - 2, 1/2)), worst-case total variation, Dobrushin
     coefficients; hence SLEM = 0 exactly and P^(K-1) = Pi, P^(K-2) != Pi
  D  independent dense check for K = 8..12: C^(K-1) = 2^((K-1)(K-2)) J entrywise, proved modulo primes whose
     product exceeds the entry bound; C^(K-2) differs (certificate prime); K = 8 also in exact Python integers
  E  floating point as a guide only: numpy eigenvalues of P_K report spurious |lambda_2| (Jordan blocks)
  F  carrier 2 for J = 1..6 exact (J = 7 modulo primes): P^J = Pi, P^(J-1) != Pi, rank P^n = 3^(J-n), the
     stationary law (Tao's mod-9 law at J = 2), max primitive |Fourier coefficient| of Syrac(Z/3^J) against THM-4520
"""
import itertools
import math
import random
import sys
import time
from fractions import Fraction

import numpy as np

KS = range(8, 21)


# ============================================================================ carrier 1: the 2-adic Ulam matrix
def v2(t):
    """2-adic valuation of a positive int64 array (exact: the lowest set bit is a power of two)."""
    low = t & -t
    return (np.frexp(low.astype(np.float64))[1] - 1).astype(np.int64)


def phi(c, m):
    """The class map: the uniform law on {x = c mod 2^m} (odd c) goes to the uniform law on {x = c' mod 2^m'}.
    If v = v_2(3c+1) < m: c' = (3c+1)/2^v mod 2^(m-v), m' = m - v.  Otherwise the image is the stationary law,
    encoded as the root class (1, 1)."""
    t = 3 * c + 1
    v = v2(t)
    det = v < m
    m2 = np.where(det, m - v, 1)
    c2 = np.where(det, (t >> np.minimum(v, 62)) & ((np.int64(1) << m2) - 1), 1)
    return c2.astype(np.int64), m2.astype(np.int64)


def node(c, m):
    """Heap-style id of the class (c mod 2^m): 2^(m-1) + (c-1)/2, in [1, 2^K)."""
    return (np.int64(1) << (m - 1)) + ((c - 1) >> 1)


def build_counting_matrix(K):
    """C_K = 2^(K-1) P_K in CSR form (indptr, indices, data); rows/columns indexed by i = (a-1)/2."""
    n = 1 << (K - 1)
    a = 2 * np.arange(n, dtype=np.int64) + 1
    t = 3 * a + 1
    v = v2(t)
    special = v >= K
    vv = np.minimum(v, K - 1)
    counts = np.where(special, n, np.int64(1) << vv)
    indptr = np.zeros(n + 1, dtype=np.int64)
    np.cumsum(counts, out=indptr[1:])
    nnz = int(indptr[-1])
    rows = np.repeat(np.arange(n, dtype=np.int64), counts)
    j = np.arange(nnz, dtype=np.int64) - indptr[rows]
    r = np.where(special, 1, (t >> np.minimum(v, 62)) & ((np.int64(1) << (K - vv)) - 1))   # odd, mod 2^(K-v)
    start = np.where(special, 0, (r - 1) >> 1)
    stride = np.where(special, 1, np.int64(1) << np.maximum(K - vv - 1, 0))
    indices = start[rows] + j * stride[rows]
    data = np.where(special, 1, np.int64(1) << (K - 1 - vv))[rows]
    return indptr, indices, data, a, special


def section_A_B(K):
    n = 1 << (K - 1)
    t0 = time.time()
    indptr, indices, data, a, special = build_counting_matrix(K)
    nnz = len(indices)
    rowsum = np.add.reduceat(data, indptr[:-1])
    colsum = np.zeros(n, dtype=np.int64)
    np.add.at(colsum, indices, data)
    assert nnz == K * n and int(special.sum()) == 1
    assert np.all(rowsum == n) and np.all(colsum == n)
    # (B1) every row is the uniform law on the class phi(a, K): support = the class, values constant
    c1, m1 = phi(a, np.full(n, K, dtype=np.int64))
    size = np.int64(1) << (K - m1)                      # residues mod 2^K in the class
    assert np.all(np.diff(indptr) == size)
    rows = np.repeat(np.arange(n, dtype=np.int64), np.diff(indptr))
    jj = np.arange(nnz, dtype=np.int64) - indptr[rows]
    expect = ((c1 - 1) >> 1)[rows] + jj * (np.int64(1) << (m1 - 1))[rows]
    assert np.array_equal(indices, expect) and np.all(data == (n // size)[rows])
    # (B2) the merge identity at every coarser level: mu(c,m) = (mu(c,m+1) + mu(c+2^m,m+1))/2, so the lemma
    # P mu(c,m) = mu(phi(c,m)) holds at level m iff the two children's images are equal to phi(c,m) or are its
    # two children
    for m in range(K - 1, 0, -1):
        c = 2 * np.arange(1 << (m - 1), dtype=np.int64) + 1
        mm = np.full(len(c), m, dtype=np.int64)
        p0c, p0m = phi(c, mm)
        pac, pam = phi(c, mm + 1)
        pbc, pbm = phi(c + (np.int64(1) << m), mm + 1)
        same = (pac == p0c) & (pam == p0m) & (pbc == p0c) & (pbm == p0m)
        kids = (pam == p0m + 1) & (pbm == p0m + 1) & (((pac - p0c) & ((np.int64(1) << p0m) - 1)) == 0) & \
               (((pbc - p0c) & ((np.int64(1) << p0m) - 1)) == 0) & (pac != pbc)
        assert np.all(same | kids), (K, m)
    return dict(n=n, nnz=nnz, secs=time.time() - t0, special=int(a[special][0]))


def tree_rank(K, members):
    """dim span{1_C : C in members} for classes C given by node ids; bottom-up over the binary tree of classes."""
    inS = np.zeros(1 << K, dtype=bool)
    inS[members] = True
    full = np.zeros(1 << K, dtype=bool)
    dim = np.zeros(1 << K, dtype=np.int64)
    lv = np.arange(1 << (K - 1), dtype=np.int64) + (1 << (K - 1))
    full[lv] = inS[lv]
    dim[lv] = inS[lv]
    for m in range(K - 1, 0, -1):
        t = np.arange(1 << (m - 1), dtype=np.int64)
        nd = (1 << (m - 1)) + t
        L = (1 << m) + t
        R = L + (1 << (m - 1))
        both = full[L] & full[R]
        dim[nd] = dim[L] + dim[R] + (inS[nd] & ~both)
        full[nd] = inS[nd] | both
    return int(dim[1])


def section_C(K):
    n = 1 << (K - 1)
    c = 2 * np.arange(n, dtype=np.int64) + 1
    m = np.full(n, K, dtype=np.int64)
    hit = np.full(n, -1, dtype=np.int64)
    ranks, worst_tv, avg_tv, dobrushin = [], [], [], []
    for step in range(0, K):
        ids = np.unique(node(c, m))
        levels = np.unique(m)
        if step == 0:   # R_0 = the 2^(K-1) fine classes
            assert len(ids) == n and np.all(levels == K)
        else:           # R_n = every class of precision 1..K-n (all 2^(K-n) - 1 of them)
            assert len(ids) == (1 << (K - step)) - 1 and np.array_equal(levels, np.arange(1, K - step + 1)), (K, step)
        r = tree_rank(K, ids)
        assert r == 1 << (K - 1 - step)
        ranks.append(r)
        mx = int(m.max())
        worst_tv.append(Fraction(1) - Fraction(2) ** (1 - mx))
        avg_tv.append(sum(Fraction(int(cnt), n) * (1 - Fraction(2) ** (1 - int(mv)))
                          for mv, cnt in zip(*np.unique(m, return_counts=True))))
        # Dobrushin: two incomparable classes exist iff some precision >= 2 (two distinct level-2 classes)
        dobrushin.append(1 if mx >= 2 else 0)
        newly = (m == 1) & (hit < 0)
        hit[newly] = step
        c, m = phi(c, m)
    assert np.all(hit >= 0) and hit.max() == K - 1
    dist = np.bincount(hit)
    assert all(int(dist[h]) == 2 * math.comb(K - 2, h - 1) for h in range(1, K)) and dist[0] == 0
    assert all(worst_tv[s] == 1 - Fraction(1, 2 ** (K - 1 - s)) for s in range(K))
    assert dobrushin == [1] * (K - 1) + [0]
    # Jordan type of N = P - Pi on W = (mean zero): blocks of size >= s = r_(s-1) - r_s, r_s = rank(P^s) - 1
    rr = [x - 1 for x in ranks] + [0]
    ge = [rr[s - 1] - rr[s] for s in range(1, K)]
    exact = {s: ge[s - 1] - (ge[s] if s < K - 1 else 0) for s in range(1, K)}
    assert exact == {**{s: 2 ** (K - 2 - s) for s in range(1, K - 1)}, K - 1: 1}
    mean_hit = Fraction(int((hit * 1).sum()), n)
    assert mean_hit == Fraction(K, 2)
    return dict(ranks=ranks, worst_tv=worst_tv, avg_tv=avg_tv, jordan=exact, mean_hit=mean_hit,
                dist=[int(x) for x in dist[1:]])


# ============================================================================ D, E: dense checks and float guide
SMALL_PRIMES = [p for p in range(2 ** 20 - 1, 2 ** 19, -2) if all(p % q for q in range(3, int(p ** 0.5) + 1, 2))][:12]


def dense_C(K):
    indptr, indices, data, a, special = build_counting_matrix(K)
    n = 1 << (K - 1)
    C = np.zeros((n, n), dtype=np.int64)
    rows = np.repeat(np.arange(n), np.diff(indptr))
    C[rows, indices] = data
    return C


def matpow_mod(C, e, p):
    A = (C % p).astype(np.float64)
    R = np.eye(len(C), dtype=np.float64)
    for _ in range(e):
        R = np.fmod(R @ A, p)          # every partial sum is an integer below 2^53: exact
    return R


def section_D(K):
    """C^(K-1) = 2^((K-1)(K-2)) * J entrywise; entries lie in [0, 2^((K-1)^2)], so agreement modulo primes with
    product above 2^((K-1)^2 + 1) proves equality.  C^(K-2) != 2^((K-1)(K-3)) J: one prime suffices."""
    C = dense_C(K)
    n = len(C)
    assert n * (SMALL_PRIMES[0] ** 2) < 2 ** 53
    need = (K - 1) ** 2 + 1
    used, prod = [], 1
    for p in SMALL_PRIMES:
        if prod > 2 ** need:
            break
        R = matpow_mod(C, K - 1, p)
        target = pow(2, (K - 1) * (K - 2), p)
        assert np.all(R == target), (K, p)
        used.append(p); prod *= p
    assert prod > 2 ** need
    p = SMALL_PRIMES[0]
    R = matpow_mod(C, K - 2, p)
    differs = bool(np.any(R != pow(2, (K - 1) * (K - 3), p)))
    assert differs
    return len(used)


def section_D_exact_K8():
    C = dense_C(8).tolist()
    n = len(C)
    M = [row[:] for row in C]
    for _ in range(6):                               # C^7
        M = [[sum(M[i][k] * C[k][j] for k in range(n) if C[k][j]) for j in range(n)] for i in range(n)]
    assert all(x == 2 ** 42 for row in M for x in row)         # (K-1)(K-2) = 42
    return True


def section_E(K):
    C = dense_C(K).astype(np.float64)
    P = C / C.sum(axis=1)[:, None]
    ev = np.linalg.eigvals(P)
    ev = sorted(ev, key=lambda z: -abs(z))
    assert abs(ev[0] - 1) < 1e-9
    return abs(ev[1])


# ============================================================================ F: carrier 2, Tao's chain mod 3^J
def dlog_table(J):
    M3 = 3 ** J
    L = 2 * 3 ** (J - 1)
    dl, x = {}, 1
    for k in range(L):
        dl[x] = k
        x = 2 * x % M3
    assert len(dl) == L                               # 2 is a primitive root mod 3^J
    return dl, L


def tao_matrix_int(J):
    """M = (2^L - 1) P as a dict-of-rows of exact integers (rows: all residues mod 3^J; columns: units)."""
    M3 = 3 ** J
    dl, L = dlog_table(J)
    units = sorted(dl)
    rows = []
    for x in range(M3):
        u = (3 * x + 1) % M3
        du = dl[u]
        row = {}
        for y in units:
            a0 = (du - dl[y]) % L
            if a0 == 0:
                a0 = L
            row[y] = 2 ** (L - a0)
        rows.append(row)
    return rows, L


def law_after(rows, start, steps, M3):
    """Exact law (integer numerators over (2^L - 1)^steps) of x_steps from x_0 = start."""
    d = {start: 1}
    for _ in range(steps):
        nd = {}
        for x, w in d.items():
            for y, m in rows[x].items():
                nd[y] = nd.get(y, 0) + w * m
        d = nd
    return d


def section_F_exact(J):
    M3 = 3 ** J
    rows, L = tao_matrix_int(J)
    den = 2 ** L - 1
    assert all(sum(r.values()) == den for r in rows)  # stochastic: sum_a 2^(L-a) over a = 1..L
    pi = law_after(rows, 0, J, M3)
    starts = list(range(M3)) if J <= 4 else [0, 1, 2, 5, M3 - 1]
    same_J = all(law_after(rows, s, J, M3) == pi for s in starts)
    l0 = law_after(rows, 0, J - 1, M3)
    differs = any(law_after(rows, s, J - 1, M3) != l0 for s in starts[1:])
    assert same_J and differs
    total = den ** J
    law = {y: Fraction(w, total) for y, w in pi.items()}
    assert sum(law.values()) == 1
    return law, L, len(starts)


def modular_matrix(J, q):
    M3 = 3 ** J
    dl, L = dlog_table(J)
    inv = pow((pow(2, L, q) - 1) % q, -1, q)
    pw = [pow(2, L - a, q) * inv % q for a in range(L + 1)]
    units = np.array(sorted(dl), dtype=np.int64)
    dlu = np.array([dl[int(y)] for y in units], dtype=np.int64)
    P = np.zeros((M3, M3), dtype=np.float64)
    for x in range(M3):
        du = dl[(3 * x + 1) % M3]
        a0 = (du - dlu) % L
        a0[a0 == 0] = L
        P[x, units] = np.array([pw[int(a)] for a in a0], dtype=np.float64)
    return P


def rank_mod(A, q):
    A = (A % q).astype(np.int64).copy()
    nr, nc = A.shape
    r = 0
    for col in range(nc):
        piv = next((i for i in range(r, nr) if A[i, col] % q), None)
        if piv is None:
            continue
        A[[r, piv]] = A[[piv, r]]
        inv = pow(int(A[r, col]), -1, q)
        A[r] = A[r] * inv % q
        nz = np.nonzero(A[:, col])[0]
        for i in nz:
            if i != r:
                A[i] = (A[i] - A[i, col] * A[r]) % q
        r += 1
        if r == nr:
            break
    return r


def section_F_modular(J, q, ranks=True):
    """Rows of P^J equal mod q (consistent with the proof), rows of P^(J-1) differ mod q (rigorous certificate
    that P^(J-1) != Pi over Q), and rank_q(P^n) = 3^(J-n): a lower bound on the rational rank that meets the
    structural upper bound (rows of P^n depend on x_0 only mod 3^(J-n)), hence the exact rank."""
    P = modular_matrix(J, q)
    n = len(P)
    assert n * q * q < 2 ** 53
    R = np.eye(n)
    powers = [R]
    for s in range(J):
        R = np.fmod(R @ P, q)
        powers.append(R)
    assert np.all(powers[J] == powers[J][0])
    assert np.any(powers[J - 1] != powers[J - 1][0])
    rk = None
    if ranks:
        rk = [rank_mod(powers[s].astype(np.int64), q) for s in range(J + 1)]
        assert rk == [3 ** (J - s) for s in range(J + 1)], rk
    return rk


def fourier_max(law, J):
    M3 = 3 ** J
    ys = np.array(list(law.keys()), dtype=np.float64)
    ws = np.array([float(w) for w in law.values()])
    best = 0.0
    for xi in range(1, M3):
        if xi % 3 == 0:
            continue
        best = max(best, abs(np.sum(ws * np.exp(2j * np.pi * xi * ys / M3))))
    return best


def section_G_tmap_3adic(J):
    """Contrast: the T-map random chain on Z/3^J, x -> x/2 or (3x+1)/2 with probability 1/2 each (2^-1 mod 3^J).
    Exact characteristic polynomial: det(lambda I - P) = (lambda - 1) prod_(d=1..J) (lambda^(L_d) - 2^(-L_d)),
    L_d = 2 3^(d-1): the depth-preserving halving branch is 1/2 times the doubling cycle on the units mod 3^d
    (THM-4631 at u = v = 1), so the second-largest eigenvalue modulus is exactly 1/2."""
    import sympy as sp
    M3 = 3 ** J
    h = pow(2, -1, M3)
    P = sp.zeros(M3, M3)
    for x in range(M3):
        P[x, x * h % M3] += sp.Rational(1, 2)
        P[x, (3 * x + 1) * h % M3] += sp.Rational(1, 2)
    lam = sp.symbols('lam')
    cp = sp.expand(P.charpoly(lam).as_expr())
    target = sp.expand((lam - 1) * sp.prod([lam ** (2 * 3 ** (d - 1)) - sp.Rational(1, 2 ** (2 * 3 ** (d - 1)))
                                            for d in range(1, J + 1)]))
    return sp.simplify(cp - target) == 0


def section_G_tmap_2adic(K):
    """Contrast: the T-map Ulam matrix on all residues mod 2^K (de Bruijn): row a is uniform on the two lifts of
    T(a) mod 2^(K-1).  Exact integers: (2P)^K = 2^(K-1) J? -> P^K = Pi, P^(K-1) != Pi (THM-4629 2(c))."""
    n = 1 << K
    C = [[0] * n for _ in range(n)]
    for a in range(n):
        t = a // 2 if a % 2 == 0 else (3 * a + 1) // 2
        base = t % (n // 2)
        C[a][base] += 1
        C[a][base + n // 2] += 1
    cols = [[k for k in range(n) if C[k][j]] for j in range(n)]
    M = [row[:] for row in C]
    powers = [None, [row[:] for row in C]]
    for _ in range(K - 1):
        M = [[sum(M[i][k] * C[k][j] for k in cols[j]) for j in range(n)] for i in range(n)]
        powers.append(M)
    allsame_K = all(row == powers[K][0] for row in powers[K])
    allsame_K1 = all(row == powers[K - 1][0] for row in powers[K - 1])
    return allsame_K and not allsame_K1


def section_C2(K):
    """l^2(pi) norms: ||(P - Pi)^n|| = sqrt(P(Bin(K-2, 1/2) >= n)) (proved in the note: P^n f = sum over levels of
    conditional expectations with level weights; the maximiser is f = +-1 on 1 mod 4 / 3 mod 4).  Exact value of
    ||P^n f0||^2 for that f0, and a float SVD as a guide that nothing larger occurs."""
    n = 1 << (K - 1)
    exact = [sum(Fraction(math.comb(K - 2, j), 2 ** (K - 2)) for j in range(s, K - 1)) for s in range(0, K)]
    # ||P^s f0||^2 exactly through the classes: f0 averages to +-1 on every class of precision >= 2, to 0 on the root
    c = 2 * np.arange(n, dtype=np.int64) + 1
    m = np.full(n, K, dtype=np.int64)
    f0vals = []
    for s in range(K):
        f0vals.append(Fraction(int(np.sum(m >= 2)), n))
        c, m = phi(c, m)
    assert f0vals == exact
    svd = None
    if K <= 12:
        C = dense_C(K).astype(np.float64)
        N = C / n - 1.0 / n
        Nk = np.eye(n)
        svd = []
        for s in range(1, K):
            Nk = Nk @ N
            svd.append(np.linalg.norm(Nk, 2) ** 2)
        assert all(abs(svd[s - 1] - float(exact[s])) < 1e-9 for s in range(1, K))
    return exact, svd


def section_E3(J, law):
    """Float guide for the 3-adic chain: numpy's |lambda_2| of P_J, and the norm of P_J - Pi on l^2(pi) restricted to
    the units (the support of pi; the non-units are transient)."""
    M3 = 3 ** J
    rows, L = tao_matrix_int(J)
    units = sorted(rows[0])
    P = np.zeros((M3, M3))
    for x, r in enumerate(rows):
        for y, w in r.items():
            P[x, y] = float(Fraction(w, 2 ** L - 1))
    ev = sorted(np.linalg.eigvals(P), key=lambda z: -abs(z))
    pu = np.array([float(law[y]) for y in units])
    Q = P[np.ix_(units, units)]
    Nq = Q - pu[None, :]
    S = np.sqrt(pu)[:, None] * Nq / np.sqrt(pu)[None, :]
    return abs(ev[1]), np.linalg.norm(S, 2)


def phi_p(c, m, p):
    """The class map for x -> (p x + 1)/2^v (p odd): p x + 1 is a bijection from odd to even residues mod 2^m."""
    t = p * c + 1
    t = np.where(t > 0, t, t + (np.int64(1) << 40))       # p = -1: keep a positive representative (2-adic value)
    v = v2(t)
    det = v < m
    m2 = np.where(det, m - v, 1)
    c2 = np.where(det, (t >> np.minimum(v, 62)) & ((np.int64(1) << m2) - 1), 1)
    return c2.astype(np.int64), m2.astype(np.int64)


def section_H_2adic(K, p):
    """Universality: for x -> (p x + 1)/2^v the reachable classes, ranks and hitting times are exactly those of p = 3."""
    n = 1 << (K - 1)
    c = 2 * np.arange(n, dtype=np.int64) + 1
    m = np.full(n, K, dtype=np.int64)
    hit = np.full(n, -1, dtype=np.int64)
    for step in range(K):
        ids = np.unique(node(c, m))
        if step >= 1:
            assert len(ids) == (1 << (K - step)) - 1
        assert tree_rank(K, ids) == 1 << (K - 1 - step)
        hit[(m == 1) & (hit < 0)] = step
        c, m = phi_p(c, m, p)
    dist = np.bincount(hit)
    return all(int(dist[h]) == 2 * math.comb(K - 2, h - 1) for h in range(1, K))


def section_H_padic(p, J, q=1048573):
    """Universality: the chain x -> (p x + 1) 2^-a on Z/p^J (a ~ Geom(1/2)): P^J = Pi, P^(J-1) != Pi and
    rank P^n = p^(J-n) (modular lower bound = structural upper bound), whether or not 2 is a primitive root."""
    M = p ** J
    L = 1
    while pow(2, L, M) != 1:
        L += 1
    inv = pow((pow(2, L, q) - 1) % q, -1, q)
    h = pow(2, -1, M)
    P = np.zeros((M, M), dtype=np.float64)
    for x in range(M):
        y = (p * x + 1) % M
        for a in range(1, L + 1):
            y = y * h % M                                 # (p x + 1) 2^-a
            P[x, y] = (P[x, y] + pow(2, L - a, q) * inv) % q
    assert M * q * q < 2 ** 53
    R = np.eye(M)
    ranks = [M]
    for s in range(1, J + 1):
        R = np.fmod(R @ P, q)
        ranks.append(rank_mod(R.astype(np.int64), q))
        if s == J - 1:
            assert np.any(R != R[0])
    assert np.all(R == R[0])
    assert ranks == [p ** (J - s) for s in range(J + 1)], ranks
    return L, ranks


def main():
    print('== A, B. the 2-adic Ulam matrices P_K: exact construction and the class lemma')
    for K in KS:
        info = section_A_B(K)
        print(f'  K = {K:2d}: {info["n"]:7d} states, nnz = K 2^(K-1) = {info["nnz"]:9d}; row sums = column sums = '
              f'2^(K-1) (doubly stochastic); special residue a* = {info["special"]} (= -1/3 mod 2^K); class lemma '
              f'checked at all {K} levels  [{info["secs"]:.1f}s]')
    print('== C. exact consequences: R_n, ranks, Jordan type, hitting times, total variation, Dobrushin')
    for K in KS:
        info = section_C(K)
        jt = info['jordan']
        print(f'  K = {K:2d}: rank P^n = 2^(K-1-n) for n = 0..K-1 (checked); Jordan type of P - Pi on the mean-zero space: '
              f'one block of size {K - 1} and 2^(K-2-s) blocks of size s for s = 1..{K - 2} '
              f'(sizes 1..4: {[jt[s] for s in range(1, min(5, K - 1))]}); '
              f'mean exact mixing time = K/2 = {info["mean_hit"]}; '
              f'worst TV at n = K-2: {info["worst_tv"][K - 2]}; average TV at n = K/2: {float(info["avg_tv"][K // 2]):.4f}')
    K = 12
    info = section_C(K)
    print(f'  K = 12 hitting-time counts (h = 1..11): {info["dist"]} = 2 C(10, h-1)')
    print(f'  K = 12 average TV profile n = 0..11: {[round(float(x), 4) for x in info["avg_tv"]]}')
    print('  -> for every K: P_K = Pi + N_K with N_K nilpotent of index exactly K - 1; spectrum {1, 0, ..., 0};')
    print('     second-largest eigenvalue modulus 0; Dobrushin coefficient of P^n is 1 for n <= K - 2 and 0 at n = K - 1')
    print('== C2. the l^2(pi) operator norms: ||(P_K - Pi)^n||^2 = P(Bin(K-2, 1/2) >= n)')
    for K in KS:
        exact, svd = section_C2(K)
        s1 = math.sqrt(float(exact[1]))
        print(f'  K = {K:2d}: ||P_K - Pi||^2 = 1 - 2^-(K-2) = {exact[1]} (norm {s1:.7f}; singular-value gap '
              f'1 - norm = {1 - s1:.3e} ~ 2^-(K-1) = {2.0 ** -(K - 1):.3e})' +
              ('; float SVD of all powers agrees' if svd else ''))
    print('  -> the eigenvalue gap is 1 for every K, but P_K is far from normal: its second singular value tends to 1')
    print('     like 1 - 2^-(K-1), and the l^2 profile of the worst function is the binomial tail (cutoff at n = K/2)')
    print('== D. independent dense checks (no class machinery)')
    for K in range(8, 13):
        k = section_D(K)
        print(f'  K = {K}: C^(K-1) = 2^((K-1)(K-2)) J proved by congruence modulo {k} primes (product > 2^((K-1)^2+1)); '
              f'C^(K-2) != 2^((K-1)(K-3)) J (certificate prime {SMALL_PRIMES[0]})')
    print('  K = 8 in exact Python integers: C^7 = 2^42 J:', section_D_exact_K8())
    print('== E. floating point (a guide only): numpy eigenvalues of P_K')
    for K in range(8, 13):
        s = section_E(K)
        print(f'  K = {K}: numpy reports |lambda_2| = {s:.4f} (exact value 0; a Jordan block of size {K - 1} perturbed '
              f'by 1e-16 has eigenvalues of modulus about 1e-16^(1/{K - 1}) = {1e-16 ** (1 / (K - 1)):.4f})')
    print('== F. carrier 2: Tao\'s Syracuse chain on Z/3^J')
    known = {2: [0.37792], 3: [0.25224], 4: [0.17700], 5: [0.12927], 6: [0.09611]}
    laws = {}
    for J in range(1, 7):
        law, L, nstarts = section_F_exact(J)
        laws[J] = law
        tv = Fraction(1, 2) * sum(abs(law.get(y, 0) - (Fraction(1, 2 * 3 ** (J - 1)) if y % 3 else 0))
                                  for y in range(3 ** J))
        line = (f'  J = {J}: {3 ** J} states, L = {L}; rows of P^J equal (exact, {nstarts} starts), rows of P^(J-1) '
                f'differ; TV(Syrac(Z/3^J), uniform on units) = {float(tv):.5f}')
        if J >= 2:
            fm = fourier_max(law, J)
            line += f'; max primitive |hat mu_J| = {fm:.5f}'
            if J in known:
                line += f' (THM-4520 (5) table: {known[J][0]})'
        print(line)
        if J == 2:
            tao = [Fraction(x, 63) for x in (0, 8, 16, 0, 11, 4, 0, 2, 22)]
            got = [law.get(y, Fraction(0)) for y in range(9)]
            print(f'    stationary law mod 9 = {[str(g) for g in got]}; Tao\'s mod-9 law: {got == tao}')
            assert got == tao
    q = 1048573
    for J in range(2, 9):
        rk = section_F_modular(J, q, ranks=(J <= 6))
        print(f'  J = {J} (mod q = {q}): rows of P^J equal mod q, rows of P^(J-1) differ mod q (certificate)' +
              (f'; rank P^n = {rk} = 3^(J-n): the modular lower bound meets the structural upper bound' if rk else ''))
    print('  Jordan type of P_J - Pi (from rank P^n = 3^(J-n), proved for all J in the note): two blocks of size J and')
    print('  4 3^(J-1-s) blocks of size s for s = 1..J-1, e.g. J = 6:', {s: 4 * 3 ** (6 - 1 - s) for s in range(1, 6)} | {6: 2})
    print('  float guide (numpy): |lambda_2| of P_J and the l^2(pi) norm of P_J - Pi on the units:')
    for J in range(2, 7):
        lam2, nrm = section_E3(J, laws[J])
        print(f'    J = {J}: numpy |lambda_2| = {lam2:.4f} (exact 0); ||P_J - Pi||_(l^2(pi)) = {nrm:.6f}')
    print('== G. contrast: the T-map (x/2 or (3x+1)/2 with fair coins)')
    for J in (1, 2, 3):
        ok = section_G_tmap_3adic(J)
        print(f'  Z/3^{J}: det(lambda I - P) = (lambda - 1) prod_(d=1..{J}) (lambda^(L_d) - 2^-(L_d)) exactly: {ok} '
              f'-> second-largest eigenvalue modulus exactly 1/2 (the doubling cycles of THM-4631 at u = v = 1)')
        assert ok
    for K in (6, 7, 8):
        ok = section_G_tmap_2adic(K)
        print(f'  Z/2^{K} (de Bruijn Ulam matrix, exact integers): P^K = Pi and P^(K-1) != Pi: {ok} (THM-4629 2(c))')
        assert ok
    print('== H. universality in the multiplier p (the finite sections cannot see p = 3)')
    for p in (1, 5, 7, -1):
        oks = [section_H_2adic(K, p) for K in range(8, 15)]
        print(f'  x -> ({p} x + 1)/2^v on odd residues mod 2^K, K = 8..14: same reachable classes, ranks 2^(K-1-n), '
              f'Jordan type and hitting times as p = 3: {all(oks)}')
        assert all(oks)
    for p, Js in ((5, (1, 2, 3)), (7, (1, 2, 3))):
        for J in Js:
            L, rk = section_H_padic(p, J)
            print(f'  x -> ({p} x + 1) 2^-a on Z/{p}^{J} (ord of 2 = {L}{", 2 not primitive" if L < (p - 1) * p ** (J - 1) else ""}): '
                  f'P^J = Pi, P^(J-1) != Pi (mod q), rank P^n = {rk} = {p}^(J-n)')
    print('ALL CHECKS PASSED')


if __name__ == '__main__':
    main()
