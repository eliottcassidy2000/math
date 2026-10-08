#!/usr/bin/env python3
"""audit H: independent vectorized hidden-state DP and exact balance test (no imports from the author's modules).

Coordinates ("position basis"): the debt lattice is spanned by the log-vectors w_j = log(m_j/m_0); take as basis the first
rho of w_1, w_2, ... (position order) that are linearly independent, and express every w_j in it (rational solve, scaled to
integers).  For the maps audited here this is the convention under which the author's integer forms are stated; for maps
whose non-unit multipliers are distinct primes in increasing position order it is the prime-exponent basis.
"""
import itertools, math
from fractions import Fraction as Fr
import numpy as np
from hcore import factor, primes_of, expvec, lattice_rank, ratio_group, std_r


def position_coords(m):
    primes = primes_of(m)
    E = [expvec(x, primes) for x in m]
    w = [tuple(a - b for a, b in zip(E[j], E[0])) for j in range(len(m))]
    basis = []
    for j in range(1, len(m)):
        if lattice_rank(basis + [w[j]]) == len(basis) + 1:
            basis.append(w[j])
    rho = len(basis)
    # solve w_j = sum_k c_k basis_k by least squares over Q (exact): normal equations
    B = [[Fr(x) for x in b] for b in basis]
    G = [[sum(B[a][t] * B[c][t] for t in range(len(primes))) for c in range(rho)] for a in range(rho)]
    Ginv = inv_frac(G)
    coords = []
    for j in range(len(m)):
        rhs = [sum(B[a][t] * w[j][t] for t in range(len(primes))) for a in range(rho)]
        c = [sum(Ginv[a][b] * rhs[b] for b in range(rho)) for a in range(rho)]
        # exactness check: reconstruct
        rec = [sum(c[a] * B[a][t] for a in range(rho)) for t in range(len(primes))]
        assert rec == [Fr(x) for x in w[j]], "not in the span"
        coords.append(c)
    den = 1
    for c in coords:
        for x in c:
            den = den * x.denominator // math.gcd(den, x.denominator)
    return rho, [tuple(int(x * den) for x in c) for c in coords], den


def inv_frac(G):
    n = len(G)
    A = [list(row) + [Fr(int(i == j)) for j in range(n)] for i, row in enumerate(G)]
    for c in range(n):
        p = next(i for i in range(c, n) if A[i][c] != 0)
        A[c], A[p] = A[p], A[c]
        pv = A[c][c]
        A[c] = [x / pv for x in A[c]]
        for i in range(n):
            if i != c and A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[c])]
    return [row[n:] for row in A]


class Tables:
    def __init__(self, d, m, r=None):
        self.d, self.m = d, list(m)
        self.r = list(r) if r is not None else std_r(m)
        for i in range(d):
            assert math.gcd(m[i], d) == 1 and (m[i] * i + self.r[i]) % d == 0
        self.rho, self.v, self.den = position_coords(m)
        n = self.rho
        self.Dab = np.zeros((d, d, n, n), dtype=np.int64)
        for a in range(d):
            if math.gcd(a, d) != 1:
                continue
            for b in range(d):
                Dm = np.zeros((n, n), dtype=np.int64)
                for j in range(d):
                    xi = np.array(self.v[(a * j + b) % d], dtype=np.int64) - np.array(self.v[j], dtype=np.int64)
                    Dm += np.outer(xi, xi)
                self.Dab[a, b] = Dm
        self.H = {}
        self.T = {}

    def Hs(self, s):
        if s not in self.H:
            self.H[s] = np.array(sorted(ratio_group(self.m, self.d ** s)), dtype=np.int64)
        return self.H[s]

    def blocks_at(self, s, Mv, Ev):
        """A_s at the nodes (Mv[i], Ev[i]) (given mod d^s), from the full level-(s-1) table."""
        d = self.d
        a = Mv % d
        b = Ev % d
        out = self.Dab[a, b] * (d ** (s - 1))
        if s == 1:
            return out
        mod, pmod = d ** s, d ** (s - 1)
        Tp = self.table(s - 1)
        Hp = self.Hs(s - 1)
        lookup = -np.ones(pmod, dtype=np.int64)
        lookup[Hp] = np.arange(len(Hp))
        m = np.array(self.m, dtype=np.int64)
        r = np.array(self.r, dtype=np.int64)
        minv = np.array([pow(int(x), -1, mod) for x in self.m], dtype=np.int64)
        for j in range(d):
            i = (a * j + b) % d
            ratio = (m[i] * minv[j]) % mod                  # m_i/m_j mod d^s
            N = (m[i] * Ev + r[i] - (((ratio * Mv) % mod) * r[j]) % mod) % mod
            assert np.all(N % d == 0)
            e2 = (N // d) % pmod
            M2 = (ratio * Mv) % pmod
            ix = lookup[M2]
            assert np.all(ix >= 0)
            out = out + Tp[ix, e2]
        return out

    def table(self, s):
        if s not in self.T:
            H = self.Hs(s)
            Mv = np.repeat(H, self.d ** s)
            Ev = np.tile(np.arange(self.d ** s, dtype=np.int64), len(H))
            self.T[s] = self.blocks_at(s, Mv, Ev).reshape(len(H), self.d ** s, self.rho, self.rho)
        return self.T[s]


def det3_obj(S):
    """exact 3x3 determinants of an (N,3,3) int64 array via Python ints (object arrays)."""
    o = S.astype(object)
    a, b, c = o[:, 0, 0], o[:, 0, 1], o[:, 0, 2]
    dd, e, f = o[:, 1, 0], o[:, 1, 1], o[:, 1, 2]
    g, h, i = o[:, 2, 0], o[:, 2, 1], o[:, 2, 2]
    return a * (e * i - f * h) - b * (dd * i - f * g) + c * (dd * h - e * g)


def pd_charpoly_vec(S):
    """Exact PD test of symmetric integer matrices S (N,n,n), n = 3 or 4, by the elementary symmetric functions
    e_1..e_n of the eigenvalues (sums of principal minors): PD iff all > 0.  Python-int arithmetic."""
    N, n, _ = S.shape
    o = S.astype(object)
    ok = np.ones(N, dtype=bool)
    for k in range(1, n + 1):
        e = np.zeros(N, dtype=object)
        for idx in itertools.combinations(range(n), k):
            sub = o[:, list(idx)][:, :, list(idx)]
            e = e + det_obj(sub)
        ok &= np.array([x > 0 for x in e], dtype=bool)
    return ok


def det_obj(sub):
    k = sub.shape[1]
    if k == 1:
        return sub[:, 0, 0]
    if k == 2:
        return sub[:, 0, 0] * sub[:, 1, 1] - sub[:, 0, 1] * sub[:, 1, 0]
    tot = 0
    for c in range(k):
        minor = np.delete(np.delete(sub, 0, axis=1), c, axis=2)
        tot = tot + ((-1) ** c) * sub[:, 0, c] * det_obj(minor)
    return tot


def adj_det(Q):
    from hcore import adj_int, det_int
    return adj_int(Q), det_int(Q)


def balance_exact(Q, A):
    """A: (N,n,n) int64 symmetric PSD blocks. Returns boolean array: tr(A Q^-1) > 2 lambda_max(A Q^-1)  exactly,
    via S' = tr(A adj Q) Q - 2 det(Q) A  positive definite (charpoly-sign test)."""
    adjQ, detQ = adj_det(Q)
    assert detQ > 0
    adjQ = np.array(adjQ, dtype=np.int64)
    Qa = np.array(Q, dtype=np.int64)
    assert np.abs(A).max() < 2 ** 24 and np.abs(adjQ).max() < 2 ** 20
    t = np.einsum('nxy,yx->n', A, adjQ)
    S = t[:, None, None] * Qa[None] - 2 * detQ * A
    assert np.abs(S).max() < 2 ** 62
    return pd_charpoly_vec(S)


def margins_float(Q, A):
    L = np.linalg.cholesky(np.array(Q, dtype=float))
    Li = np.linalg.inv(L)
    W = Li @ A.astype(float) @ Li.T
    ev = np.linalg.eigvalsh(W)
    tr = ev.sum(axis=1)
    with np.errstate(divide='ignore', invalid='ignore'):
        mg = np.where(tr > 0, (tr - 2 * ev[:, -1]) / np.where(tr > 0, tr, 1), np.inf)
        al = np.where(tr > 0, tr / np.where(ev[:, -1] > 0, ev[:, -1], 1) - 2, np.inf)
    return mg, al


def exact_rank_vec(A):
    """exact rank (0..3) of (N,3,3) integer matrices."""
    N = len(A)
    nz = np.any(A.reshape(N, -1) != 0, axis=1)
    det = det3_obj(A)
    r3 = np.array([x != 0 for x in det], dtype=bool)
    o = A.astype(object)
    r2 = np.zeros(N, dtype=bool)
    for (p, q) in itertools.combinations(range(3), 2):
        for (s, t) in itertools.combinations(range(3), 2):
            mnr = o[:, p, s] * o[:, q, t] - o[:, p, t] * o[:, q, s]
            r2 |= np.array([x != 0 for x in mnr], dtype=bool)
    return np.where(r3, 3, np.where(r2, 2, np.where(nz, 1, 0)))
