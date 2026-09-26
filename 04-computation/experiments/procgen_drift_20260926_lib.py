#!/usr/bin/env python3
"""
procgen_drift_20260926_lib.py -- shared exact helpers for the drift lane (session collatz-procgen-20260922,
2026-09-26): distance of q n + 1 (q odd; mainly q = 5) to the bounded-lookahead-provable class (i) of the
strategy cube, by sign flips.

Setting (THM-4474): a level-k sign strategy sigma : odd residues mod 2^k -> {+1,-1};
    T_sigma(n) = n/2 (n even),  (q n + sigma(n mod 2^k))/2 (n odd).
Flip set R = odd residues with sigma = -1 (sigma = + elsewhere); q n + 1 itself is R = {}.
Parity graph G_sigma on Z/2^k: s -> the two lifts of T(s) mod 2^(k-1); a cycle with a odd nodes and length p is
expanding iff q^a > 2^p.  Class (i) <=> every cycle of G_sigma is contracting (THM-4474 Theorem A; its proof
uses only that q is odd).

Everything used as a proof step is exact (integers / Fractions).  The C engine (procgen_drift_20260926_engine.c,
loaded with ctypes) is used for speed; every certificate it returns is re-checked here:
  * a class-(i) certificate is an integer potential psi with psi(t) <= psi(s) - w(s) on every edge, w = fd - fn on
    odd and -fn on even nodes, for a fraction F = fn/fd < log_q 2 (Lemma P of the cube-distance note);
  * a failure certificate is an explicit closed walk of G_sigma with q^a > 2^p.
"""
import os
import sys
import ctypes
import hashlib
import subprocess
from fractions import Fraction
from math import comb, gcd, log

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
SCR = os.path.join(ROOT, 'scratch', 'procgen_drift')
ENGINE_SRC = os.path.join(HERE, 'procgen_drift_20260926_engine.c')


def engine_so_path():
    """the compiled engine is named after the source hash, so a rebuild never overwrites a library that another
    running process has loaded"""
    with open(ENGINE_SRC, 'rb') as f:
        tag = hashlib.sha256(f.read()).hexdigest()[:12]
    return os.path.join(SCR, f'engine_{tag}.so')


class CheckFailed(Exception):
    pass


def check(cond, msg):
    if not cond:
        raise CheckFailed("CHECK FAILED: " + msg)


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        h.update(f.read())
    return h.hexdigest()


# ----------------------------------------------------------------------------- thresholds
def cf_of_log2_base(q, terms=60):
    """continued fraction of log_q 2 (high precision float via mpmath); only used to PROPOSE fractions, every
    property used is re-verified exactly by integer comparisons q^a vs 2^p."""
    import mpmath
    mpmath.mp.dps = 120
    x = mpmath.log(2) / mpmath.log(q)
    out = []
    for _ in range(terms):
        a = int(mpmath.floor(x))
        out.append(a)
        fr = x - a
        if fr == 0:
            break
        x = 1 / fr
    return out


def below(q, a, p):
    """a/p < log_q 2  <=>  q^a < 2^p (exact; equality impossible for p >= 1)"""
    return q ** a < 2 ** p


def best_lower_approx(D, q):
    """F = a/p = the largest fraction with p <= D and a/p < log_q 2, with an exact certificate:
    an upper Farey neighbour G = a'/p' > log_q 2 with a'p - a p' = 1 and p + p' > D, p' <= D
    (no fraction with denominator <= D lies strictly between Farey neighbours whose denominators sum > D)."""
    cf = cf_of_log2_base(q)
    # exact Stern-Brocot descent: every mediant is placed by an exact comparison q^a vs 2^p; the continued fraction
    # only bounds the number of descent steps (block i has cf[i] steps); stop when the next mediant exceeds D
    lo = (0, 1)   # 0/1 < c
    hi = (1, 0)   # "infinity"
    i = 0
    while True:
        a_i = cf[i]
        for _ in range(a_i):
            med = (lo[0] + hi[0], lo[1] + hi[1])
            if med[1] > D:
                break
            if below(q, med[0], med[1]):
                lo = med
            else:
                hi = med
        else:
            i += 1
            continue
        break
    a, p = lo
    a2, p2 = hi
    check(below(q, a, p), "lower bound not below log_q 2")
    check(p2 == 0 or not below(q, a2, p2), "upper neighbour not above log_q 2")
    check(a2 * p - a * p2 == 1, "not Farey neighbours")
    check(p <= D and (p2 == 0 or p + p2 > D), "Farey certificate fails")
    return Fraction(a, p)


# ----------------------------------------------------------------------------- basic maps
def step(n, q, sig):
    return n // 2 if n % 2 == 0 else (q * n + sig) // 2


def parity_word(r, L, q):
    """first L parities of the all-plus orbit of r (depends only on r mod 2^L)"""
    w = []
    n = r
    for _ in range(L):
        w.append(n & 1)
        n = step(n, q, 1)
    return w


def targets(k, q, flip):
    """numpy array t0[s] = T_sigma(s) mod 2^(k-1); flip = uint8 array of length 2^k"""
    N = 1 << k
    H = N >> 1
    s = np.arange(N, dtype=np.int64)
    sig = np.where(flip.astype(bool), -1, 1)
    t = np.where(s % 2 == 0, s // 2, (q * s + sig) // 2)
    return (t % H).astype(np.int64)


def flip_array(k, flipset):
    a = np.zeros(1 << k, dtype=np.uint8)
    for r in flipset:
        check(r & 1 and 0 < r < (1 << k), "flip residue must be odd and < 2^k")
        a[r] = 1
    return a


def verify_potential(k, q, flip, psi, F):
    """exact edge-by-edge check: psi(t) <= psi(s) - w(s) for all edges, w = fd-fn (odd), -fn (even)"""
    fn, fd = F.numerator, F.denominator
    N = 1 << k
    H = N >> 1
    psi = np.asarray(psi, dtype=np.int64)
    check(len(psi) == N, "potential has wrong length")
    check(int(psi.max()) < 2 ** 62 and int(psi.min()) > -2 ** 62, "potential out of range")
    t0 = targets(k, q, flip)
    s = np.arange(N)
    w = np.where(s % 2 == 1, fd - fn, -fn).astype(np.int64)
    ok = np.all(psi[t0] <= psi - w) and np.all(psi[t0 + H] <= psi - w)
    return bool(ok)


def verify_expanding_walk(k, q, flip, cyc):
    """cyc is a closed walk of G_sigma; returns (a, p, q^a > 2^p)"""
    N = 1 << k
    H = N >> 1
    t0 = targets(k, q, flip)
    p = len(cyc)
    check(p >= 1, "empty walk")
    for i in range(p):
        s, t = cyc[i], cyc[(i + 1) % p]
        check(0 <= t < N and t % H == t0[s], "not a closed walk of G_sigma")
    a = sum(1 for v in cyc if v & 1)
    return a, p, q ** a > 2 ** p


# ----------------------------------------------------------------------------- engine
class Engine:
    """one live instance at a time: the C engine keeps a single global state, so creating a new Engine invalidates
    the previous one (any later use of a stale instance raises)."""
    _lib = None
    _owner = None

    @classmethod
    def lib(cls):
        if cls._lib is None:
            os.makedirs(SCR, exist_ok=True)
            so = engine_so_path()
            if not os.path.exists(so):
                tmp = so + f'.tmp{os.getpid()}'
                subprocess.run(['cc', '-O2', '-shared', '-fPIC', '-o', tmp, ENGINE_SRC], check=True)
                os.replace(tmp, so)
            lib = ctypes.CDLL(so)
            ll = ctypes.c_longlong
            lib.eng_init.argtypes = [ctypes.c_int, ll, ctypes.c_void_p, ll, ll]
            lib.eng_init.restype = ctypes.c_int
            lib.eng_solve.restype = ctypes.c_int
            lib.eng_toggle.argtypes = [ctypes.c_int]
            lib.eng_toggle.restype = ctypes.c_int
            lib.eng_get_psi.argtypes = [ctypes.c_void_p]
            lib.eng_get_flip.argtypes = [ctypes.c_void_p]
            lib.eng_get_cycle.argtypes = [ctypes.c_void_p]
            lib.eng_get_cycle.restype = ctypes.c_int
            lib.eng_karp.argtypes = [ctypes.POINTER(ll), ctypes.POINTER(ll)]
            lib.eng_kill.argtypes = [ctypes.c_void_p, ctypes.c_int]
            lib.eng_short_cycles.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_void_p,
                                             ctypes.c_int, ctypes.c_void_p]
            lib.eng_short_cycles.restype = ctypes.c_int
            cls._lib = lib
        return cls._lib

    def __init__(self, k, q, flip, F):
        self.k, self.q, self.F = k, q, F
        self.N = 1 << k
        self.flip0 = np.ascontiguousarray(flip, dtype=np.uint8)
        L = self.lib()
        L.eng_init(k, q, self.flip0.ctypes.data, F.numerator, F.denominator)
        Engine._owner = self

    def _live(self):
        check(Engine._owner is self, "stale Engine instance used (the C engine holds one global state)")

    def solve(self):
        self._live()
        return bool(self.lib().eng_solve())

    def toggle(self, r):
        self._live()
        return bool(self.lib().eng_toggle(int(r)))

    def psi(self):
        self._live()
        out = np.zeros(self.N, dtype=np.int64)
        self.lib().eng_get_psi(out.ctypes.data)
        return out

    def flip(self):
        self._live()
        out = np.zeros(self.N, dtype=np.uint8)
        self.lib().eng_get_flip(out.ctypes.data)
        return out

    def cycle(self):
        self._live()
        out = np.zeros(self.N + 1, dtype=np.int32)
        n = self.lib().eng_get_cycle(out.ctypes.data)
        # pointer cycle s -> best[s] is a walk of G_sigma in the forward direction
        return [int(v) for v in out[:n]]

    def kill(self, nodes):
        self._live()
        arr = np.ascontiguousarray(np.array(nodes, dtype=np.int32))
        self.lib().eng_kill(arr.ctypes.data, len(arr))

    def revive_all(self):
        self._live()
        self.lib().eng_revive_all()

    def short_cycles(self, P, cap, needflip=True):
        """simple cycles of length <= P with density > F (F must be the best lower approximation of log_q 2 with
        denominator >= P for 'density > F' to mean 'expanding'), rooted at their minimum node"""
        self._live()
        outcap = cap * P
        out = np.zeros(outcap, dtype=np.int32)
        lens = np.zeros(cap, dtype=np.int32)
        n = self.lib().eng_short_cycles(int(P), int(cap), 1 if needflip else 0, out.ctypes.data, int(outcap),
                                        lens.ctypes.data)
        res = []
        pos = 0
        for i in range(n):
            res.append([int(v) for v in out[pos:pos + lens[i]]])
            pos += lens[i]
        return res

    def karp(self):
        self._live()
        a = ctypes.c_longlong()
        b = ctypes.c_longlong()
        self.lib().eng_karp(ctypes.byref(a), ctypes.byref(b))
        return Fraction(a.value, b.value)


def classify(k, q, flip, F=None):
    """exact decision of class (i) for the strategy with flip array `flip`:
    returns ('I', F, psi) with a verified potential at threshold F (F = best lower approximation of log_q 2 with
    denominator <= 2^k unless given), or ('X', cycle, (a, p)) with a verified expanding closed walk."""
    if F is None:
        F = best_lower_approx(1 << k, q)
    E = Engine(k, q, flip, F)
    if E.solve():
        psi = E.psi()
        check(verify_potential(k, q, flip, psi, F), "engine potential fails the exact check")
        return 'I', F, psi
    cyc = E.cycle()
    a, p, exp_ = verify_expanding_walk(k, q, flip, cyc)
    # the pointer cycle is simple (p <= 2^k) and has density > F, hence > log_q 2 when F is the best lower
    # approximation with denominator <= 2^k
    check(Fraction(a, p) > F, "engine cycle not above threshold")
    check(exp_ or F != best_lower_approx(1 << k, q), "simple cycle above the exact threshold is not expanding")
    return 'X', cyc, (a, p)


def disjoint_cycles(k, q, flip, M, F=None):
    """up to M node-disjoint expanding cycles of G_sigma (each verified; simple cycles above the exact threshold)"""
    if F is None:
        F = best_lower_approx(1 << k, q)
    E = Engine(k, q, flip, F)
    out = []
    for _ in range(M):
        if E.solve():
            break
        cyc = E.cycle()
        a, p, ex = verify_expanding_walk(k, q, flip, cyc)
        check(ex and Fraction(a, p) > F, "harvested cycle not expanding")
        out.append(cyc)
        E.kill(cyc)
    return out


def rho_max(k, q, flip):
    """exact maximum cycle density by Dinkelbach iteration with verified certificates:
    returns (F, cyc, psi): cyc a closed walk of density F, psi a verified potential at threshold F."""
    F = Fraction(0, 1)
    cyc0 = None
    while True:
        E = Engine(k, q, flip, F)
        if E.solve():
            psi = E.psi()
            check(verify_potential(k, q, flip, psi, F), "Dinkelbach potential fails")
            if cyc0 is None:
                # F = 0: need a cycle of density 0 -> the loop at 0 (0 -> 0 is an edge: T(0)=0)
                cyc0 = [0]
            a, p, _ = verify_expanding_walk(k, q, flip, cyc0)
            check(Fraction(a, p) == F, "Dinkelbach witness density mismatch")
            return F, cyc0, psi
        cyc = E.cycle()
        a, p, _ = verify_expanding_walk(k, q, flip, cyc)
        check(Fraction(a, p) > F, "Dinkelbach cycle does not improve")
        F = Fraction(a, p)
        cyc0 = cyc


# ----------------------------------------------------------------------------- necklaces
def necklaces_with_ones(n, j):
    tot = 0
    for t in range(n):
        g = gcd(n, t) if t else n
        if j % (n // g) == 0:
            tot += comb(g, j // (n // g))
    check(tot % n == 0, "Burnside")
    return tot // n


def critical_ones(L, q):
    a = 0
    while not q ** a > 2 ** L:
        a += 1
    return a


def necklace_lower_bound(k, q):
    """N_k(q): binary necklaces of length k with q^(ones) > 2^k (expanding rotation cycles of B(2,k))"""
    a0 = critical_ones(k, q)
    return sum(necklaces_with_ones(k, j) for j in range(a0, k + 1))


# ----------------------------------------------------------------------------- stationary laws
def closed_classes(k, q, flip):
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import connected_components
    N = 1 << k
    H = N >> 1
    t0 = targets(k, q, flip)
    rows = np.repeat(np.arange(N), 2)
    cols = np.stack([t0, t0 + H], axis=1).reshape(-1)
    A = csr_matrix((np.ones(2 * N), (rows, cols)), shape=(N, N))
    ncomp, lab = connected_components(A, directed=True, connection='strong')
    leaving = np.zeros(ncomp, dtype=bool)
    src_lab = lab[rows]
    dst_lab = lab[cols]
    np.logical_or.at(leaving, src_lab[src_lab != dst_lab], True)
    out = []
    for c in range(ncomp):
        if not leaving[c]:
            out.append(np.nonzero(lab == c)[0])
    return out, t0


def stationary(k, q, flip):
    """stationary laws of the uniform-lift chain P_sigma on each closed class (float64)"""
    from scipy.sparse import csr_matrix, identity
    from scipy.sparse.linalg import spsolve
    N = 1 << k
    H = N >> 1
    classes, t0 = closed_classes(k, q, flip)
    res = []
    for nodes in classes:
        n = len(nodes)
        pos = -np.ones(N, dtype=np.int64)
        pos[nodes] = np.arange(n)
        r = np.repeat(np.arange(n), 2)
        c = np.stack([pos[t0[nodes]], pos[t0[nodes] + H]], axis=1).reshape(-1)
        check(np.all(c >= 0), "closed class not closed")
        P = csr_matrix((np.full(2 * n, 0.5), (r, c)), shape=(n, n))
        M = (P.T - identity(n, format='csr')).tolil()
        M[n - 1, :] = np.ones(n)
        b = np.zeros(n)
        b[n - 1] = 1.0
        pi = spsolve(M.tocsr(), b)
        check(np.all(pi > -1e-12) and abs(pi.sum() - 1) < 1e-9, "stationary solve failed")
        res.append((nodes, np.maximum(pi, 0.0)))
    return res


def log_c(q):
    return log(2) / log(q)


# ----------------------------------------------------------------------------- drift-lane analysis helpers
def parity_words_all(k, q):
    """uint8 array W[s, j] = j-th parity of the all-plus orbit of s (s = 0..2^k-1, j < k)"""
    N = 1 << k
    W = np.zeros((N, k), dtype=np.uint8)
    cur = np.arange(N, dtype=np.int64)
    for j in range(k):
        b = cur & 1
        W[:, j] = b
        cur = np.where(b == 1, (q * cur + 1) // 2, cur // 2)
    return W


def gains(k, q):
    """g[s] = ones(Phi_{k-1}((q s + 1)/2)) - ones(Phi_{k-1}((q s - 1)/2)) for odd s (0 for even s)"""
    N = 1 << k
    H = N >> 1
    ones = np.zeros(H, dtype=np.int64)
    cur = np.arange(H, dtype=np.int64)
    for _ in range(k - 1):
        b = cur & 1
        ones += b
        cur = np.where(b == 1, (q * cur + 1) // 2, cur // 2)
    s = np.arange(N, dtype=np.int64)
    g = ones[((q * s + 1) // 2) % H] - ones[((q * s - 1) // 2) % H]
    g[s % 2 == 0] = 0
    return g


def last_letter(k, q):
    """last parity letter Phi_k(t)_{k-1} of every node t (the appended letter)"""
    return parity_words_all(k, q)[:, k - 1]


def stationary_theta(k, q, flip, theta):
    """stationary laws of the theta-chain Q_sigma^theta (the lift whose parity word ends in 1 is taken with
    probability theta) on every closed class (float64); theta = 1/2 is the uniform-lift chain P_sigma"""
    from scipy.sparse import csr_matrix, identity
    from scipy.sparse.linalg import spsolve
    N = 1 << k
    H = N >> 1
    classes, t0 = closed_classes(k, q, flip)
    last = last_letter(k, q)
    res = []
    for nodes in classes:
        n = len(nodes)
        pos = -np.ones(N, dtype=np.int64)
        pos[nodes] = np.arange(n)
        a = t0[nodes]
        b = t0[nodes] + H
        pa = np.where(last[a] == 1, theta, 1 - theta)
        pb = np.where(last[b] == 1, theta, 1 - theta)
        check(np.allclose(pa + pb, 1.0), "the two lifts must carry different last letters")
        r = np.concatenate([np.arange(n), np.arange(n)])
        c = np.concatenate([pos[a], pos[b]])
        check(np.all(c >= 0), "closed class not closed")
        P = csr_matrix((np.concatenate([pa, pb]), (r, c)), shape=(n, n))
        M = (P.T - identity(n, format='csr')).tolil()
        M[n - 1, :] = np.ones(n)
        rhs = np.zeros(n)
        rhs[n - 1] = 1.0
        pi = spsolve(M.tocsr(), rhs)
        check(np.all(pi > -1e-12) and abs(pi.sum() - 1) < 1e-9, "stationary solve failed")
        full = np.zeros(N)
        full[nodes] = np.maximum(pi, 0.0)
        res.append((nodes, full))
    return res


def h2(p):
    from math import log2
    return 0.0 if p <= 0 or p >= 1 else -p * log2(p) - (1 - p) * log2(1 - p)


def partner(s, k, q):
    """s* = s - 2 q^{-1} mod 2^k: the unique odd node whose plus-target pair equals the minus-target pair of s"""
    N = 1 << k
    return (s - 2 * pow(q, -1, N)) % N


def merge_entropy(k, q, flipset, pi):
    """sum over flips s with unflipped partner s* of (pi(s)+pi(s*)) h(pi(s)/(pi(s)+pi(s*))); also returns
    pi(R), pi(R*) (R* = the unflipped partners)"""
    S = set(flipset)
    M = 0.0
    mR = 0.0
    mRs = 0.0
    for s in S:
        sp = partner(s, k, q)
        check(sp & 1 == 1, "partner must be odd")
        if sp in S:
            continue
        a, b = float(pi[s]), float(pi[sp])
        if a + b > 0:
            M += (a + b) * h2(a / (a + b))
        mRs += b
    mR = float(sum(pi[s] for s in S))
    return M, mR, mRs


def p0_root(digits=40):
    """the root p0 in (0, 1/2) of p + h(p) = 1, bracketed rigorously (interval sign change at high precision)"""
    import mpmath
    mpmath.mp.dps = digits
    f = lambda p: p - p * mpmath.log(p, 2) - (1 - p) * mpmath.log(1 - p, 2) - 1
    lo, hi = mpmath.mpf('0.2'), mpmath.mpf('0.25')
    check(f(lo) < 0 < f(hi), "p0 bracket")
    for _ in range(120):
        mid = (lo + hi) / 2
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
    return lo, hi
