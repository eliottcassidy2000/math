#!/usr/bin/env python3
"""
procgen_floor_20260926_lib.py -- exact helpers for the floor lane (session collatz-procgen-20260922,
2026-09-26): the min-max cycle density rho*(q,k) of the q n +- 1 strategy cube, computed as the value
of a two-player mean-payoff game, with certificates for both bounds.

Setting (THM-4474, THM-4481).  Level k, odd q, N = 2^k, H = 2^(k-1).  A sign strategy sigma on the odd
residues mod N gives T(n) = n/2 (n even), (q n + sigma(n mod N))/2 (n odd).  The parity graph G_sigma
has nodes Z/N and edges s -> both lifts of T(s) mod H.  rho_max(sigma) = the largest odd density of a
cycle of G_sigma; rho*(q,k) = min over sigma of rho_max(sigma).

Game (this lane).  Min picks the sign at every odd node (a positional strategy = a sign strategy),
Max picks the lift at every pair P = {P, P+H} (a positional strategy tau : Z/H -> {0,1}).
Certificates at a threshold F = fn/fd, e(s) = fd - fn (s odd), -fn (s even):
  * UPPER  rho*(q,k) <= F : a flip set R and an integer potential psi with
        psi(t) + e(s) <= psi(s)   for every edge s -> t of G_sigma;
    every cycle then has e-weight <= 0, i.e. odd density <= F.
  * LOWER  rho*(q,k) >= F : a nonempty node set W, a lift strategy tau and integers f on W with,
    for every s in W and every option P of s (1 option if s even, both signs if s odd),
        t = P + tau(P) H  in W   and   f(t) <= f(s) + e(s).
    For every sign strategy the (sigma, tau) play from a node of W stays in W and eventually closes
    a cycle of G_sigma; summing f(t) - f(s) <= e(s) around it gives e-weight >= 0, i.e. density >= F.
Both are checked here edge by edge with exact integers (numpy int64 with range checks).
"""
import os
import sys
import ctypes
import hashlib
import subprocess
from fractions import Fraction
from math import log

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
SCR = os.path.join(ROOT, 'scratch', 'procgen_floor')
GAME_SRC = os.path.join(HERE, 'procgen_floor_20260926_game.c')


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


def log_c(q):
    return log(2) / log(q)


def below(q, a, p):
    """a/p < log_q 2  <=>  q^a < 2^p (exact)"""
    return q ** a < 2 ** p


# ----------------------------------------------------------------------------- engine
class Game:
    """ctypes wrapper of the C arena solver; one global state per process (a new Game invalidates the
    previous one)"""
    _lib = None
    _owner = None

    @classmethod
    def lib(cls):
        if cls._lib is None:
            os.makedirs(SCR, exist_ok=True)
            with open(GAME_SRC, 'rb') as f:
                tag = hashlib.sha256(f.read()).hexdigest()[:12]
            so = os.path.join(SCR, f'game_{tag}.so')
            if not os.path.exists(so):
                tmp = so + f'.tmp{os.getpid()}'
                subprocess.run(['cc', '-O2', '-shared', '-fPIC', '-o', tmp, GAME_SRC], check=True)
                os.replace(tmp, so)
            lib = ctypes.CDLL(so)
            ll = ctypes.c_longlong
            lib.g_init.argtypes = [ctypes.c_int, ll]
            lib.g_init.restype = ctypes.c_int
            lib.g_set_allow.argtypes = [ctypes.c_void_p]
            lib.g_set_tau.argtypes = [ctypes.c_void_p]
            lib.g_solve.argtypes = [ctypes.c_int, ll, ll, ll, ll]
            lib.g_solve.restype = ll
            lib.g_get_f.argtypes = [ctypes.c_void_p]
            lib.g_extract_sigma.argtypes = [ctypes.c_void_p]
            lib.g_extract_tau.argtypes = [ctypes.c_void_p]
            lib.g_min_mean_tau.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.POINTER(ll),
                                           ctypes.POINTER(ll)]
            lib.g_exhaustive_minmax.argtypes = [ctypes.POINTER(ll), ctypes.POINTER(ll), ctypes.POINTER(ll)]
            lib.g_exhaustive_minmax.restype = ctypes.c_int
            cls._lib = lib
        return cls._lib

    def __init__(self, k, q, allow=None):
        self.k, self.q = k, q
        self.N = 1 << k
        self.H = self.N >> 1
        L = self.lib()
        L.g_init(k, q)
        if allow is not None:
            self.allow = np.ascontiguousarray(allow, dtype=np.uint8)
            L.g_set_allow(self.allow.ctypes.data)
        else:
            self.allow = None
        Game._owner = self

    def _live(self):
        check(Game._owner is self, "stale Game instance used")

    def freeze_tau(self, tau):
        """freeze Max to the lift tau[P] in {0,1} at every pair (2 = free)"""
        self._live()
        self.tauf = np.ascontiguousarray(tau, dtype=np.uint8)
        check(len(self.tauf) == self.H, "tau length")
        self.lib().g_set_tau(self.tauf.ctypes.data)

    def solve(self, mode, F, cap, maxwork=10 ** 12):
        """mode 0 = Min energy, 1 = Max energy; returns the number of TOP nodes (-1: work limit)"""
        self._live()
        return int(self.lib().g_solve(int(mode), F.numerator, F.denominator, int(cap), int(maxwork)))

    def f(self):
        self._live()
        out = np.zeros(self.N, dtype=np.int64)
        self.lib().g_get_f(out.ctypes.data)
        return out

    def sigma(self):
        self._live()
        out = np.zeros(self.N, dtype=np.uint8)
        self.lib().g_extract_sigma(out.ctypes.data)
        return out

    def tau(self):
        self._live()
        out = np.zeros(self.H, dtype=np.uint8)
        self.lib().g_extract_tau(out.ctypes.data)
        return out

    def min_mean_tau(self, tau, start):
        self._live()
        t = np.ascontiguousarray(tau, dtype=np.uint8)
        a = ctypes.c_longlong()
        b = ctypes.c_longlong()
        self.lib().g_min_mean_tau(t.ctypes.data, int(start), ctypes.byref(a), ctypes.byref(b))
        return Fraction(a.value, b.value) if b.value > 0 and a.value >= 0 else None

    def exhaustive(self):
        self._live()
        a = ctypes.c_longlong()
        b = ctypes.c_longlong()
        c = ctypes.c_longlong()
        r = self.lib().g_exhaustive_minmax(ctypes.byref(a), ctypes.byref(b), ctypes.byref(c))
        check(r == 0, "exhaustive only for k <= 5")
        return Fraction(a.value, b.value), c.value


# ----------------------------------------------------------------------------- arena (numpy, exact)
def arena(k, q):
    """numpy arrays: even target pair te[s] (s even), plus/minus target pairs tp[s], tm[s] (s odd)"""
    N = 1 << k
    H = N >> 1
    s = np.arange(N, dtype=np.int64)
    te = (s // 2) % H
    tp = ((q * s + 1) // 2) % H
    tm = ((q * s - 1) // 2) % H
    return te, tp, tm


def targets_of(k, q, flip):
    """target pair of every node under the flip array (flip[s] = 1 means sign - at odd s)"""
    te, tp, tm = arena(k, q)
    N = 1 << k
    s = np.arange(N)
    fl = np.asarray(flip).astype(bool)
    return np.where(s % 2 == 0, te, np.where(fl, tm, tp))


def weights(k, F):
    N = 1 << k
    s = np.arange(N, dtype=np.int64)
    fn, fd = F.numerator, F.denominator
    return np.where(s % 2 == 1, fd - fn, -fn).astype(np.int64)


def verify_upper(k, q, flip, psi, F, allow=None, chunk=1 << 20):
    """exact: psi(t) + e(s) <= psi(s) on every edge of G_sigma  ==>  rho_max(sigma) <= F   (chunked)"""
    N = 1 << k
    H = N >> 1
    flip = np.asarray(flip, dtype=np.uint8)
    psi = np.asarray(psi)
    check(len(flip) == N and len(psi) == N, "upper certificate: wrong lengths")
    check(not flip[0::2].any(), "upper certificate: flips must sit on odd residues")
    check(int(np.abs(psi).max()) < 2 ** 30, "potential out of range")
    fn, fd = F.numerator, F.denominator
    ok = True
    for lo in range(0, N, chunk):
        s = np.arange(lo, min(N, lo + chunk), dtype=np.int64)
        fl = flip[lo:lo + len(s)].astype(bool)
        if allow is not None:
            a = np.asarray(allow)[lo:lo + len(s)]
            od = (s % 2 == 1)
            check(np.all(np.where(od & fl, a & 2, 1) != 0) and np.all(np.where(od & ~fl, a & 1, 1) != 0),
                  "upper cert uses a forbidden sign")
        t = np.where(s % 2 == 0, (s // 2) % H, np.where(fl, ((q * s - 1) // 2) % H, ((q * s + 1) // 2) % H))
        e = np.where(s % 2 == 1, fd - fn, -fn).astype(np.int64)
        ps = psi[lo:lo + len(s)].astype(np.int64)
        ok &= bool(np.all(psi[t].astype(np.int64) + e <= ps)) and bool(np.all(psi[t + H].astype(np.int64) + e <= ps))
        if not ok:
            return False
    return bool(ok)


def verify_lower(k, q, W, tau, f, F, allow=None, chunk=1 << 20):
    """exact: W nonempty; for s in W and every option P of s (restricted by allow if given):
    t = P + tau[P] H in W and f(t) <= f(s) + e(s)  ==>  every sign strategy (using allowed signs) has a
    cycle of density >= F   (chunked)"""
    N = 1 << k
    H = N >> 1
    W = np.asarray(W, dtype=bool)
    tau = np.asarray(tau)
    f = np.asarray(f)
    check(len(W) == N and len(tau) == H and len(f) == N, "lower certificate: wrong lengths")
    check(bool(W.any()), "lower certificate: W empty")
    check(int(tau.max()) <= 1, "tau must be 0/1")
    check(int(np.abs(f[W]).max()) < 2 ** 30, "f out of range")
    fn, fd = F.numerator, F.denominator
    for lo in range(0, N, chunk):
        s = np.arange(lo, min(N, lo + chunk), dtype=np.int64)
        s = s[W[lo:lo + len(s)]]
        if len(s) == 0:
            continue
        e = np.where(s % 2 == 1, fd - fn, -fn).astype(np.int64)
        fs = f[s].astype(np.int64)
        ev = (s % 2 == 0)
        opts = [(np.where(ev, (s // 2) % H, ((q * s + 1) // 2) % H), 1),
                (np.where(ev, (s // 2) % H, ((q * s - 1) // 2) % H), 2)]
        for P, bit in opts:
            sel = np.ones(len(s), dtype=bool)
            if allow is not None:
                sel = ev | ((np.asarray(allow)[s] & bit) != 0)
            t = P + tau[P].astype(np.int64) * H
            if not (bool(np.all(W[t[sel]])) and bool(np.all(f[t[sel]].astype(np.int64) <= fs[sel] + e[sel]))):
                return False
    return True


# ----------------------------------------------------------------------------- the search
def mediant(a, b):
    return Fraction(a.numerator + b.numerator, a.denominator + b.denominator)


def test_upper(G, F, cap):
    """Min energy at F; returns (flip, psi) if Min wins everywhere, else None"""
    ntop = G.solve(0, F, cap)
    if ntop != 0:
        return None
    f = G.f()
    check(int(f.max()) < 2 ** 30, "potential too large for int32")
    return G.sigma(), f.astype(np.int32)


def test_lower(G, F, cap):
    """Max energy at F; returns (W, tau, f) if Max wins somewhere, else None"""
    ntop = G.solve(1, F, cap)
    if ntop < 0 or ntop == G.N:
        return None
    f = G.f()
    check(int(f.max()) < 2 ** 30, "potential too large for int32")
    W = f >= 0
    f = np.where(W, f, 0).astype(np.int32)
    return W, G.tau(), f


def rho_star(k, q, allow=None, cap=None, log=None, lo=Fraction(0), hi=Fraction(1), tries=4, known_upper=None):
    """exact rho*(q,k) with both certificates, by a Stern-Brocot search on the value.
    lo < hi must be Farey neighbours bracketing the value (0/1, 1/1 always do).  known_upper (e.g. the value
    at level k-1, an upper bound since lifting keeps the map) lets the search skip the tests at mediants
    above it (they are certainly '<').
    Both energy tests run with a moderate cap (caps only shrink the certified regions, so a too-small cap
    can only make the search fail, never produce a false certificate); on failure the caps grow 8x.
    Returns dict(rho, upper=(flip, psi), lower=(W, tau, f), path=[...])."""
    G = Game(k, q, allow)
    N = 1 << k
    check(hi.numerator * lo.denominator - lo.numerator * hi.denominator == 1, "lo, hi must be Farey neighbours")
    lo0, hi0 = lo, hi
    mult = 1
    for attempt in range(tries):
        lo, hi = lo0, hi0
        path = []
        while True:
            m = mediant(lo, hi)
            if m.denominator > N:
                break                                   # search failed: caps too small
            if known_upper is not None and m > known_upper:
                path.append((m, '<*'))
                hi = m
                continue
            fd = m.denominator
            cmax = mult * (8 * fd + (k * fd) // 2 + 8) if cap is None else mult * cap   # Max credits seen: <= 8.3 fd
            cmin = mult * (4 * fd + (k * fd) // 4 + 8) if cap is None else mult * cap   # Min credits seen: <= 3 fd
            c = max(cmax, cmin)
            low = test_lower(G, m, cmax)
            if low is not None:
                up = test_upper(G, m, cmin)
                if up is not None:
                    path.append((m, 'value'))
                    return {'rho': m, 'upper': up, 'lower': low, 'path': path, 'allow': allow, 'cap': c}
                path.append((m, '>'))
                lo = m
            else:
                path.append((m, '<'))
                hi = m
            if log:
                log(f"      k={k} q={q}: {path[-1][0]} {path[-1][1]}")
        mult *= 8
        if log:
            log(f"      k={k} q={q}: search failed with cap multiplier {mult // 8}; retrying")
    check(False, f"rho_star search failed for k={k}, q={q} after {tries} cap increases")


def certify(k, q, res):
    """re-check both certificates of a rho_star result"""
    F = res['rho']
    flip, psi = res['upper']
    W, tau, f = res['lower']
    check(verify_upper(k, q, flip, psi, F, res.get('allow')), f"upper certificate fails (k={k}, q={q})")
    check(verify_lower(k, q, W, tau, f, F, res.get('allow')), f"lower certificate fails (k={k}, q={q})")
    return True


# ----------------------------------------------------------------------------- misc exact tools
def parity_word(x, L, q, sig=1):
    w = []
    n = x
    for _ in range(L):
        w.append(n & 1)
        n = n // 2 if n % 2 == 0 else (q * n + sig) // 2
    return w


def rho_max_exact(k, q, flip):
    """exact max odd density of a cycle of G_sigma (independent pure-python Karp; small k only)"""
    N = 1 << k
    H = N >> 1
    t = targets_of(k, q, flip)
    NEG = -10 ** 9
    D = [[NEG] * N for _ in range(N + 1)]
    for v in range(N):
        D[0][v] = 0
    for j in range(1, N + 1):
        Dp, Dc = D[j - 1], D[j]
        for u in range(N):
            if Dp[u] == NEG:
                continue
            val = Dp[u] + (u & 1)
            a = int(t[u])
            for b in (a, a + H):
                if val > Dc[b]:
                    Dc[b] = val
    best = None
    for v in range(N):
        if D[N][v] == NEG:
            continue
        w = None
        for j in range(N):
            if D[j][v] == NEG:
                continue
            fr = Fraction(D[N][v] - D[j][v], N - j)
            if w is None or fr < w:
                w = fr
        if w is not None and (best is None or w > best):
            best = w
    return best


def value_of_tau(k, q, tau, log=None):
    """exact value of a frozen Max strategy tau (tau[P] in {0,1}): max over start nodes of the least odd
    density of a cycle Min can reach in G^tau; returned with its two certificates (Max-energy lower
    certificate W,f; Min-energy potential showing Min reaches a cycle of density <= value from everywhere)"""
    G = Game(k, q)
    G.freeze_tau(tau)
    N = 1 << k
    lo, hi = Fraction(0), Fraction(1)
    if G.solve(0, lo, 16 * k) == 0:
        return lo, None                      # Min reaches a density-0 cycle from every node
    while True:
        m = mediant(lo, hi)
        check(m.denominator <= N, "value_of_tau search left the Farey range")
        c = 16 * m.denominator + k * m.denominator
        low = test_lower(G, m, c * 8)
        if low is not None:
            ntop = G.solve(0, m, c * 8)
            if ntop == 0:
                return m, low
            lo = m
        else:
            hi = m


def pair_words(k, q):
    """uint8 array Wd[P, j] = j-th letter of the parity word Phi_{k-1}(P) of pair P under the + map"""
    H = 1 << (k - 1)
    Wd = np.zeros((H, k - 1), dtype=np.uint8)
    cur = np.arange(H, dtype=np.int64)
    for j in range(k - 1):
        b = cur & 1
        Wd[:, j] = b
        cur = np.where(b == 1, (q * cur + 1) // 2, cur // 2)
    return Wd


# ----------------------------------------------------------------------------- symmetric certificates
def nu_symmetric_upper(k, q, F, cap=None):
    """Min-energy least fixed point at F (finite everywhere iff some strategy has rho_max <= F); returns
    (flip, f) with f nu-invariant (f(-s) = f(s)) and flip an argmin strategy, or None"""
    G = Game(k, q)
    c = cap if cap is not None else 64 * F.denominator * k
    if G.solve(0, F, c) != 0:
        return None
    f = G.f().astype(np.int64)
    return G.sigma(), f


def is_nu_invariant(k, f):
    N = 1 << k
    s = np.arange(N)
    return bool(np.all(f[(-s) % N] == f))


def nu_symmetric_lower(k, q, F, cap=None):
    """Max-energy least fixed point at F; returns (W, tau, f) with W = -W, f nu-invariant and tau
    nu-equivariant (the tau-lift of -P is minus the tau-lift of P; ties broken consistently), or None"""
    G = Game(k, q)
    N = 1 << k
    H = N >> 1
    c = cap if cap is not None else 64 * F.denominator * k
    ntop = G.solve(1, F, c)
    if ntop == N:
        return None
    f = G.f().astype(np.int64)
    W = f >= 0
    fv = np.where(W, f, 2 ** 40)
    P = np.arange(H)
    tau = (fv[P + H] < fv[P]).astype(np.int64)
    tie = fv[P + H] == fv[P]
    # nu maps lift P + bH to -(P + bH) mod N; equivariant tie-break: on a tie at P, choose b=0 if P <= (-P mod H), and at
    # -P choose the lift -(P) = N - P ... expressed as a lift bit of pair (-P mod H)
    negP = (-P) % H
    for p in np.nonzero(tie)[0]:
        p = int(p)
        m = int(negP[p])
        if p <= m:
            tau[p] = 0
            node = (-(p + 0 * H)) % N          # image of lift (p, 0)
            tau[m] = 1 if node >= H else 0
    return W, tau.astype(np.uint8), np.where(W, f, 0)


def tau_equivariant(k, tau):
    N = 1 << k
    H = N >> 1
    P = np.arange(H)
    lift = P + tau.astype(np.int64) * H
    img = (-lift) % N                         # image node of the tau-lift of P
    negP = (-P) % H
    ok = img == negP + tau[negP].astype(np.int64) * H
    ok[H // 2] = True                         # nu swaps the two lifts of pair H/2 (f ties there automatically)
    return bool(np.all(ok))


# ----------------------------------------------------------------------------- the negative-integer adversary
def U_edges(k, q):
    """U_k: nodes u = 1..H; even u -> u/2; odd u -> rho((q u - 1)/2), rho((q u + 1)/2), rho(m) = the
    representative of m mod H in [1, H].  Returns (u, v1, v2) arrays (v1 = v2 for even u)."""
    H = 1 << (k - 1)
    u = np.arange(1, H + 1, dtype=np.int64)
    def rho(m):
        r = m % H
        return np.where(r == 0, H, r)
    v1 = np.where(u % 2 == 0, u // 2, rho((q * u - 1) // 2))
    v2 = np.where(u % 2 == 0, u // 2, rho((q * u + 1) // 2))
    return u, v1, v2


def U_matches_tau_one(k, q):
    """U_k is the Min graph of the arena under tau = 1 (lift with top bit 1), via node x = N - u"""
    N = 1 << k
    H = N >> 1
    te, tp, tm = arena(k, q)
    u, v1, v2 = U_edges(k, q)
    x = N - u
    succ_p = np.where(x % 2 == 0, te[x], tp[x]) + H
    succ_m = np.where(x % 2 == 0, te[x], tm[x]) + H
    # sign + at x (odd) gives u' = (q u - 1)/2 reduced; sign - gives (q u + 1)/2 reduced
    return bool(np.all(N - succ_p == v1) and np.all(N - succ_m == v2))


def karp_min_density(nodes_succ, par):
    """min odd density of a cycle of a small digraph (succ lists as a 2-column array, parities), exact
    rational, Karp with a virtual source (numpy)"""
    n = len(par)
    INF = 10 ** 9
    D = np.full((n + 1, n), INF, dtype=np.int64)
    D[0, :] = 0
    for j in range(1, n + 1):
        val = D[j - 1] + par
        for c in range(nodes_succ.shape[1]):
            np.minimum.at(D[j], nodes_succ[:, c], val)
    best = None
    for v in range(n):
        if D[n, v] >= INF:
            continue
        w = None
        for j in range(n):
            if D[j, v] >= INF:
                continue
            fr = Fraction(int(D[n, v] - D[j, v]), n - j)
            if w is None or fr > w:
                w = fr
        if best is None or w < best:
            best = w
    return best


# ----------------------------------------------------------------------------- stationary laws (uniform-lift chain)
def closed_classes_of(k, q, flip):
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import connected_components
    N = 1 << k
    H = N >> 1
    t = targets_of(k, q, flip)
    rows = np.repeat(np.arange(N), 2)
    cols = np.stack([t, t + H], axis=1).reshape(-1)
    A = csr_matrix((np.ones(2 * N), (rows, cols)), shape=(N, N))
    nc, lab = connected_components(A, directed=True, connection='strong')
    leaving = np.zeros(nc, dtype=bool)
    np.logical_or.at(leaving, lab[rows][lab[rows] != lab[cols]], True)
    return [np.nonzero(lab == c)[0] for c in range(nc) if not leaving[c]], t


def stationary_odd_upper_certificate(k, q, flip, c, scale=2 ** 24):
    """for every closed class C of the uniform-lift chain of sigma: an integer function h on C with
    scale*[s odd] + (h(t1) + h(t2))/2 <= scale*c + h(s) for all s in C  (exact integer check, c a Fraction);
    averaging against the stationary law gives pi_C(odd) <= c.  Returns (ok, list of float pi_C(odd))."""
    from scipy.sparse import csr_matrix, identity
    from scipy.sparse.linalg import spsolve
    N = 1 << k
    H = N >> 1
    classes, t = closed_classes_of(k, q, flip)
    out = []
    ok = True
    for C in classes:
        n = len(C)
        pos = -np.ones(N, dtype=np.int64)
        pos[C] = np.arange(n)
        a, b = pos[t[C]], pos[t[C] + H]
        check(np.all(a >= 0) and np.all(b >= 0), "closed class not closed")
        P = csr_matrix((np.full(2 * n, 0.5), (np.concatenate([np.arange(n), np.arange(n)]), np.concatenate([a, b]))),
                       shape=(n, n))
        # stationary law
        M = (P.T - identity(n, format='csr')).tolil()
        M[n - 1, :] = np.ones(n)
        rhs = np.zeros(n)
        rhs[n - 1] = 1.0
        pi = spsolve(M.tocsr(), rhs)
        odd = (C % 2 == 1).astype(float)
        g = float(pi @ odd)
        out.append(g)
        # Poisson equation h - P h = odd - g, h[0] = 0 (float), then integer rounding with slack
        M2 = (identity(n, format='csr') - P).tolil()
        M2[0, :] = 0
        M2[0, 0] = 1.0
        r2 = odd - g
        r2[0] = 0.0
        h = spsolve(M2.tocsr(), r2)
        hi = np.round(h * scale).astype(object)
        cs = Fraction(c) * scale
        lhs_ok = True
        for i in range(n):
            lhs = Fraction(int(odd[i]) * scale) + Fraction(int(hi[a[i]]) + int(hi[b[i]]), 2)
            if lhs > cs + int(hi[i]):
                lhs_ok = False
                break
        ok &= lhs_ok
    return ok, out


def critical_cycle_tau1(k, q):
    """value v of the frozen strategy tau = 1 and one cycle of U_k of density exactly v (found inside a strongly
    connected component of the tight subgraph of the Max certificate); returns (v, cycle as u-values, #wraps)"""
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import connected_components
    H = 1 << (k - 1)
    N = 2 * H
    v, low = value_of_tau(k, q, np.ones(H, dtype=np.uint8))
    W, tau, f = low
    f = f.astype(np.int64)
    u, v1, v2 = U_edges(k, q)
    fu = f[N - u]
    ew = np.where(u % 2 == 1, v.denominator - v.numerator, -v.numerator)
    src, dst = [], []
    for vv in (v1, v2):
        t = fu[vv - 1] == fu + ew
        src.append(u[t] - 1)
        dst.append(vv[t] - 1)
    src = np.concatenate(src)
    dst = np.concatenate(dst)
    A = csr_matrix((np.ones(len(src)), (src, dst)), shape=(H, H))
    nc, lab = connected_components(A, directed=True, connection='strong')
    sizes = np.bincount(lab)
    for c in np.argsort(-sizes):
        nodes = np.nonzero(lab == c)[0] + 1
        S = set(int(z) for z in nodes)
        adj = {}
        for a, b in zip(src + 1, dst + 1):
            if int(a) in S and int(b) in S:
                adj.setdefault(int(a), []).append(int(b))
        if not adj:
            continue
        x = int(nodes[0])
        seen = {}
        path = []
        while x not in seen:
            seen[x] = len(path)
            path.append(x)
            x = adj[x][0]
        cyc = path[seen[x]:]
        a = sum(1 for z in cyc if z % 2)
        check(Fraction(a, len(cyc)) == v, "critical cycle density differs from the value")
        # verify it is a closed walk of U_k
        for i, z in enumerate(cyc):
            nz = cyc[(i + 1) % len(cyc)]
            check(nz in (int(v1[z - 1]), int(v2[z - 1])), "critical cycle is not a walk of U_k")
        wraps = sum(1 for i, z in enumerate(cyc) if z % 2 and cyc[(i + 1) % len(cyc)] not in ((q * z - 1) // 2, (q * z + 1) // 2))
        return v, cyc, wraps
    return v, None, None
