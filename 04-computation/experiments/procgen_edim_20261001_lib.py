"""procgen_edim_20261001_lib.py -- helpers for the edge multiset dimension lane (2026-10-01).

Definitions (Allikvere, arXiv:2608.09983, Sec. 1-2): for an edge e=uv and a vertex s,
d(e,s)=min(d(u,s),d(v,s)); the edge multiset representation of e w.r.t. S is the multiset
{d(e,s): s in S}, equivalently the histogram H^S_e(r)=#{s in S: d(e,s)=r}, r=0..d-1.
S (nonempty) is edge-multiset resolving iff all edge histograms are pairwise distinct.
Everything here is independent of the C search programs (pure Python / numpy).
"""
import itertools
from collections import Counter
from math import comb, lgamma, log, log2, exp, sqrt
import numpy as np

# ---------------------------------------------------------------- exact verifier (no shortcuts)
def edges(d):
    return [(u, u | (1 << i)) for u in range(1 << d) for i in range(d) if not (u >> i) & 1]
def hd(a, b): return bin(a ^ b).count('1')
def reps(d, S):
    return [tuple(sorted(min(hd(u, s), hd(v, s)) for s in S)) for (u, v) in edges(d)]
def is_resolving(d, S):
    R = reps(d, S); return len(S) > 0 and len(set(R)) == len(R)
def defect(d, S):
    R = reps(d, S); return len(R) - len(set(R))
def colliding_pairs(d, S):
    return sum(c * (c - 1) // 2 for c in Counter(reps(d, S)).values())
def mask_to_set(mask): return [v for v in range(mask.bit_length()) if (mask >> v) & 1]
def set_to_mask(S): return sum(1 << s for s in S)

def is_resolving_np(d, S):
    """vectorised exact check for larger d (histograms as integer rows, compared exactly)"""
    N = 1 << d; S = np.asarray(sorted(set(S)), dtype=np.int64)
    U = np.array([u for u in range(N) for i in range(d) if not (u >> i) & 1], dtype=np.int64)
    I = np.array([i for u in range(N) for i in range(d) if not (u >> i) & 1], dtype=np.int64)
    H = np.zeros((len(U), d), dtype=np.int64)
    pc = np.array([bin(x).count('1') for x in range(N)], dtype=np.int64)
    for s in S:
        dist = pc[(U ^ s) & ~(1 << I)]          # projection lemma: d(e,s)=|u'-s'|
        H[np.arange(len(U)), dist] += 1
    return len(S) > 0 and len(np.unique(H, axis=0)) == len(U)

# ---------------------------------------------------------------- hypercube automorphisms
def aut_perms(n):
    """all 2^n n! automorphisms of Q_n as vertex permutations (x -> pi(x xor t))"""
    P = []
    for pi in itertools.permutations(range(n)):
        for t in range(1 << n):
            P.append([sum(1 << pi[i] for i in range(n) if ((v ^ t) >> i) & 1) for v in range(1 << n)])
    return P
def stabilizer_size(n, S, P=None):
    P = P or aut_perms(n); T = set(S)
    return sum(1 for g in P if set(g[s] for s in S) == T)
def burnside_subset_orbits(n):
    """number of orbits of a-subsets of V(Q_n) under Aut(Q_n), a=0..2^n (cycle-index / Burnside)"""
    N = 1 << n; tot = [0] * (N + 1); P = aut_perms(n)
    for perm in P:
        seen = [False] * N; poly = [1] + [0] * N
        for v in range(N):
            if not seen[v]:
                L = 0; x = v
                while not seen[x]: seen[x] = True; x = perm[x]; L += 1
                new = poly[:]
                for i in range(N + 1 - L):
                    if poly[i]: new[i + L] += poly[i]
                poly = new
        for i in range(N + 1): tot[i] += poly[i]
    assert all(x % len(P) == 0 for x in tot)
    return [x // len(P) for x in tot]

def check_q5_reps(a, path, burnside):
    """rep file is a complete irredundant orbit system: sizes, distinct, each = min of its orbit, count = Burnside"""
    R = np.array([int(l, 16) for l in open(path)], dtype=np.int64)
    assert len(np.unique(R)) == len(R)
    pc = np.zeros(len(R), dtype=np.int64)
    for v in range(32): pc += (R >> v) & 1
    assert (pc == a).all()
    G = np.array(aut_perms(5), dtype=np.int64)
    bits = [((R >> v) & 1) for v in range(32)]
    for g in G:
        img = np.zeros(len(R), dtype=np.int64)
        for v in range(32): img |= bits[v] << int(g[v])
        assert (img >= R).all(), 'rep not minimal in its orbit'
    assert len(R) == burnside[a], (len(R), burnside[a])
    return len(R)

# ---------------------------------------------------------------- leaf counting (numpy DP)
SIG5 = np.array([[1 if not (v >> i) & 1 else -1 for i in range(5)] for v in range(32)])
def _dpB(b):
    W = 2 * b + 1; dp = [np.zeros((W,) * 5, dtype=np.int64) for _ in range(b + 1)]; dp[0][(b,) * 5] = 1
    for v in range(32):
        s = SIG5[v]
        for j in range(b, 0, -1):
            dp[j] += np.roll(dp[j - 1], shift=tuple(int(x) for x in s), axis=(0, 1, 2, 3, 4))
    return dp[b]
_BETA_CACHE = {}
def _betas(path):
    if path not in _BETA_CACHE:
        R = np.array([int(l, 16) for l in open(path)], dtype=np.int64)
        Bt = np.zeros((len(R), 5), dtype=np.int64)
        for v in range(32): Bt += ((R >> v) & 1)[:, None] * SIG5[v][None, :]
        _BETA_CACHE[path] = Counter(map(tuple, Bt.tolist()))
    return _BETA_CACHE[path]
def leafcount(k, a, method, path):
    """number of (A,B) leaves of searchA (method 'A': |beta_i|<=a-b) / searchB ('B': |beta_i|>=a-b)"""
    b = k - a; dl = a - b; dp = _dpB(b); vals = np.arange(-b, b + 1); tot = 0
    betas = _betas(path)
    for beta, cA in betas.items():
        idx = []
        for i in range(5):
            f = beta[i] + vals
            ok = (np.abs(f) <= dl) if method == 'A' else (np.abs(f) >= dl)
            idx.append(np.nonzero(ok)[0])
        tot += cA * int(dp[np.ix_(*idx)].sum())
    return tot

def normal_form(S, rule, canon5):
    """map a subset S of Q6 to the (A,B) normal form of method A (rule='max') or B (rule='min');
    canon5(mask32)->(rep_mask, g) with g an Aut(Q5) vertex permutation sending mask to rep."""
    beta = [sum(1 if not (s >> i) & 1 else -1 for s in S) for i in range(6)]
    ab = [abs(x) for x in beta]
    i0 = ab.index(max(ab)) if rule == 'max' else ab.index(min(ab))
    def sw(x):   # transpose coordinates i0 and 5
        bi, b5 = (x >> i0) & 1, (x >> 5) & 1
        x &= ~((1 << i0) | (1 << 5)); return x | (bi << 5) | (b5 << i0)
    T = [sw(s) for s in S]
    if sum(1 if not (s >> 5) & 1 else -1 for s in T) < 0: T = [s ^ 32 for s in T]
    A = sum(1 << s for s in T if s < 32); Bv = [s - 32 for s in T if s >= 32]
    rep, g = canon5(A)
    return rep, sorted(g[v] for v in Bv), beta

# ---------------------------------------------------------------- tournaments (tiling cube, n=5)
def tiling_tournaments(n=5):
    tiles = [(a, b) for a in range(1, n + 1) for b in range(1, n + 1) if a >= b + 2]
    def tour(t):
        A = set((i + 1, i) for i in range(1, n))
        for j, (a, b) in enumerate(tiles): A.add((b, a) if (t >> j) & 1 else (a, b))
        return A
    return tiles, [tour(t) for t in range(1 << len(tiles))]
def tour_canon(A, n=5):
    return min(tuple(sorted((p[x - 1], p[y - 1]) for (x, y) in A)) for p in itertools.permutations(range(1, n + 1)))
def ham_paths(A, n=5):
    return sum(all((p[i], p[i + 1]) in A for i in range(n - 1)) for p in itertools.permutations(range(1, n + 1)))

# ---------------------------------------------------------------- entropy lower bound
def g_ent(mu): return 0.0 if mu <= 0 else (mu + 1) * log2(mu + 1) - mu * log2(mu)
def entropy_budget(d, m):
    s = 0.0; L = log(m) - (d - 1) * log(2)
    for r in range(d - 1):           # level d-1 dropped: determined by the others (levels sum to m)
        lm = L + lgamma(d) - lgamma(r + 1) - lgamma(d - r)
        if lm < -60: continue
        s += g_ent(exp(lm))
    return s
def entropy_mmin(d):
    need = log2(d) + (d - 1)
    lo, hi = 1, 2
    while entropy_budget(d, hi) < need: hi *= 2
    while lo < hi:
        mid = (lo + hi) // 2
        if entropy_budget(d, mid) >= need: hi = mid
        else: lo = mid + 1
    return lo

# ---------------------------------------------------------------- cells, orbits, union bound
def Nh(n, h, a, b):
    t = 0
    for j in range(0, h + 1):
        if b == a + h - 2 * j and 0 <= a - j <= n - h: t += comb(h, j) * comb(n - h, a - j)
    return t
def cells_parallel(d, h): return [[2 * Nh(d - 1, h, a, b) for b in range(d)] for a in range(d)]
def cells_cross(d, h):
    return [[sum(Nh(d - 2, h, a - x, b - y) for x in (0, 1) for y in (0, 1) if a - x >= 0 and b - y >= 0)
             for b in range(d)] for a in range(d)]
def edge_pair_orbits(d):
    out = []
    for h in range(1, d): out.append(('par', h, d * 2 ** (d - 1) * comb(d - 1, h) // 2, cells_parallel(d, h)))
    for h in range(0, d - 1): out.append(('crs', h, d * (d - 1) * 2 ** (d - 1) * comb(d - 2, h), cells_cross(d, h)))
    return out
def _maxatom_bin(n, q):
    """max_k P(Bin(n,q)=k), attained at k = floor((n+1)q) (or its neighbour); scipy's pmf is accurate for large n"""
    from scipy.stats import binom
    if n == 0: return 1.0
    k = int((n + 1) * q)
    return float(max(binom.pmf(kk, n, q) for kk in (k - 1, k, k + 1) if 0 <= kk <= n))
def _collision(n, q):
    from scipy.stats import binom
    if n == 0: return 1.0
    mu = n * q; s = sqrt(n * q * (1 - q)); lo = max(0, int(mu - 14 * s - 20)); hi = min(n, int(mu + 14 * s + 20))
    k = np.arange(lo, hi + 1); p = binom.pmf(k, n, q); tail = 1 - p.sum()
    return float((p * p).sum() + max(tail, 0.0) + 1e-15)
def atom_diff(n1, n2, q):
    """upper bound on max_y P(X-X'=y), X~Bin(n1,q), X'~Bin(n2,q) independent"""
    from scipy.stats import binom
    if q == 0.5:
        N = n1 + n2
        return exp(lgamma(N + 1) - lgamma(N // 2 + 1) - lgamma(N - N // 2 + 1) - N * log(2)) if N > 0 else 1.0
    if n1 + n2 <= 3000:
        p1 = binom.pmf(np.arange(n1 + 1), n1, q); p2 = binom.pmf(np.arange(n2 + 1), n2, q)
        return float(np.convolve(p1, p2[::-1]).max()) * (1 + 1e-6)
    return min(_maxatom_bin(n1, q), _maxatom_bin(n2, q), sqrt(_collision(n1, q) * _collision(n2, q))) * (1 + 1e-6)
def forest_bound(M, q):
    d = len(M); ed = []
    for a in range(d):
        for b in range(a + 1, d):
            if M[a][b] + M[b][a] > 0: ed.append((-log(atom_diff(M[a][b], M[b][a], q)), a, b))
    ed.sort(reverse=True); par = list(range(d))
    def f(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    W = 0.0
    for w, a, b in ed:
        ra, rb = f(a), f(b)
        if ra != rb: par[ra] = rb; W += w
    return exp(-W)
def union_bound(d, q): return sum(c * forest_bound(M, q) for typ, h, c, M in edge_pair_orbits(d))
def size_bound(d, q):
    """if U<1: smallest M with P(Bin(2^d,q)>M) < 1-U (then edim_m(Q_d)<=M); else None"""
    from scipy.stats import binom
    u = union_bound(d, q)
    if u >= 1: return None, u
    n = 2 ** d; t = (1 - u) * 0.999; M = int(binom.ppf(1 - t, n, q))
    while binom.sf(M, n, q) >= t: M += 1
    return M, u
