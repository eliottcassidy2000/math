"""procgen_edim2_20261001_lib.py -- library for the edim2 lane (2026-10-01).

Edge multiset dimension of hypercubes (Allikvere, arXiv:2608.09983v1): growth rate, uniform estimates,
certified sparse union bounds.  Everything here is independent of the previous lane's code.

Conventions: Q_d, n = d-1.  Pair types: ('par', h), 1 <= h <= n, and ('crs', h), 0 <= h <= n-1.
cells(typ, d, h)[a][b] = n_ab = #{w : d(e,w)=a, d(f,w)=b} for the paper's representative pair.
Rigorous numerics: exact integers / fractions.Fraction, and mpmath.iv (outward-rounded intervals).
"""
from fractions import Fraction as Fr
from math import comb, isqrt, log, sqrt, pi, exp, lgamma
from functools import lru_cache
import itertools

# ----------------------------------------------------------------------------- cells and pair types
def cells(typ, d, h):
    """ordered cell sizes n[a][b] (paper Prop. 14, eqs. (3),(4))"""
    n = d - 1
    M = [[0] * d for _ in range(d)]
    if typ == 'par':
        g = n - h
        for p in range(h + 1):
            for r in range(g + 1):
                M[p + r][h - p + r] += 2 * comb(h, p) * comb(g, r)
    else:
        g = n - 1 - h
        for p in range(h + 1):
            for r in range(g + 1):
                c = comb(h, p) * comb(g, r)
                for eps in (0, 1):
                    for eta in (0, 1):
                        M[eta + p + r][eps + h - p + r] += c
    return M

def type_count(typ, d, h):
    n = d - 1
    return d * 2 ** (d - 2) * comb(n, h) if typ == 'par' else comb(d, 2) * 2 ** d * comb(n - 1, h)

def pair_types(d):
    n = d - 1
    return [('par', h) for h in range(1, n + 1)] + [('crs', h) for h in range(0, n)]

def brute_cells(d, e, f):
    def dist(edge, w):
        u, v = edge
        return min(bin(u ^ w).count('1'), bin(v ^ w).count('1'))
    M = [[0] * d for _ in range(d)]
    for w in range(1 << d):
        M[dist(e, w)][dist(f, w)] += 1
    return M

def representative(typ, d, h):
    """the paper's representative edge pair (bits: coordinate i = bit i); edges in direction d-1 (and d-2)"""
    v = (1 << h) - 1
    e = (0, 1 << (d - 1))
    f = (v, v | (1 << (d - 1))) if typ == 'par' else (v, v | (1 << (d - 2)))
    return e, f

def Nmat(M):
    d = len(M)
    return {(a, b): M[a][b] + M[b][a] for a in range(d) for b in range(a + 1, d) if M[a][b] + M[b][a] > 0}

def is_forest(F, d):
    par = list(range(d))
    def f(x):
        while par[x] != x:
            x = par[x]
        return x
    for a, b in F:
        ra, rb = f(a), f(b)
        if ra == rb:
            return False
        par[ra] = rb
    return True

def kruskal(edges, d, weight):
    """maximum-weight spanning forest (any forest is valid for the forest lemma)"""
    ed = sorted(((weight(e), e) for e in edges), reverse=True)
    par = list(range(d))
    def f(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    F = []
    for _, (a, b) in ed:
        ra, rb = f(a), f(b)
        if ra != rb:
            par[ra] = rb
            F.append((a, b))
    return F

def potential_forest(Ne, d, c):
    """each level a picks the level b with |b-c| < |a-c| maximizing N_ab (a forest by the potential argument)"""
    F = []
    for a in range(d):
        best = None
        for b in range(d):
            if abs(b - c) < abs(a - c):
                key = (min(a, b), max(a, b))
                N = Ne.get(key, 0)
                if N > 0 and (best is None or N > best[0]):
                    best = (N, key)
        if best:
            F.append(best[1])
    return F

# ----------------------------------------------------------------------------- binomial helpers
def b_(M, j):
    return comb(M, j) / 2 ** M if 0 <= j <= M else 0.0

@lru_cache(maxsize=None)
def beta_fr(N):
    """paper's beta(N) = C(N, floor(N/2)) / 2^N, exact"""
    return Fr(comb(N, N // 2), 2 ** N)

# ----------------------------------------------------------------------------- exact max atom
def maxatom_exact(n1, n2, q):
    """max_t P(Bin(n1,q) - Bin(n2,q) = t), exact Fraction (q a Fraction).  O(N^2) big-int ops."""
    a, b = q.numerator, q.denominator
    c = b - a
    # W = U + (n2 - V) ~ Bin(n1, q) * Bin(n2, 1-q); numerators over b^(n1+n2)
    p1 = [comb(n1, k) * a ** k * c ** (n1 - k) for k in range(n1 + 1)]
    p2 = [comb(n2, k) * c ** k * a ** (n2 - k) for k in range(n2 + 1)]
    best = 0
    for t in range(n1 + n2 + 1):
        s = 0
        for k in range(max(0, t - n2), min(n1, t) + 1):
            s += p1[k] * p2[t - k]
        if s > best:
            best = s
    return Fr(best, b ** (n1 + n2))

# ----------------------------------------------------------------------------- interval arithmetic
def iv_ctx(prec=160):
    from mpmath import iv
    iv.prec = prec
    return iv

def iv_frac(iv, x):
    """interval containing the Fraction x"""
    return iv.mpf(x.numerator) / iv.mpf(x.denominator)

def lemmaA_iv(iv, N, q):
    """interval upper bound of Lemma A: (1 + 1/(4x))/sqrt(2 pi x) + exp(-x)/2, x = N q (1-q)"""
    x = iv.mpf(N) * iv_frac(iv, q * (1 - q))
    return (1 + 1 / (4 * x)) / iv.sqrt(2 * iv.pi * x) + iv.exp(-x) / 2

def upper(ivx):
    return ivx.b

# ----------------------------------------------------------------------------- Fourier-Hoelder J values (q = 1/2)
def logJ_float(s):
    """ln of (1/2pi) int_{-pi}^{pi} |cos(t/2)|^s dt = Gamma((s+1)/2) / (sqrt(pi) Gamma(s/2+1)), s >= 0 real"""
    return lgamma((s + 1) / 2) - lgamma(s / 2 + 1) - 0.5 * log(pi)

def J_iv(iv, s_num, s_den=1):
    """interval for (1/2pi) int |cos(t/2)|^s dt with s = s_num/s_den rational >= 0 (exact Gamma ratio)"""
    s = iv.mpf(s_num) / iv.mpf(s_den)
    return iv.gamma((s + 1) / 2) / (iv.sqrt(iv.pi) * iv.gamma(s / 2 + 1))

def J_exact_int(iv, s):
    """(1/2pi) int_{-pi}^{pi} |cos(t/2)|^s dt for integer s >= 0, as an interval:
    s even: C(s, s/2)/2^s ; s odd: (2/pi) (s-1)!!/s!!  (Wallis integrals).  J is decreasing in s."""
    if s % 2 == 0:
        return iv_frac(iv, beta_fr(s))
    num = 1; den = 1
    for i in range(1, s + 1):
        if i % 2 == 0: num *= i
        else: den *= i
    # (s-1)!!/s!! = (2*4*...*(s-1)) / (1*3*...*s)
    return 2 / iv.pi * iv.mpf(num) / iv.mpf(den)

def J_float_int(s):
    return exp(logJ_float(s))

# ----------------------------------------------------------------------------- forest selections
def forests_greedy(Ne, d, k):
    """k pairwise edge-disjoint forests, greedily: max-weight spanning forest of the remaining edges"""
    rest = set(Ne)
    out = []
    for _ in range(k):
        F = kruskal(rest, d, lambda e: Ne[e])
        if not F:
            break
        out.append(F)
        rest -= set(F)
    return out

# ----------------------------------------------------------------------------- certified union bounds
def certified_U_half(d, iv, kmax=3, NEXACT=None):
    """certified upper bound (interval upper end) for the union bound at q = 1/2, per type the minimum of
       (i) the paper's single-forest bound prod beta(N) (exact rationals), and
       (ii) the Fourier-Hoelder bound with k = 2..kmax greedy edge-disjoint forests and weights 1/k:
            prod_j prod_{e in F_j} J(k N_e)^(1/k)   (J exact for integer arguments).
       returns (U_upper, per-type list)."""
    tot = iv.mpf(0)
    per = []
    for typ, h in pair_types(d):
        Ne = Nmat(cells(typ, d, h))
        cnt = type_count(typ, d, h)
        F1 = kruskal(Ne, d, lambda e: Ne[e])
        p1 = Fr(1)
        for e in F1:
            p1 *= beta_fr(Ne[e])
        best = iv_frac(iv, p1); bestk = 1
        for k in range(2, kmax + 1):
            Fs = forests_greedy(Ne, d, k)
            if len(Fs) < k:
                continue
            lp = iv.mpf(0)
            for F in Fs:
                for e in F:
                    lp += iv.log(J_exact_int(iv, k * Ne[e])) / k
            val = iv.exp(lp)
            if val.b < best.b:
                best = val; bestk = k
        tot += cnt * best
        per.append((typ, h, cnt, float(best.b), bestk))
    return tot.b, per

def atom_bound_iv(iv, n1, n2, q, NEXACT=300):
    N = n1 + n2
    if N <= NEXACT:
        return iv_frac(iv, maxatom_exact(n1, n2, q))
    return lemmaA_iv(iv, N, q)

def certified_U_q(d, q, iv, kmax=1, NEXACT=300):
    """certified union bound for density q (Fraction): per type min over
       (i) single max-weight forest with exact atoms (N <= NEXACT) or Lemma A, and
       (ii) Fourier-Hoelder with k = 2..kmax greedy forests, J_q(kN) <= LemmaA(k N q(1-q))."""
    tot = iv.mpf(0)
    for typ, h in pair_types(d):
        M = cells(typ, d, h)
        Ne = Nmat(M)
        cnt = type_count(typ, d, h)
        F1 = kruskal(Ne, d, lambda e: Ne[e])
        val = iv.mpf(1)
        for (a, b) in F1:
            val *= atom_bound_iv(iv, M[a][b], M[b][a], q, NEXACT)
        best = val
        for k in range(2, kmax + 1):
            Fs = forests_greedy(Ne, d, k)
            if len(Fs) < k:
                continue
            lp = iv.mpf(0)
            for F in Fs:
                for e in F:
                    x = iv.mpf(k * Ne[e]) * iv_frac(iv, q * (1 - q))
                    la = (1 + 1 / (4 * x)) / iv.sqrt(2 * iv.pi * x) + iv.exp(-x) / 2
                    if la.b < 1:
                        lp += iv.log(la) / k
            v2 = iv.exp(lp)
            if v2.b < best.b:
                best = v2
        tot += cnt * best
    return tot

def binom_tail_upper(iv, N, q, M):
    """interval upper bound for P(Bin(N,q) > M), via P(X=M+1)/(1-r), r = (N-M-1)q/((M+2)(1-q)) < 1"""
    k = M + 1
    r = Fr(N - k, 1) * q / (Fr(k + 1) * (1 - q))
    assert r < 1
    if k * k * 1000 < N:
        # C(N,k) <= N^k / k!  and  ln k! >= k ln k - k + (1/2) ln(2 pi k)  (Stirling, rigorous lower bound)
        lnC = k * iv.log(iv.mpf(N)) - (k * iv.log(iv.mpf(k)) - k + iv.log(2 * iv.pi * k) / 2)
    else:
        lnC = iv.log(iv.mpf(comb(N, k)))
    lnp = lnC + k * iv.log(iv_frac(iv, q)) + (N - k) * iv.log(iv_frac(iv, 1 - q))
    return iv.exp(lnp) / iv_frac(iv, 1 - r)

def certified_size(d, q, iv, kmax=1, NEXACT=300):
    """returns (M, U_upper, tail_upper) with U + tail < 1 (then edim_m(Q_d) <= M), or (None, U, None)"""
    U = certified_U_q(d, q, iv, kmax, NEXACT)
    if U.b >= 1:
        return None, U.b, None
    N = 2 ** d
    mu = N * q
    M = int(mu) + 1
    # increase M until the tail bound fits
    step = max(1, int(float(mu) ** 0.5 / 8))
    while True:
        if Fr(N - M - 1) * q < Fr(M + 2) * (1 - q):
            t = binom_tail_upper(iv, N, q, M)
            if (U + t).b < 1:
                break
        M += step
    # tighten downward
    while M - 1 > mu:
        if Fr(N - M) * q < Fr(M + 1) * (1 - q):
            t = binom_tail_upper(iv, N, q, M - 1)
            if (U + t).b < 1:
                M -= 1
                continue
        break
    t = binom_tail_upper(iv, N, q, M)
    return M, U.b, t.b

# ----------------------------------------------------------------------------- exact verification of landmark sets
def popcount_table(bits):
    import numpy as np
    t = np.zeros(1 << bits, dtype=np.int8)
    for b in range(bits):
        t += ((np.arange(1 << bits) >> b) & 1).astype(np.int8)
    return t

def is_resolving_exact(d, S, return_pairs=False):
    """exact check that the landmark set S (vertex integers) is edge-multiset resolving in Q_d.
    Projection lemma: for e = {u, u + e_i}, d(e,s) = popcount((u ^ s) & ~(1<<i)).  Histograms are built
    exactly (int32 counts), rows are compared exactly via lexicographic sorting.  Memory ~ E*d*4 bytes."""
    import numpy as np
    N = 1 << d
    S = np.array(sorted(set(int(s) for s in S)), dtype=np.int64)
    assert len(S) > 0 and S.min() >= 0 and S.max() < N
    pc = popcount_table(d).astype(np.int64)
    rows = []
    for i in range(d):
        U = np.array([u for u in range(N) if not (u >> i) & 1], dtype=np.int64)
        mask = ~(1 << i)
        H = np.zeros((len(U), d), dtype=np.int32)
        for s in S:
            np.add.at(H, (np.arange(len(U)), pc[(U ^ s) & mask]), 1) if False else None
            dist = pc[(U ^ s) & mask]
            H[np.arange(len(U)), dist] += 1
        rows.append(H)
    H = np.concatenate(rows, axis=0)
    order = np.lexsort(H.T[::-1])
    Hs = H[order]
    eq = np.all(Hs[1:] == Hs[:-1], axis=1)
    if return_pairs:
        # number of colliding unordered pairs
        from collections import Counter
        c = Counter(map(tuple, H.tolist()))
        return int(sum(v * (v - 1) // 2 for v in c.values()))
    return not bool(eq.any())

# ----------------------------------------------------------------------------- sharper atom bound: e^{-x} I_0(x)
def expI0_iv(iv, x):
    """interval upper bound for e^{-x} I_0(x), x an interval > 0, via the positive power series
    I_0(x) = sum_k (x/2)^(2k)/(k!)^2 with a geometric tail bound (ratio <= 1/4 beyond K >= x)."""
    xb = float(x.b)
    K = int(2 * xb) + 20
    y = (x / 2) ** 2
    term = iv.mpf(1); S = iv.mpf(1)
    for k in range(1, K + 1):
        term = term * y / (k * k)
        S += term
    # tail: next term t_{K+1} <= term * y/(K+1)^2, ratio <= y/(K+2)^2 <= 1/4
    nxt = term * y / ((K + 1) ** 2)
    S += nxt * iv.mpf(4) / 3
    return S * iv.exp(-x)

def atom_fourier_iv(iv, N, q, mult=1):
    """interval upper bound for J_q(mult*N)^(1/mult) via e^{-x}I_0(x), x = mult*N*q(1-q);
    for mult = 1 this bounds max_t P(Bin(n,q)-Bin(N-n,q)=t).  Uses Lemma A when x > 60."""
    x = iv.mpf(mult * N) * iv_frac(iv, q * (1 - q))
    if float(x.a) > 60:
        v = (1 + 1 / (4 * x)) / iv.sqrt(2 * iv.pi * x) + iv.exp(-x) / 2      # Lemma A
    else:
        v = expI0_iv(iv, x)
    if mult == 1:
        return v
    return iv.exp(iv.log(v) / mult)

def certified_U(d, q, iv, kmax=2, NEXACT=150):
    """certified union bound at density q (Fraction).  Per pair type: minimum over
       (i) max-weight forest, atoms exact (N <= NEXACT) or e^{-x}I_0(x) [x = N q(1-q)];
       (ii) Fourier-Hoelder with k = 2..kmax greedy edge-disjoint forests, weights 1/k:
            prod_j prod_{e in F_j} (e^{-kx_e} I_0(k x_e))^(1/k).
       Returns (interval U, list of per-type (typ, h, count, upper bound, k used))."""
    tot = iv.mpf(0)
    per = []
    for typ, h in pair_types(d):
        M = cells(typ, d, h)
        Ne = Nmat(M)
        cnt = type_count(typ, d, h)
        F1 = kruskal(Ne, d, lambda e: Ne[e])
        val = iv.mpf(1)
        for (a, b) in F1:
            N = M[a][b] + M[b][a]
            if N <= NEXACT:
                av = iv_frac(iv, maxatom_exact(M[a][b], M[b][a], q))
            else:
                av = atom_fourier_iv(iv, N, q)
            if av.b < 1:
                val *= av
        best = val; bk = 1
        for k in range(2, kmax + 1):
            Fs = forests_greedy(Ne, d, k)
            if len(Fs) < k:
                continue
            v2 = iv.mpf(1)
            for F in Fs:
                for e in F:
                    av = atom_fourier_iv(iv, Ne[e], q, mult=k)
                    if av.b < 1:
                        v2 *= av
            if v2.b < best.b:
                best = v2; bk = k
        tot += cnt * best
        per.append((typ, h, cnt, float(best.b), bk))
    return tot, per

# ----------------------------------------------------------------------------- Theorem D (closed-form estimate, q = 1/2)
def thmD_funcs(iv):
    """interval versions of the explicit functions of Theorem D"""
    pi_ = iv.pi
    ln2 = iv.log(iv.mpf(2))
    def gam_lo(m):   # >= ln C(m, floor(m/2))
        return m * ln2 - iv.log(pi_ * (m + 2) / 2) / 2
    def gam_hi(m):   # >= ln C(m, floor(m/2)) from above
        return m * ln2 - iv.log(pi_ * (iv.mpf(2 * m + 1) / 2) / 2) / 2
    def ell_lo(m):   # <= ln prod_r C(m, r)
        if m <= 1:
            return iv.mpf(0)
        return iv.mpf(m * m - 1) / 2 - iv.mpf(m + 1) / 2 * iv.log(iv.mpf(m))
    def P_lo(n, h):  # path bound at (h, g = n-1-h)
        g = n - 1 - h
        return ((g + 1) * (iv.log(pi_) + gam_lo(h)) + ell_lo(g)) / 2
    def Zt_lo(n, h):  # zigzag bound with h/2 in place of ceil(h/2)
        g = n - 1 - h
        return iv.mpf(h) / 2 * (iv.log(pi_) + gam_lo(g)) + (ell_lo(h) - gam_hi(h)) / 2
    def P_tan0(n):   # tangent-line lower bound of P_lo at h = 0 (tangent at h_c)
        hc = (n - 1) // 2
        tau0 = iv.log(pi_ * (hc + 2) / 2) - iv.mpf(hc) / (hc + 2)
        return (n * (iv.log(pi_) - tau0 / 2) + ell_lo(n - 1)) / 2
    def Z_tanfar(n):  # tangent-line lower bound of Zt_lo at h = n-1 (tangent at g_* = n-2-h_c)
        hc = (n - 1) // 2
        gs = n - 2 - hc
        sig0 = iv.log(pi_ * (gs + 2) / 2) - iv.mpf(gs) / (gs + 2)
        return iv.mpf(n - 1) / 2 * (iv.log(pi_) - sig0 / 2) + (ell_lo(n - 1) - gam_hi(n - 1)) / 2
    def Phi(n):
        hc = (n - 1) // 2
        vals = [P_lo(n, hc), P_tan0(n), Zt_lo(n, hc + 1), Z_tanfar(n)]
        lo = min(v.a for v in vals)
        return lo, vals
    def Psi(n):
        return iv.mpf(n) / 4 * iv.log(2 * pi_) + (ell_lo(n) - gam_hi(n)) / 4
    return dict(gam_lo=gam_lo, gam_hi=gam_hi, ell_lo=ell_lo, P_lo=P_lo, Zt_lo=Zt_lo, P_tan0=P_tan0,
                Z_tanfar=Z_tanfar, Phi=Phi, Psi=Psi)

def thmD_bound(iv, d):
    """upper end of  d^2 2^(2d-3) e^(-Phi(n)) + d 2^(d-2) e^(-Psi(n))"""
    f = thmD_funcs(iv)
    n = d - 1
    lo, _ = f['Phi'](n)
    t1 = iv.mpf(d * d) * iv.mpf(2) ** (2 * d - 3) * iv.exp(-iv.mpf(lo))
    t2 = iv.mpf(d) * iv.mpf(2) ** (d - 2) * iv.exp(-f['Psi'](n))
    return (t1 + t2).b, t1.b, t2.b

def thmD_crude(n):
    """closed-form lower bounds E1..E4 (for the four endpoint values) and E5 (for Psi), float; used for n >= 19"""
    from math import log, pi
    l2 = log(2)
    E1 = 0.5 * ((n * n - n - 2) / 4 * l2 + (n * n - 2 * n - 3) / 8 - (n + 2) / 4 * log(n * (n + 3) / (8 * pi)))
    E2 = (n * (n - 3) / 8) * l2 + (n * n - 4) / 16 - (n + 1) / 8 * log((n + 2) / (4 * pi)) \
        - (n + 3) / 8 * log((n + 1) / 2) - (n + 1) / 4 * l2
    E3 = 0.5 * ((n / 2) * log(4 * pi / (n + 3)) + (n * n - 2 * n) / 2 - (n / 2) * log(n - 1))
    E4 = (n - 1) / 2 * (log(pi) - 0.5 * log(pi * (n + 2) / 4)) + 0.5 * (((n - 1) ** 2 - 1) / 2 - (n / 2) * log(n - 1) - (n - 1) * l2)
    E5 = (n / 4) * log(2 * pi) + 0.25 * ((n * n - 1) / 2 - (n + 1) / 2 * log(n) - n * l2)
    return E1, E2, E3, E4, E5

def certified_size_v2(d, q, iv, NEXACT=150):
    """(M, U_upper, tail_upper): certified U (certified_U with kmax=2) + binomial tail < 1  =>  edim_m(Q_d) <= M"""
    U, _ = certified_U(d, q, iv, kmax=2, NEXACT=NEXACT)
    if U.b >= 1:
        return None, float(U.b), None
    N = 2 ** d
    mu = N * q
    # smallest M with U + P(Bin(N,q) > M) < 1, searching upward from mu
    M = int(mu) + 1
    while True:
        if Fr(N - M - 1) * q < Fr(M + 2) * (1 - q):
            t = binom_tail_upper(iv, N, q, M)
            if (U + t).b < 1:
                break
        M += max(1, int(float(mu) ** 0.5 / 20))
    while M - 1 > mu and Fr(N - M) * q < Fr(M + 1) * (1 - q):
        t = binom_tail_upper(iv, N, q, M - 1)
        if (U + t).b < 1:
            M -= 1
        else:
            break
    t = binom_tail_upper(iv, N, q, M)
    return M, float(U.b), float(t.b)
