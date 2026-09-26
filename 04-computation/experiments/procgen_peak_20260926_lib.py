#!/usr/bin/env python3
"""procgen_peak_20260926_lib -- core exact/float machinery for the peak-discounted provability price.

Setting.  T_q(n) = n/2 (n even), (q n + 1)/2 (n odd), q odd >= 3.  A parity word u in {0,1}^L has
e_j(u) = number of ones among u_0..u_{j-1}, prefix slopes w_j(u) = q^(e_j)/2^j (w_0 = 1).
    Bad_L       = { u : w_j(u) > 1 for every 1 <= j <= L }        (w_j = 1 is impossible for j >= 1)
    w*(u)       = max_{0 <= j <= L-1} w_j(u)                       (the peak slope; edits act at times <= L-1)
    rho_L       = |Bad_L| / 2^L
    rho^peak_L  = 2^-L * sum_{u in Bad_L} 1 / w*(u)
With c = log_q 2 and S_j = e_j - j c we have w_j = q^(S_j), so w* = q^M with M = max_{j<L} S_j.

Layer cake.  F(s) = #{u in Bad_L : M(u) <= s} is a right-continuous step function whose jumps sit at the
breakpoints s = e - j c (0 <= j <= L-1, integer e).  Then (Abel summation, exact)
    sum_{u in Bad_L} q^(-M(u)) = sum_k (F(s_k) - F(s_{k-1})) q^(-s_k),
and F(s) is a lattice-path count in the strip  j c < e_j <= j c + s  (1 <= j <= L-1),  e_L > L c.
All floors floor(m c) are computed exactly by integer comparison q^n <= 2^m.
"""
import math
from fractions import Fraction
import numpy as np


def check(cond, msg):
    """Every printed claim of the runner goes through here; a failed claim aborts the run."""
    if not cond:
        raise SystemExit("CHECK FAILED: " + msg)
    print("  ok:", msg, flush=True)


def H2(x):
    return -x * math.log2(x) - (1 - x) * math.log2(1 - x)


def cq(q):
    return math.log(2) / math.log(q)


def kappa(q):
    """Mogulskii constant of the second-order term: rho^peak = 2^{-(1-H(c))L} exp(-kappa L^{1/3}(1+o(1)))."""
    c = cq(q)
    s2 = c * (1 - c)
    lam = (1 - c) / c
    z = q / max(1.0, lam)
    return 1.5 * (math.pi ** 2 * s2) ** (1 / 3) * math.log(z) ** (2 / 3), z, lam, s2


# ---------------------------------------------------------------------------------------------------
# exact floors
# ---------------------------------------------------------------------------------------------------
def floor_table(q, M):
    """A[k + M] = floor(k * log_q 2) for -M <= k <= M (exact: largest n with q^n <= 2^k; negative k via
    floor(-x) = -floor(x) - 1 for irrational x = |k| c)."""
    vals = [0] * (M + 1)
    n, qn1, p2 = 0, q, 1
    for m in range(M + 1):
        while qn1 <= p2:
            n += 1
            qn1 *= q
        vals[m] = n
        p2 <<= 1
    A = np.empty(2 * M + 1, dtype=np.int64)
    for m in range(M + 1):
        A[M + m] = vals[m]
        if m:
            A[M - m] = -(vals[m] + 1)
    return A


# ---------------------------------------------------------------------------------------------------
# brute force (enumeration of all 2^L words); exact rationals
# ---------------------------------------------------------------------------------------------------
def brute(q, L):
    """Return dict with exact rho (Fraction), rho_peak (Fraction), strata {j: |S_j|} where
    S_j = {u in Bad_L : 2^j <= w*(u) < 2^(j+1)}, peak-weighted strata {j: 2^-L sum_{S_j} 1/w*}."""
    lq, l2 = math.log(q), math.log(2)
    N = 1 << L
    idx = np.arange(N, dtype=np.int64)
    bits = ((idx[:, None] >> np.arange(L)[None, :]) & 1).astype(np.int16)
    e = np.zeros((N, L + 1), dtype=np.int16)
    e[:, 1:] = np.cumsum(bits, axis=1)
    j = np.arange(L + 1)
    # exact badness: q^e_j > 2^j  <=>  e_j >= floor(j c) + 1
    fl = floor_table(q, L + 1)
    lo = np.array([0] + [int(fl[L + 1 + jj]) + 1 for jj in range(1, L + 1)])
    bad = np.all(e[:, 1:] >= lo[None, 1:], axis=1)
    x = e * lq - j[None, :] * l2               # distinct values for distinct (e,j): float argmax is safe at L<=20
    jstar = np.argmax(x[:, :L], axis=1)
    estar = e[np.arange(N), jstar]
    jb, eb = jstar[bad], estar[bad]
    rho = Fraction(int(bad.sum()), N)
    rp = Fraction(0)
    strata, pstrata = {}, {}
    if bad.any():
        pairs, counts = np.unique(np.stack([jb, eb], axis=1), axis=0, return_counts=True)
        for (jj, ee), cnt in zip(pairs.tolist(), counts.tolist()):
            w = Fraction(q ** ee, 2 ** jj)
            rp += Fraction(cnt, 1) / w
            # stratum index floor(log2 w) exactly
            t = 0
            while w >= 2 ** (t + 1):
                t += 1
            strata[t] = strata.get(t, 0) + cnt
            pstrata[t] = pstrata.get(t, Fraction(0)) + Fraction(cnt, N) / w
    return {"rho": rho, "rho_peak": rp / N, "strata": strata, "pstrata": pstrata}


# ---------------------------------------------------------------------------------------------------
# strip DP (vectorised over breakpoints); exact (object ints) or float (probabilities)
# ---------------------------------------------------------------------------------------------------
def strip_counts(q, L, fl, M, jk, ek, exact=False):
    """For each breakpoint (jk[k], ek[k]) with strip height s_k = ek - jk c >= 0 return
    #{u in Bad_L : S_j(u) <= s_k for 0 <= j <= L-1}  (exact=True, python ints)  or that count / 2^L (float).
    Column j (1 <= j <= L-1) allows fl[j]+1 <= e <= ek + fl[j - jk]; column L only needs e >= fl[L]+1."""
    jk = np.asarray(jk, dtype=np.int64)
    ek = np.asarray(ek, dtype=np.int64)
    K = len(jk)

    def D(j):
        return 0 if j == 0 else int(fl[M + j]) + 1

    def W(j):
        return ek + fl[M + j - jk] - D(j) + 1

    Wmax = max(int(W(j).max()) for j in range(L)) if L > 0 else 1
    Wmax = max(Wmax, 1)
    if exact:
        f = np.zeros((K, Wmax + 2), dtype=object)
        f[:, :] = 0
        f[:, 0] = 1
    else:
        f = np.zeros((K, Wmax + 2))
        f[:, 0] = 1.0
    cols = np.arange(Wmax + 2)[None, :]
    for j in range(0, L - 1):
        d = D(j + 1) - D(j)
        new = np.empty_like(f)
        if d == 0:
            new[:, 0] = f[:, 0]
            new[:, 1:] = f[:, 1:] + f[:, :-1]
        else:
            new[:, :-1] = f[:, 1:] + f[:, :-1]
            new[:, -1] = f[:, -1]
        if not exact:
            new *= 0.5
        new[cols >= W(j + 1)[:, None]] = 0
        f = new
    d = D(L) - D(L - 1)
    if exact:
        wts = np.array([2] * (Wmax + 2), dtype=object)
        if d == 1:
            wts[0] = 1
        return [int(v) for v in f.dot(wts)]
    wts = np.full(Wmax + 2, 1.0)
    if d == 1:
        wts[0] = 0.5
    return f @ wts


def rho_exact(q, L, fl=None, M=None):
    """exact |Bad_L| by a 1-D DP (python ints)."""
    if fl is None:
        M = L + 2
        fl = floor_table(q, M)
    f = {0: 1}
    for j in range(L):
        lo = int(fl[M + j + 1]) + 1
        nf = {}
        for e, v in f.items():
            for b in (0, 1):
                if e + b >= lo:
                    nf[e + b] = nf.get(e + b, 0) + v
        f = nf
    return sum(f.values())


def rho_float(q, L, fl, M):
    f = np.zeros(L + 2)
    f[0] = 1.0
    for j in range(L):
        new = 0.5 * f
        new[1:] += 0.5 * f[:-1]
        new[: int(fl[M + j + 1]) + 1] = 0.0
        f = new
    return float(f.sum())


def rho_band_exact(q, L, K_num, K_den, fl=None, M=None, include_L=True):
    """THM-4478 band population: #{u : 1 <= w_j <= K for 0<=j<=L} (include_L) or with the upper
    condition only for j <= L-1 (include_L=False); K = K_num/K_den given as integers; exact."""
    f = {0: 1}
    for j in range(1, L + 1):
        nf = {}
        for e, v in f.items():
            for b in (0, 1):
                e2 = e + b
                if q ** e2 < 2 ** j:
                    continue
                if (j <= L - 1 or include_L) and q ** e2 * K_den > K_num * 2 ** j:
                    continue
                nf[e2] = nf.get(e2, 0) + v
        f = nf
    return sum(f.values())


def breakpoints(q, L, fl, M, smax=None):
    """(j, e) with 0<=j<=L-1, e >= fl[j]+1 (j>=1), e <= j, s = e - j c <= smax; plus (0,0)."""
    c = cq(q)
    out = [(0, 0)]
    for j in range(1, L):
        for e in range(int(fl[M + j]) + 1, j + 1):
            if smax is not None and e - j * c > smax:
                break
            out.append((j, e))
    return out


def rho_peak_exact(q, L):
    """exact rho_L, rho^peak_L (Fractions) and exact strata counts, all breakpoints (small/moderate L)."""
    M = 2 * L + 4
    fl = floor_table(q, M)
    bp = breakpoints(q, L, fl, M)
    bp.sort(key=lambda p: Fraction(q ** p[1], 2 ** p[0]))
    F = strip_counts(q, L, fl, M, [p[0] for p in bp], [p[1] for p in bp], exact=True)
    S = Fraction(0)
    prev = 0
    for (j, e), Fk in zip(bp, F):
        S += Fraction((Fk - prev) * 2 ** j, q ** e)
        prev = Fk
    nbad = rho_exact(q, L)
    assert F[-1] == nbad
    # strata: F at s = t c (virtual breakpoints (-t, 0)), t = 1..T with 2^T > max w*
    T = int(math.ceil((L - 1) * math.log2(q / 2))) + 2
    Fv = strip_counts(q, L, floor_table(q, L + T + 4), L + T + 4, [-t for t in range(1, T + 1)], [0] * T, exact=True)
    strata = {}
    prevv = 0
    for t in range(1, T + 1):
        cnt = Fv[t - 1] - prevv
        if cnt:
            strata[t - 1] = cnt
        prevv = Fv[t - 1]
    assert prevv == nbad
    return {"rho": Fraction(nbad, 2 ** L), "rho_peak": S / 2 ** L, "strata": strata, "nbad": nbad}


def rho_peak_float(q, L, smax, chunk=20000, want_strata=False):
    """float layer cake with truncation at s <= smax; returns dict with rigorous bracket [lo, hi]
    (up to float rounding) for rho^peak, rho_L, #breakpoints, and optionally peak-weighted strata."""
    c = cq(q)
    M = 2 * L + 8 + (int(smax / c) + 4 if want_strata else 0)
    fl = floor_table(q, M)
    bp = breakpoints(q, L, fl, M, smax)
    s = np.array([e - j * c for (j, e) in bp])
    order = np.argsort(s, kind="stable")
    bp = [bp[i] for i in order]
    s = s[order]
    gaps = np.diff(s)
    assert len(gaps) == 0 or gaps.min() > 1e-9, "breakpoint order not certified"
    jk = np.array([p[0] for p in bp])
    ek = np.array([p[1] for p in bp])
    G = np.empty(len(bp))
    for a in range(0, len(bp), chunk):
        G[a:a + chunk] = strip_counts(q, L, fl, M, jk[a:a + chunk], ek[a:a + chunk])
    invw = np.exp(-s * math.log(q))
    rho = rho_float(q, L, fl, M)
    dG = np.diff(np.concatenate([[0.0], G]))
    # Abel form with positive terms only: sum_{k<K} G_k (q^-s_k - q^-s_{k+1}) + G_K q^-s_K
    dinv = invw[:-1] * (-np.expm1(-np.diff(s) * math.log(q)))
    A = float(np.sum(G[:-1] * dinv) + G[-1] * invw[-1])
    out = {"lo": A, "hi": A + max(rho - G[-1], 0.0) * invw[-1], "rho": rho, "nbp": len(bp),
           "G_last": float(G[-1])}
    if want_strata:
        # peak-weighted strata: stratum t collects breakpoints with 2^t <= q^s < 2^(t+1), i.e. t c <= s < (t+1) c
        tk = np.floor(s / c + 1e-12).astype(np.int64)
        ps = {}
        for t in np.unique(tk):
            ps[int(t)] = float(np.sum((dG * invw)[tk == t]))
        out["pstrata"] = ps
        T = int(smax / c)
        Fv = strip_counts(q, L, fl, M, [-t for t in range(1, T + 1)], [0] * T)
        cnt = {}
        prevv = 0.0
        for t in range(1, T + 1):
            cnt[t - 1] = float(Fv[t - 1] - prevv)
            prevv = Fv[t - 1]
        out["strata"] = cnt
    return out


# ---------------------------------------------------------------------------------------------------
# capacity constants
# ---------------------------------------------------------------------------------------------------
def Mstar(q, L, fl=None, M=None):
    """single-scale capacity: sum_{k=0}^{L-1} sum_{e=e_min(k)}^{k} (floor(e/q)+1), e_min(0)=0, e_min(k)=fl[k]+1."""
    if fl is None:
        M = L + 2
        fl = floor_table(q, M)
    tot = 1  # k = 0, e = 0
    for k in range(1, L):
        for e in range(int(fl[M + k]) + 1, k + 1):
            tot += e // q + 1
    return tot


def Mstar_closed_bound(q, L):
    return Fraction(L * (L + 1) * (2 * L - 2 + 3 * q), 6 * q)


def R_L(q, L):
    return sum(k // q + 1 for k in range(L))


def window(q, K_num, K_den=1):
    """floor(log_q K) + 1 exactly for K = K_num/K_den >= 1."""
    t = 0
    while q ** (t + 1) * K_den <= K_num:
        t += 1
    return t + 1


def M_L(q, L, K_num, K_den=1):
    return window(q, K_num, K_den) * R_L(q, L)


# ---------------------------------------------------------------------------------------------------
# T2: Chernoff majorant and block minorant
# ---------------------------------------------------------------------------------------------------
def chernoff_sum_exact(q, L):
    """(q/2) * sum_{e > cL} C(L,e) q^(-e)  (exact Fraction); majorant of rho^peak_L."""
    tot = Fraction(0)
    for e in range(L + 1):
        if q ** e > 2 ** L:
            tot += Fraction(math.comb(L, e), q ** e)
    return Fraction(q, 2) * tot


def chernoff_sum_log2(q, L):
    """log2 of (q/2) sum_{e > cL} C(L,e) q^-e via log-sum-exp (large L); e > cL decided exactly (q^e > 2^L)."""
    p2 = 2 ** L
    terms = [math.lgamma(L + 1) - math.lgamma(e + 1) - math.lgamma(L - e + 1) - e * math.log(q)
             for e in range(L + 1) if q ** e > p2]
    m = max(terms)
    return (math.log(q / 2) + m + math.log(sum(math.exp(t - m) for t in terms))) / math.log(2)


def block_lower_log2(q, L, b):
    """log2 of 2^-L (binom(b,k_b)/b)^t q^-(t+b+r): a minorant of rho^peak_L (THM-4478 section 4 with c_q)."""
    c = cq(q)
    kb = int(math.floor(c * b)) + 1  # ceil(c b), c b irrational
    t, r = divmod(L, b)
    return (-L + t * (math.log2(math.comb(b, kb)) - math.log2(b)) - (t + b + r) * math.log2(q))


# ---------------------------------------------------------------------------------------------------
# T3: elementary two-sided confinement bounds (Lemma C) evaluated numerically
# ---------------------------------------------------------------------------------------------------
def elem_upper_norm(q, L):
    """rigorous majorant of 2^{(1-H)L} rho^peak_L:  C_lam * sum_{m>=0} z^-m (13/14)^floor((L-1)/n_{m+1}),
    n_a = ceil(2 a^2 / sigma^2)."""
    c = cq(q)
    s2 = c * (1 - c)
    lam = (1 - c) / c
    z = q / max(1.0, lam)
    Clam = max(1.0, lam ** (1 - c))
    tot = 0.0
    for m in range(0, 100000):
        a = m + 1
        na = math.ceil(2 * a * a / s2)
        term = z ** (-m) * (13 / 14) ** ((L - 1) // na)
        tot += term
        if z ** (-m) < 1e-300 or (m > 10 and z ** (-m) < tot * 1e-18):
            break
    return Clam * tot


def elem_lower_norm(q, L):
    """rigorous minorant of 2^{(1-H)L} rho^peak_L: max over integers a >= 38 of
    (q/min(1,lam))^-a * c^ceil(a/(2(1-c))) * 32^-ceil(L/floor(a^2/(1152 sigma^2)))  (log-scale, natural log)."""
    c = cq(q)
    s2 = c * (1 - c)
    lam = (1 - c) / c
    best = -float("inf")
    for a in range(38, 4000):
        n = int(a * a // (1152 * s2))
        if n < 1 or n * s2 < 1:
            continue
        r = math.ceil(a / (2 * (1 - c)))
        val = (-a * (math.log(q) - math.log(min(1.0, lam))) + r * math.log(c)
               - math.ceil(L / n) * math.log(32))
        best = max(best, val)
    return best  # natural log
