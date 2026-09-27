#!/usr/bin/env python3
"""procgen_robin3_20260926_lib -- helpers for the robin3 lane (session collatz-procgen-20260922).

Goal: the Robin inequality N_m(L) <= Gamma * A_(m+c0)(L) for ALL m with a SMALL shift c0 (HYP-9142; THM-4513 has c0 = 15).
Coordinates and objects are those of the robin2 note (V_t = U_t - floor(ct), S_t = V_t - {ct}); the robin2 library is
imported read-only for the Sturmian word delta and the exact counters N_m(L), A_M(L).

Contents
  * interval-arithmetic helpers (mpmath.iv): the larger root eta_+ of p e^{eta d} + (1-p) e^{-eta c} = e^{-kappa},
    the analytic EM criterion (section EM of the note) and the analytic top-hit bound TB-D, both at a threshold m with
    monotone majorants for larger m;
  * float-rigorous dynamic programs (float64, nonnegative sums, exact halving; relative error <= 1.01 (steps+size) u,
    u = 2^-53): the EM kernel comparison over all landing Sturmian factors, and the decoupled top-hit bound TB-B';
  * exact small-case utilities for sanity checks (enumeration of bridges, exact likelihood ratios).
"""
import math
import os
import sys
from fractions import Fraction

import numpy as np
from mpmath import iv

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import procgen_robin2_20260926_lib as R2  # noqa: E402  (read-only import)

C = R2.C
D = R2.D
MU = R2.MU
U_ROUND = 2.0 ** -53

iv.dps = 40
I_C = iv.log(2) / iv.log(3)
I_D = 1 - I_C
I_MU = I_C - iv.mpf(1) / 2
ROBBINS = iv.mpf('0.6753')          # <= sqrt(2/pi) e^{-1/6} = 0.675395...


def ivq(x):
    """exact rational/int -> point interval"""
    x = Fraction(x)
    return iv.mpf(x.numerator) / x.denominator


def hi(x):
    return float(x.b)


def lo(x):
    return float(x.a)


# ------------------------------------------------------------------------------------------------ parameters
def n1_of(m, v0):
    """EM bridge length n1(m) = ceil((m+1)/v0), v0 a Fraction"""
    v0 = Fraction(v0)
    q = Fraction(m + 1) / v0
    return -((-q.numerator) // q.denominator)


def k_b(m):
    """robin2's k_b(m): largest k >= 0 with m - 1 - mu k >= sqrt(k)"""
    return R2.k_b(m)


def k_b_certified(m):
    """k_b(m) recomputed with interval arithmetic at the boundary (guards against float rounding)"""
    k = R2.k_b(m)
    ok_k = (k == 0) or ((m - 1 - I_MU * k - iv.sqrt(k)).a >= 0)
    ok_k1 = (m - 1 - I_MU * (k + 1) - iv.sqrt(k + 1)).b < 0
    if not (ok_k and ok_k1):
        raise AssertionError("k_b certification failed at m=%d" % m)
    return k


# ------------------------------------------------------------------------------------------------ eta_+
def phi_iv(p, h):
    return p * iv.exp(h * I_D) + (1 - p) * iv.exp(-h * I_C)


def eta_plus_float(p, kappa):
    """float: larger root of p e^{h d} + (1-p) e^{-h c} = e^{-kappa} (0 if the minimum exceeds the target)"""
    target = math.exp(-kappa)
    f = lambda h: p * math.exp(h * D) + (1 - p) * math.exp(-h * C)
    hm = math.log((1 - p) * C / (p * D))
    if hm <= 0 or f(hm) > target:
        return 0.0
    lo_, hi_ = hm, 20.0
    for _ in range(200):
        mid = 0.5 * (lo_ + hi_)
        if f(mid) <= target:
            lo_ = mid
        else:
            hi_ = mid
    return lo_


def eta_lower(p_hi, kappa):
    """certified lower bound eta_lo for eta_+(p, kappa), valid for every p <= p_hi (p_hi: interval or number;
    phi_p(h) is increasing in p for h > 0, and {phi_p <= e^{-kappa}} is an interval containing eta_+ as right end).
    Certificate: phi_{p_hi}(eta_lo) <= e^{-kappa} in interval arithmetic. Returns (eta_lo as interval point, float)."""
    P = iv.mpf(iv.mpf(p_hi).b)
    K = iv.mpf(kappa) if not isinstance(kappa, Fraction) else ivq(kappa)
    e0 = eta_plus_float(float(P.b), float(K.a))
    if e0 <= 0:
        raise AssertionError("eta_+ not positive")
    for shrink in (1e-12, 1e-9, 1e-6, 1e-4):
        eta = iv.mpf(repr(e0 * (1 - shrink)))
        eta = iv.mpf(eta.a)
        if phi_iv(P, eta).b <= iv.exp(-K).a:
            return eta, float(eta.a)
    raise AssertionError("could not certify eta lower bound")


# ------------------------------------------------------------------------------------------------ analytic EM
def em_analytic(mE, c0, v0):
    """Certified upper bound for the EM criterion at every m >= mE (bridge length n1(m) = ceil((m+1)/v0)):
        (n-e)/(e+1) + P_e(TH) <= Ef + Top + Late,
        Ef   = (d+v0)/(c-v0)                                   (uniform in m),
        Top  = exp(1/(6 n)) exp(-eta_lo * y_lo)                (Lemma LR + Ville, main part t <= n/2),
        Late = exp(-v*^2 n - 4 v* y_lo) / (1 - exp(-4 v*^2))   (Hoeffding, t > n/2),
    with n = n1(mE), v* = v0 (mE-1-c)/(mE+1+v0) <= chord speed, y_lo = c0 + 2c - 1 <= top gap,
    eta_lo <= eta_+(c - v*, 1/n). Every term is nonincreasing in m (see note), so the bound at mE covers m >= mE.
    Also checks the side conditions of Lemma LR at mE (they only improve with m)."""
    v0q = Fraction(v0)
    n = n1_of(mE, v0q)
    V0 = ivq(v0q)
    Ef = (I_D + V0) / (I_C - V0)
    vstar = V0 * (mE - 1 - I_C) / (mE + 1 + V0)
    ylo = c0 + 2 * I_C - 1
    yhi = c0 + I_C
    eta_i, eta_f = eta_lower(I_C - vstar, Fraction(1, n))
    top = iv.exp(iv.mpf(1) / (6 * n)) * iv.exp(-eta_i * ylo)
    late = iv.exp(-vstar ** 2 * n - 4 * vstar * ylo) / (1 - iv.exp(-4 * vstar ** 2))
    total = Ef + top + late
    # side conditions of Lemma LR (for every m >= mE): n >= 64; p in [0.4, 0.6]; w <= 2(y+1)/n + v <= 1/4; y >= 2;
    # e - i >= 1 and (n-e) - (t-i) >= 1 for t <= n/2 (n (c/2 - v0) >= y + 2)
    p_lo = I_C - V0
    p_hi = I_C - vstar
    w_hi = 2 * (yhi + 1) / n + V0
    side = (n >= 64 and p_lo.a >= 0.4 and p_hi.b <= 0.6 and w_hi.b <= 0.25 and ylo.a >= 2
            and (n * (I_C / 2 - V0)).a >= (yhi + 2).b)
    return dict(n=n, Ef=Ef, vstar=vstar, eta=eta_f, top=top, late=late, total=total, side=side,
                ok=(total.b < 1) and side)


# ------------------------------------------------------------------------------------------------ analytic TB (regime D)
def hoeffding_series(z, v, smax=None):
    """certified upper bound (interval) for sum_{s>=1} exp(-2 (z + v s)^2 / s), z, v > 0 intervals.
    Partial sum to S plus the tail bound sum_{s>S} exp(-4 z v - 2 v^2 s) = e^{-4zv} e^{-2v^2 (S+1)}/(1-e^{-2v^2})."""
    zf, vf = lo(z), lo(v)
    S = int(4 * zf / vf + 40.0 / (vf * vf)) + 10 if smax is None else smax
    tot = iv.mpf(0)
    # vectorised float evaluation of the partial sum with an explicit relative safety factor
    s = np.arange(1, S + 1, dtype=np.float64)
    terms = np.exp(-2.0 * (zf + vf * s) ** 2 / s)
    part = float(np.sum(terms))
    # each term is an upper bound up to relative error 1e-13 (exp and arithmetic in float); inflate
    tot = iv.mpf(part) * (1 + iv.mpf('1e-10')) + iv.mpf(S) * iv.mpf('1e-300')
    tail = iv.exp(-4 * z * v - 2 * v ** 2 * (S + 1)) / (1 - iv.exp(-2 * v ** 2))
    return tot + tail


def tb_regime_D(mX, c0, v0, vmin, vsmall):
    """Certified upper bound q_D for q = P(TH'|BOK_k), valid for every m >= mX and k_b(m) < k <= k1(m) = n1(m)+1.
    q <= max(Q_L, Q_H) + R_B  (endpoint classes L: v_e >= vmin; H: vsmall <= v_e < vmin; B: v_e < vsmall)."""
    v0q, vminq, vsq = Fraction(v0), Fraction(vmin), Fraction(vsmall)
    V0, VMIN, VS = ivq(v0q), ivq(vminq), ivq(vsq)
    kb = k_b_certified(mX)
    kmin = kb + 1
    k1 = n1_of(mX, v0q) + 1
    y = iv.mpf(c0)
    pref = iv.exp(iv.mpf(1) / (6 * kmin))

    def top(v):
        eta_i, eta_f = eta_lower(I_C - v, Fraction(1, kmin))
        main = pref * iv.exp(-eta_i * y)
        late = iv.exp(-v ** 2 * kmin - 4 * v * y) / (1 - iv.exp(-4 * v ** 2))
        return main + late, eta_f

    topL, etaL = top(VMIN)
    QL = topL * I_C / VMIN
    topH, etaH = top(VS)
    # z_* = m - 1 - vmin ((m+1)/v0 + 2) <= smallest class-H endpoint height (k <= k1 <= (m+1)/v0 + 2)
    zstar = (mX - 1) - VMIN * ((mX + 1) / V0 + 2)
    dip = hoeffding_series(zstar, VS)
    if zstar.a <= 0 or dip.b >= 1:
        QH = iv.mpf(math.inf)
    else:
        QH = topH / (1 - dip)
    # class B ratio, worst over k in (k_b, k1] and m >= mX. With K1(m) = (m+1)/v0 + 2 >= k1(m), Kb(m) = k_b(m)+1:
    #   R_B <= rho(m) = (sqrt(K1)/0.6753) (c K1/(m-2)) exp(-Kb [KL_lo - KL_hi(mX)]),
    # KL_lo = 2 (c - vsmall - 1/2)^2 <= KL(c - vsmall || 1/2) (Pinsker), KL_hi = 2 y^2/(1-4y^2) >= KL(1/2+y || 1/2),
    # y_max = max(mu - (m-2)/K1(m), (1+sqrt(Kb))/Kb) is nonincreasing in m. rho(m+1)/rho(m) <= (1+1/(m+1))^1.5 e^{-Delta}
    # since k_b(m+1) >= k_b(m)+1; so rho is decreasing once 1.5/(mX+1) < Delta := KL_lo - KL_hi(mX) (checked).
    K1 = (mX + 1) / V0 + 2
    ymax = iv.mpf(max(hi(I_MU - (mX - 2) / K1), hi((1 + iv.sqrt(kmin)) / kmin)))
    kl_hi = 2 * ymax ** 2 / (1 - 4 * ymax ** 2)
    kl_lo = 2 * (I_C - VS - iv.mpf(1) / 2) ** 2
    Delta = kl_lo - kl_hi
    RB = iv.sqrt(K1) / ROBBINS * (I_C * K1 / (mX - 2)) * iv.exp(-kmin * Delta)
    rb_monotone = (iv.mpf(3) / (2 * (mX + 1))).b < Delta.a
    qD = iv.mpf(max(hi(QL), hi(QH))) + RB
    # side conditions: e* in class L ((m-2)/k1 >= vmin); Lemma LR conditions for classes L, H at every k >= kmin:
    # p = c - v_e in [0.4, 0.6] (v_e <= m/(k_b+1) <= mu + (1+sqrt(kmin))/kmin); w <= 2(c0+2)/kmin + vmax <= 1/4.
    vmax = I_MU + (1 + iv.sqrt(kmin)) / kmin
    side = ((iv.mpf(mX - 2) / k1).a >= VMIN.b and (I_C - vmax).a >= 0.4 and (I_C - VS).b <= 0.6
            and (2 * (c0 + 2) / iv.mpf(kmin) + vmax).b <= 0.25 and kmin >= 64 and zstar.a > 0
            and (kmin * (I_C / 2 - vmax)).a >= c0 + 3 and (VS.b < VMIN.a) and (VMIN.b < V0.a)
            and dip.b < 1 and rb_monotone and ((mX - 2) / K1).a >= VMIN.b)
    return dict(kb=kb, k1=k1, QL=QL, QH=QH, RB=RB, dip=dip, zstar=zstar, etaL=etaL, etaH=etaH, qD=qD, side=side,
                Delta=Delta,
                ok=side and qD.b < 1)


def tb_regime_A(c0):
    """robin2 regime A: k <= k_b(m): P(BOK) >= 3/4 (Kolmogorov), P(TH') <= 3^{-c0} (Ville): q <= (4/3) 3^{-c0}"""
    return Fraction(4, 3) / 3 ** c0


# ------------------------------------------------------------------------------------------------ Sturmian factors
def factor_set(dl_bytes, n, scan):
    seen = {}
    for s in range(scan):
        w = dl_bytes[s:s + n]
        if w not in seen:
            seen[w] = s
    return seen


def landing_words(dl_bytes, n, scan):
    """all factors w = delta_s..delta_{s+n-1} with delta_{s-2} delta_{s-1} = 0 1; complete iff all n+3 factors of
    length n+2 were found (Morse-Hedlund)"""
    allf = factor_set(dl_bytes, n + 2, scan)
    complete = (len(allf) == n + 3)
    out = {}
    for w, s in allf.items():
        if w[0] == 0 and w[1] == 1:
            out[w[2:]] = s + 2
    return out, complete


# ------------------------------------------------------------------------------------------------ float-rigorous DPs
def em_kernel_check(words, m, M):
    """For each word (delta-factor of length n): Pm = P_(s,s+n)(m,1)/2^n, Pm1 = P_(s,s+n)(m-1,1)/2^n (V in [1,M] at
    all n steps). float64 with exact halving; relative error of each entry <= 1.01 n u. Returns arrays (Pm, Pm1)."""
    W = np.array([list(w) for w in words], dtype=np.int8)
    F, n = W.shape
    A = np.zeros((F, M + 3))
    B = np.zeros((F, M + 3))
    A[:, m - 1] = 1.0
    B[:, m] = 1.0
    for i in range(n):
        d1 = (W[:, i] == 1)[:, None]
        for X in (A, B):
            dn = np.zeros_like(X)
            up = np.zeros_like(X)
            dn[:, :-1] = X[:, 1:]
            up[:, 1:] = X[:, :-1]
            X[:] = 0.5 * np.where(d1, dn + X, X + up)
            X[:, 0] = 0.0
            X[:, M + 1:] = 0.0
    return B[:, 1].copy(), A[:, 1].copy()


def survival_lower(h, K, dl, cap_extra=80):
    """beta_h(k), k = 0..K: survival probability (V >= 1 at times 1..k) from V = h at time 0 (S-height h), with an
    extra kill above h + cap_extra (a valid LOWER bound). float64, relative error <= 1.01 (K + size) u."""
    top = h + cap_extra
    f = np.zeros(top + 3)
    f[h] = 1.0
    out = np.empty(K + 1)
    out[0] = 1.0
    for t in range(K):
        g = np.zeros(top + 3)
        if dl[t] == 1:
            g[:-1] += f[1:]
            g += f
        else:
            g += f
            g[1:] += f[:-1]
        g *= 0.5
        g[0] = 0.0
        g[top + 1:] = 0.0
        f = g
        out[t + 1] = f.sum()
    return out


def survival_upper(h, K, dl, cap_extra=80):
    """betabar_h(r), r = 0..K: survival probability from V = h at time 0, with mass above h + cap_extra lumped into a
    state counted as surviving forever (a valid UPPER bound)."""
    top = h + cap_extra
    f = np.zeros(top + 3)
    f[h] = 1.0
    safe = 0.0
    out = np.empty(K + 1)
    out[0] = 1.0
    for t in range(K):
        g = np.zeros(top + 3)
        if dl[t] == 1:
            g[:-1] += f[1:]
            g += f
        else:
            g += f
            g[1:] += f[:-1]
        g *= 0.5
        g[0] = 0.0
        safe += g[top + 1] + g[top + 2]
        g[top + 1:] = 0.0
        f = g
        out[t + 1] = f.sum() + safe
    return out


def tb_Bprime(m, c0, k1, dl, thetas):
    """Regime B' (decoupled top-hit bound), valid for every start time u >= 1 and 1 <= k <= k1:
        q <= e^{-th c0} g^k max_{1<=r<=k-1} [g^{-r} betabar_{M+1}(r)] / beta_{m-1}(k),   g = (e^{th d} + e^{-th c})/2.
    Returns (q_bound as float with safety factor, best theta). Float error: survival values carry relative error
    <= 1.01 (K+size) u < 1e-10; exp/log in float carry < 1e-13; the returned value is inflated by (1 + 1e-8)."""
    M = m + c0
    beta = survival_lower(m - 1, k1, dl)
    bbar = survival_upper(M + 1, k1, dl)
    if np.any(beta[1:] <= 0):
        raise AssertionError("zero survival")
    lb = np.log(beta)
    lbb = np.log(bbar)
    r = np.arange(k1 + 1, dtype=np.float64)
    best = (math.inf, None)
    for th in thetas:
        g = 0.5 * (math.exp(th * D) + math.exp(-th * C))
        lg = math.log(g)
        vals = lbb - r * lg
        vals[0] = -math.inf
        Bk = np.concatenate(([-math.inf], np.maximum.accumulate(vals[:-1])))   # Bk[k] = max_{1<=r<=k-1}
        crit = -th * c0 + np.max(r[2:] * lg + Bk[2:] - lb[2:])
        if crit < best[0]:
            best = (crit, th)
    return math.exp(best[0]) * (1 + 1e-8), best[1]


# ------------------------------------------------------------------------------------------------ exact small cases
def exact_LR(n, e, t, i):
    """LR_t(i) = C(n-t, e-i) / (C(n,e) p^i q^(t-i)), p = e/n (exact Fraction)"""
    if i < 0 or i > t or e - i < 0 or e - i > n - t:
        return Fraction(0)
    p = Fraction(e, n)
    q = 1 - p
    return Fraction(math.comb(n - t, e - i), math.comb(n, e)) / (p ** i * q ** (t - i))


def vpath(x0, bits, word):
    """V-path from site x0 driven by letters bits over the delta-word word"""
    out = [x0]
    x = x0
    for b, d in zip(bits, word):
        x = x + b - d
        out.append(x)
    return out


# ================================================================================================ shift-2 extension
# (Lemma LR-, Lemma TOP-LOW, Proposition EM**, room-weighted regime D'; see note section 2.10-2.12)
_FLOORS = None


def floors_cached(T):
    global _FLOORS
    if _FLOORS is None or len(_FLOORS) < T + 1:
        _FLOORS = R2.floors(max(T, 10 ** 5))
    return _FLOORS


def iid_hit_cdf_lower(p, c0, jshift, t0):
    """Lower bound for P_p(tau <= t), t = 0..t0, where tau = first t >= 1 with U_t > y' + c t and
    y' = c0 + c (jshift + 1) (the supremum of the top distance after jshift leading zeros; the CDF is nonincreasing in
    y'). U_t = number of ones among t i.i.d. Bernoulli(p) letters. Exact boundary: U_t > c0 + c(j+1+t) iff
    U_t >= c0 + floor(c(j+1+t)) + 1. Float DP; each value is a sum of products of nonnegative numbers; relative error
    < 3 t0 u; returned values are multiplied by (1 - 1e-9)."""
    F = floors_cached(jshift + t0 + 2)
    q = 1.0 - p
    dist = np.zeros(t0 + 2)
    dist[0] = 1.0
    cdf = np.zeros(t0 + 1)
    acc = 0.0
    for t in range(1, t0 + 1):
        nd = np.zeros(t0 + 2)
        nd[:] = dist * q
        nd[1:] += dist[:-1] * p
        bnd = c0 + F[jshift + 1 + t] + 1        # U_t >= bnd means the top is hit at time t
        if bnd <= t0 + 1:
            acc += float(nd[bnd:].sum())
            nd[bnd:] = 0.0
        dist = nd
        cdf[t] = acc
    return cdf * (1 - 1e-9)


def lr_lower_log(nprime, t, delta_max, p_lo, p_hi):
    """Lemma LR-: lower bound for log LR_t(i) of a bridge (n', e), p' = e/n' in [p_lo, p_hi], at a hitting
    configuration with 0 < delta <= delta_max, i.e.
       log LR >= -(delta^2/(2(n'-t))) (1/(p'-w) + 1/q') - w (p'-q')^+/(2p'q') - 1/(12(e-i)) - 1/(12(f'-j')),
    worst case over the parameter box. Returns -inf if a side condition fails."""
    w = delta_max / (nprime - t)
    if w >= p_lo - 0.02:
        return -math.inf
    q_lo = 1.0 - p_hi
    pq_min = min(p_lo * (1 - p_lo), p_hi * (1 - p_hi))
    e_minus_i = p_lo * nprime - (delta_max + p_hi * t)     # e - i = p'n' - p't - delta
    f_minus_j = q_lo * nprime - (1 - p_lo) * t             # f' - j' = q'(n'-t) + delta > q'(n'-t)
    if e_minus_i < 1 or f_minus_j < 1:
        return -math.inf
    val = -(delta_max ** 2 / (2.0 * (nprime - t))) * (1.0 / (p_lo - w) + 1.0 / q_lo)
    val -= w * max(0.0, 2 * p_hi - 1) / (2 * pq_min)
    val -= 1.0 / (12 * e_minus_i) + 1.0 / (12 * f_minus_j)
    return val - 1e-12


def em_refined(mE, c0, v0, J, t0):
    """Proposition EM**: certified-by-margins upper bound on E_T[f] for every landing bridge at every m >= mE:
         [ sum_{j=1}^J (d+v0)^j (1 - h_j_lo) + (d+v0)^(J+1)/(1-d-v0) ] / (1 - h0_up).
    h0_up: Lemma TOP (interval arithmetic, via em_analytic's Top + Late with y_lo = c0+2c-1, valid for y >= 2).
    h_j_lo: Lemma TOP-LOW with the iid CDF at (p_lo = c - v0, y'_hi = c0 + c(j+1)) and the Lemma LR- bound at
    n' = n - j, p' in [c - v0, (c - v*)(1 + J/(n-J))], delta <= c0 + c(j+1) + v0 t + 1 (running minimum in t)."""
    r = em_analytic(mE, c0, v0)
    n = r['n']
    h0 = r['top'] + r['late']
    v0f = float(Fraction(v0))
    vstar = lo(r['vstar'])
    p_lo = C - v0f
    p_hi = (C - vstar) * (1 + J / (n - J))
    q0 = D + v0f
    num = 0.0
    hs = []
    for j in range(1, J + 1):
        cdf = iid_hit_cdf_lower(p_lo, c0, j, t0)
        nprime = n - j
        Ls = []
        for t in range(1, t0 + 1):
            dmax = c0 + C * (j + 1) + v0f * t + 1
            ll = lr_lower_log(nprime, t, dmax, p_lo, p_hi)
            Ls.append(math.exp(ll) if ll > -math.inf else 0.0)
        for i in range(1, len(Ls)):
            Ls[i] = min(Ls[i], Ls[i - 1])
        hj = 0.0
        for t in range(1, t0 + 1):
            Lt = Ls[t - 1]
            Ln = Ls[t] if t < t0 else 0.0
            hj += (Lt - Ln) * cdf[t]
        hj *= (1 - 1e-9)
        hs.append(hj)
        num += q0 ** j * (1 - hj)
    num += q0 ** (J + 1) / (1 - q0)
    val = num / (1 - hi(h0))
    return dict(n=n, h0=hi(h0), hs=hs, value=val * (1 + 1e-9), ok=(val * (1 + 1e-9) < 1) and r['side'])


def tb_regime_Dw(mX, c0, v0, vmin, vsmall, J):
    """Room-weighted regime D' (valid for every m >= mX and k_b(m) < k <= k1(m)):
       q <= max(Q_L', Q_H) + R_B,  Q_L' = Top_L * Phi,  Phi = sum_{j<=J} r^j / sum_{j<=J} r^j D_j,
    D_j = max(vmin/c, 1 - Dip_j) (nondecreasing), Dip_j = e^{1/(6k)} e^{-eta_L j} + e^{-vmin^2 k - 4 vmin j}/(1-e^{-4vmin^2})
    for j >= 3 (room z > j, Lemma LR with delta >= 3 for the reversed walk), r = r_lo(mX) the fastest weight decay."""
    base = tb_regime_D(mX, c0, v0, vmin, vsmall)
    v0q, vminq = Fraction(v0), Fraction(vmin)
    V0, VMIN = ivq(v0q), ivq(vminq)
    kb = base['kb']
    kmin = kb + 1
    K1 = (mX + 1) / V0 + 2
    eta_i, eta_f = eta_lower(I_C - VMIN, Fraction(1, kmin))
    pref = iv.exp(iv.mpf(1) / (6 * kmin))
    Dj = []
    for j in range(0, J + 1):
        dj = VMIN / I_C
        if j >= 3:
            dip = pref * iv.exp(-eta_i * j) + iv.exp(-VMIN ** 2 * kmin - 4 * VMIN * j) / (1 - iv.exp(-4 * VMIN ** 2))
            cand = 1 - dip
            if cand.a > dj.a:
                dj = iv.mpf(cand.a)
        Dj.append(iv.mpf(dj.a))
    for j in range(1, len(Dj)):
        if Dj[j].a < Dj[j - 1].a:
            Dj[j] = Dj[j - 1]
    pstar_hi = I_C - (mX - 2) / K1
    rlo = (1 - pstar_hi - iv.mpf(J) / kmin) / (pstar_hi + iv.mpf(J + 1) / kmin)
    R = iv.mpf(rlo.a)
    numer = sum((R ** j for j in range(J + 1)), iv.mpf(0))
    denom = sum((R ** j * Dj[j] for j in range(J + 1)), iv.mpf(0))
    Phi = iv.mpf((numer / denom).b)
    QLw = base['QL'] * VMIN / I_C * Phi          # Top_L * Phi  (base QL = Top_L * c/vmin)
    qDw = iv.mpf(max(hi(QLw), hi(base['QH']))) + base['RB']
    # side: J <= floor(z_*) - 1 (class L contains j = 0..J), Lemma LR side condition for the dip (w <= 1/4)
    vmax = I_MU + (1 + iv.sqrt(kmin)) / kmin
    side = (base['side'] and (J + 1 <= base['zstar'].a) and (2 * (J + 2) / iv.mpf(kmin) + vmax).b <= 0.25
            and (kmin * (I_C / 2 - vmax)).a >= J + 4 and rlo.a > 0)
    return dict(base=base, Phi=Phi, rlo=rlo, QLw=QLw, qDw=qDw, side=side, ok=side and qDw.b < 1)
