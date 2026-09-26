#!/usr/bin/env python3
"""procgen_family27_20260926_theory.py

Exact and high-precision computations for the family-27 note (no integer orbits here).

  * the Moran function of the owner's backward recursion
        g(s) = 2^(-s) + (1/3)(3/2)^s      (Lagarias-Weiss M_BP(-s))
    and the forward pressure phi(t) = ((3/2)^t + (1/2)^t)/2, with g(s) = phi(s-1);
    its roots, derivatives, minimum, and the tangent slope beta_BP = max_s (-ln g(s))/s;
  * the Lagarias-Weiss delay spectrum 1 - a*gLW(1/a) and its entropy form;
  * W_k = |Bad_k| (words with 3^(o_j) > 2^j for all j <= k), exact big integers;
  * the exact word-level rise law: N_k(W) = #{u in {0,1}^k : max_j 3^(o_j)/2^j >= W} and
    eps_k(W) = 2^-k sum_{tau > k} M_k, both exact (integers / Fractions), for W = 2^(i/2);
  * P(W) = P(sup_j M_j >= W) for the fair walk, by a pruned floating DP with a Ville tail bound.

Every function is deterministic; the runner imports this module and wraps each claim in check().
"""
from fractions import Fraction
import math

import mpmath as mp

mp.mp.dps = 40

LN2 = mp.log(2)
LN3 = mp.log(3)
ALPHA = LN3 / LN2            # log_2 3
RHO = LN2 / LN3              # log_3 2


def g(s):
    s = mp.mpf(s)
    return mp.power(2, -s) + mp.power(mp.mpf(3) / 2, s) / 3


def phi(t):
    t = mp.mpf(t)
    return (mp.power(mp.mpf(3) / 2, t) + mp.power(mp.mpf(1) / 2, t)) / 2


def hbin(p):
    p = mp.mpf(p)
    if p <= 0 or p >= 1:
        return mp.mpf(0)
    return -(p * mp.log(p) + (1 - p) * mp.log(1 - p)) / LN2


def Hnat(p):
    p = mp.mpf(p)
    return -(p * mp.log(p) + (1 - p) * mp.log(1 - p))


def g_exact(s_int):
    """g at an integer s as an exact Fraction."""
    return Fraction(1, 2 ** s_int) + Fraction(3 ** s_int, 3 * 2 ** s_int) if s_int >= 0 else None


def moran_constants():
    """All constants read off g.  Returns a dict of mpf values."""
    c = {}
    c["g1"] = g(1)
    c["g2"] = g(2)
    dg = lambda s: mp.diff(g, s)
    c["dg1"] = dg(1)
    c["dg2"] = dg(2)
    c["R1"] = -1 / c["dg1"]                           # residue at s = 1 (typical delay)
    c["R1_closed"] = 1 / mp.log(2 / mp.sqrt(3))
    c["R2"] = 1 / c["dg2"]                            # residue at s = 2 (time to peak)
    c["R2_closed"] = 1 / (mp.mpf(3) / 4 * LN3 - LN2)
    c["sstar"] = mp.findroot(dg, 1.5)                 # argmin g
    c["sstar_closed"] = 1 + mp.log(LN2 / mp.log(mp.mpf(3) / 2)) / LN3
    c["gmin"] = g(c["sstar"])
    c["h"] = hbin(RHO)
    c["gmin_closed"] = mp.power(2, -(1 - c["h"]))
    # beta_BP = max_s (-ln g(s))/s ; stationary point: s g'(s)/g(s) = ln g(s)
    f = lambda s: -mp.log(g(s)) / s
    sopt = mp.findroot(lambda s: mp.diff(f, s), 1.4)
    c["s_opt"] = sopt
    c["beta_BP"] = f(sopt)
    c["gamma_BP"] = 1 / c["beta_BP"]
    # entropy form: H(p) = p ln 3 on (1/2, log_3 2)
    pstar = mp.findroot(lambda p: Hnat(p) - p * LN3, 0.609)
    c["pstar"] = pstar
    c["gamma_entropy"] = 1 / (LN2 - pstar * LN3)
    # tilted odd fraction at s_opt: (1/3)(3/2)^s / g(s)
    c["p_at_sopt"] = (mp.power(mp.mpf(3) / 2, sopt) / 3) / g(sopt)
    # Lagarias-Weiss RRW form: gamma solves gamma * gLW(1/gamma) = 1
    c["gamma_LW"] = mp.findroot(lambda a: a * gLW(1 / a) - 1, 41.6)
    c["glide_const"] = 1 / (1 - c["h"])               # glide records, T-steps per log2 n
    c["glide_const_std"] = (1 + RHO) / (1 - c["h"])   # standard-map steps per log2 n
    c["zeta"] = mp.mpf(3) / 4 * ALPHA - 1             # window rise junction beta - 1
    c["one_minus_h34"] = 1 - hbin(mp.mpf(3) / 4)
    c["tilted_drift_bits"] = mp.diff(phi, 1) / LN2
    c["drift"] = mp.log(2 / mp.sqrt(3))
    return c


def lnM_RRW(t):
    t = mp.mpf(t)
    return mp.log(mp.power(2, t) / 2 + mp.power(mp.mpf(2) / 3, t) / 2)


def gLW(a):
    """Lagarias-Weiss rate function sup_theta (theta a - ln M_RRW(theta))."""
    a = mp.mpf(a)
    th = mp.findroot(lambda t: mp.diff(lnM_RRW, t) - a, 0 if a > 0.1 else -0.4)
    return th * a - lnM_RRW(th)


def gLW_theta(a):
    a = mp.mpf(a)
    return mp.findroot(lambda t: mp.diff(lnM_RRW, t) - a, 0 if a > 0.1 else -0.4)


def delay_spectrum_LW(alpha):
    alpha = mp.mpf(alpha)
    return 1 - alpha * gLW(1 / alpha)


def delay_spectrum_entropy(alpha):
    alpha = mp.mpf(alpha)
    p = (LN2 - 1 / alpha) / LN3
    return alpha * (Hnat(p) - p * LN3)


def window_rise_spectrum(beta):
    """E_win(beta): exponent of #{n in [2^k,2^(k+1)) : max_{j<=k} T^j n >= n^beta} (PROVED)."""
    beta = mp.mpf(beta)
    zeta = mp.mpf(3) / 4 * ALPHA - 1
    if beta <= 1:
        return mp.mpf(1)
    if beta <= 1 + zeta:
        return 2 - beta
    if beta <= ALPHA:
        return hbin(beta / ALPHA)
    return None


def window_rise_spectrum_bruteforce(beta, grid=4000):
    """max_t [1 - t + t h((t+b)/(t alpha))] over t in (0,1], b = beta-1 (word counting exponent)."""
    b = mp.mpf(beta) - 1
    best = mp.mpf(-1)
    for i in range(1, grid + 1):
        t = mp.mpf(i) / grid
        q = (t + b) / (t * ALPHA)
        if q > 1:
            continue
        val = 1 - t + t * (hbin(q) if q >= mp.mpf(1) / 2 else 1)
        if val > best:
            best = val
    return best


# ---------------------------------------------------------------- exact word counts
def W_counts(K):
    """W_k = #{u in {0,1}^k : 3^(o_j) > 2^j for 1 <= j <= k}, k = 0..K (W_0 = 1)."""
    W = [1]
    f = {0: 1}                 # o -> number of admissible prefixes of current length
    for j in range(1, K + 1):
        nf = {}
        pw = 2 ** j
        for o, c in f.items():
            for o2 in (o, o + 1):
                if 3 ** o2 > pw:
                    nf[o2] = nf.get(o2, 0) + c
        f = nf
        W.append(sum(f.values()))
    return W


def _ge_threshold(o, j, i, k=None, shift=False):
    """exact test 3^o / 2^j >= 2^(i/2)  (or >= 2^(i/2) - (3/4)^k if shift)."""
    if not shift:
        if i % 2 == 0:
            return 3 ** o >= 2 ** (j + i // 2)
        return 3 ** (2 * o) >= 2 ** (2 * j + i)
    # 3^o/2^j >= 2^(i/2) - (3/4)^k  <=>  3^o 4^k + 2^j 3^k >= 2^(j + i/2) 4^k
    lhs = 3 ** o * 4 ** k + 2 ** j * 3 ** k
    if i % 2 == 0:
        return lhs >= 2 ** (j + i // 2) * 4 ** k
    return lhs * lhs >= 2 ** (2 * j + i) * 16 ** k


def rise_words(k, i, shift=False):
    """Exact (N_k(W), eps_k(W)) for W = 2^(i/2) (or W - (3/4)^k if shift):
    N_k = #{u : max_{0<=j<=k} M_j(u) >= W},  eps_k = 2^-k sum_{u: tau>k} M_k(u) (Fraction)."""
    f = {0: 1}    # o -> number of prefixes (length j) not yet absorbed
    absorbed = 0
    for j in range(1, k + 1):
        nf = {}
        for o, c in f.items():
            for o2 in (o, o + 1):
                if _ge_threshold(o2, j, i, k, shift):
                    absorbed += c * 2 ** (k - j)
                else:
                    nf[o2] = nf.get(o2, 0) + c
        f = nf
    eps = Fraction(sum(c * 3 ** o for o, c in f.items()), 4 ** k)
    return absorbed, eps


def cond_mean_overshoot(k, i):
    """E[M_tau | tau <= k] exactly (Fraction), W = 2^(i/2)."""
    f = {0: 1}
    s = Fraction(0)
    cnt = 0
    for j in range(1, k + 1):
        nf = {}
        for o, c in f.items():
            for o2 in (o, o + 1):
                if _ge_threshold(o2, j, i):
                    w = c * 2 ** (k - j)
                    cnt += w
                    s += Fraction(w * 3 ** o2, 2 ** j)
                else:
                    nf[o2] = nf.get(o2, 0) + c
        f = nf
    return (s / cnt) if cnt else None, cnt


# ---------------------------------------------------------------- P(W) numerically
def riser_probability(W, jmax=6000, prune=1e-40):
    """P(sup_j M_j >= W) for the fair walk M_j = 3^(o_j)/2^j, by a pruned DP.
    Returns (lower, upper): lower = absorbed mass by jmax, upper = lower + (Ville) remaining
    mass bound sum_{alive} P(state) M/W + pruned bound."""
    lnW = math.log(W)
    l3, l2 = math.log(3.0), math.log(2.0)
    f = {0: 1.0}
    absorbed = 0.0
    pruned_bound = 0.0
    for j in range(1, jmax + 1):
        nf = {}
        for o, p in f.items():
            half = 0.5 * p
            for o2 in (o, o + 1):
                x = o2 * l3 - j * l2
                if x >= lnW:
                    absorbed += half
                else:
                    nf[o2] = nf.get(o2, 0.0) + half
        # prune states whose Ville bound P*M/W is negligible
        f = {}
        for o, p in nf.items():
            M = math.exp(o * l3 - j * l2)
            if p * M / W < prune:
                pruned_bound += p * M / W
            else:
                f[o] = p
        if not f:
            break
    remaining = sum(p * math.exp(o * l3 - jmax * l2) / W for o, p in f.items())
    return absorbed, absorbed + remaining + pruned_bound


if __name__ == "__main__":
    c = moran_constants()
    for k_, v in c.items():
        print(k_, mp.nstr(v, 15))
