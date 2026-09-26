#!/usr/bin/env python3
"""procgen_robin_20260926_lib -- machinery for the robin lane (session collatz-procgen-20260922, HYP-9142).

Setting (pairpeak note, sections 0-5).  c = log_3 2, sigma^2 = c(1-c), lambda = (1-c)/c, H = H_2(c) (bits).
Q = iid Bernoulli(c) letters; the letter walk S_t = e_t - c t has steps +(1-c) (prob c) and -c (prob 1-c).
Tilt identity (peak note, section 4): for every word w of length L, 1 = Q(w) 2^(H L) lambda^(S_L).

  N_m(L)  = survivors of the reflected barrier Pi^(m) (Robin walk S' in (0, m-c), zone [m-1, m-c) forced down,
            counted with weight 2 per zone step) = 2^(HL) E_Q[lambda^(S'_L + J_L); survive].
  A_M(L)  = #{u : S_t > 0 (1<=t<=L), S_t < M (0<=t<=L-1)} = 2^(HL) E_Q[lambda^(S_L); same event]   (M real).

Main objects of this lane.
* Free eigenfunctions (Lemma E).  For 0 < theta < pi/c put
      beta(theta)  = ln[ (1-c) sin(c theta) / (c sin((1-c) theta)) ]      (< 0),
      kappa(theta) = c e^(beta(1-c)) cos((1-c) theta) + (1-c) e^(-beta c) cos(c theta).
  Then f(s) = e^(beta s) sin(theta s + phi) satisfies  c f(s+1-c) + (1-c) f(s-c) = kappa(theta) f(s)  for all real s.
* Robin supersolution (Theorem 1): g(s) = e^(beta s) sin(theta (s + c)), theta = pi/(m+2).
* Dirichlet subsolution (Theorem 2): psi(s) = e^(beta s) sin(pi s / M).
* End lemma (Theorem 2): explicit two-phase bridge bound q(M) (cycle lemma + Hoeffding + Robbins).
* RM reduction (Theorem 6): hard-wall counts V_k, F_k.
"""
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import procgen_pairpeak_20260926_lib as PP   # read-only reuse: exact floors, N_m(L), A_M(L) DPs
import procgen_peak_20260926_lib as PK       # read-only reuse: H2, kappa_3, rho^peak float DP

C = math.log(2) / math.log(3)
SIG2 = C * (1 - C)
LAM = (1 - C) / C
HB = PK.H2(C)                               # H in bits (0.949955...)
ETA = 1 - HB                                # 1 - H = 0.050044
KAPPA3 = PK.kappa(3)[0]                     # (3/2)(pi^2 sigma^2)^(1/3) (ln 3)^(2/3) = 2.10758
V_UP = (2 * C - 0.5) * (1.5 - 2 * C)        # min of x(1-x) on [c, 2c - 1/2]  (= 0.18139)
THETA_STAR = (math.pi ** 2 * SIG2 / math.log(3)) ** (1 / 3)


def check(cond, msg):
    """every printed claim goes through here; a failed claim aborts the run"""
    if not cond:
        raise SystemExit("CHECK FAILED: " + msg)
    print("  ok:", msg, flush=True)


# ---------------------------------------------------------------------------------------------------
# 1. free eigenfunctions
# ---------------------------------------------------------------------------------------------------
def beta(th):
    return math.log((1 - C) * math.sin(th * C) / (C * math.sin(th * (1 - C))))


def kappa(th):
    b = beta(th)
    return C * math.exp(b * (1 - C)) * math.cos(th * (1 - C)) + (1 - C) * math.exp(-b * C) * math.cos(th * C)


def zone_F(th, d=2.0):
    """F(theta) = kappa sin(d theta) - 2(1-c) e^(-beta c) sin((d+c) theta); the Robin zone condition at the top
    zone point for the supersolution of width m + d (theta = pi/(m+d)) is F > 0"""
    b = beta(th)
    return kappa(th) * math.sin(d * th) - 2 * (1 - C) * math.exp(-b * C) * math.sin((d + C) * th)


def delta_R(m):
    """smallest width shift delta such that the zone condition holds at theta = pi/(m+delta) (bisection);
    the supersolution width is m + delta_R(m)"""
    def ok(d):
        th = math.pi / (m + d)
        return zone_F(th, d) > 0
    lo, hi = 0.5, 3.0
    for _ in range(80):
        mid = (lo + hi) / 2
        if ok(mid):
            hi = mid
        else:
            lo = mid
    return hi


# ---------------------------------------------------------------------------------------------------
# 2. Theorem 1 (Robin supersolution bound) and the Dirichlet supersolution (upper bound for confinement)
# ---------------------------------------------------------------------------------------------------
def robin_bound_log(L, m):
    """ln of the Theorem 1 bound  N_m(L) <= 2^(HL) e^(|beta|(m-c)) kappa(pi/(m+2))^L"""
    th = math.pi / (m + 2)
    return HB * L * math.log(2) + abs(beta(th)) * (m - C) + L * math.log(kappa(th))


def dirichlet_upper_log(n, W):
    """ln of the bound  Q_0(0 < S_j < W, 1 <= j <= n) <= (c/(1-c)) e^(|beta| W) kappa(pi/(W+1))^n"""
    th = math.pi / (W + 1)
    return math.log(C / (1 - C)) + abs(beta(th)) * W + n * math.log(kappa(th))


def theorem_A_log(L, m):
    """ln of the pairpeak note's Theorem A bound for N_m(L)/2^L (for comparison)"""
    return PP.block_bound_log(L, m)


# ---------------------------------------------------------------------------------------------------
# 3. Theorem 2 (Dirichlet lower bound): end lemma constants
# ---------------------------------------------------------------------------------------------------
def n_of_M(M):
    """largest integer n >= 1 with 4 n ln(2n) <= (M-2)^2  (Hoeffding condition of the end lemma); 0 if none"""
    t = (M - 2) ** 2
    if t < 4 * math.log(2):
        return 0
    lo, hi = 1, 1
    while 4 * hi * math.log(2 * hi) <= t:
        hi *= 2
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if 4 * mid * math.log(2 * mid) <= t:
            lo = mid
        else:
            hi = mid
    return lo


def nmin_of_M(M):
    """n_min(M) = ceil((M/2 + 1)/(c - 1/2)): keeps the binomial parameter u/n inside [1/2, 2c - 1/2]"""
    return math.ceil((M / 2 + 1) / (C - 0.5))


def q_end_log(M, n):
    """ln q(M, n) = -ln(16 n^3) - (M/2+1)^2 (1/v_up + 1/sigma^2)/(2n)   (end lemma)"""
    return -math.log(16 * n ** 3) - (M / 2 + 1) ** 2 * (1 / V_UP + 1 / SIG2) / (2 * n)


def M1_threshold(Mmax=3000):
    """smallest integer M1 such that n(M) >= n_min(M) for every real M in [M1, Mmax]: checked on the grid
    M_a = M1 + j/4 in the interval-safe form n(M_a) >= n_min(M_a + 1/4) (both functions are non-decreasing in M).
    Beyond Mmax = 3000: n(M) >= (M-2)^2/(8 ln(M-2)) >= 4M >= n_min(M), so n(M) >= n_min(M) analytically."""
    last_bad = None
    j = 16
    while j / 4 <= Mmax:
        Ma = j / 4
        if n_of_M(Ma) < nmin_of_M(Ma + 0.25):
            last_bad = Ma + 0.25
        j += 1
    return math.floor(last_bad) + 1 if last_bad is not None else 4


def P2_log(M, M1):
    """ln P_2(M): A_M(L) >= 2^(HL) kappa(pi/M)^L / P_2(M) for all L >= 1 (Theorem 2)"""
    th = math.pi / M
    b = abs(beta(th))
    crude = math.log(M) + b - math.log(2 * C * (1 - C)) - M * math.log(LAM)
    if M < M1:
        return crude
    n = n_of_M(M)
    long_ = math.log(M) + b - math.log(2 * C * (1 - C) * LAM) - q_end_log(M, n)
    short = math.log(4) + 1.5 * math.log(2 * n) + 1 / (16 * V_UP) - math.log(C * LAM)
    tiny = 8 * HB * math.log(2)
    return min(crude, max(long_, short, tiny))


def P_log(m, M1):
    """ln P(m): N_m(L) <= P(m) A_(m+2)(L)  (Corollary 3)"""
    th = math.pi / (m + 2)
    return abs(beta(th)) * (m - C) + P2_log(m + 2, M1)


# ---------------------------------------------------------------------------------------------------
# 4. exact / float counts (reuse of the pairpeak DPs) and Q-probability DPs
# ---------------------------------------------------------------------------------------------------
def floor_table(Mx):
    return PP.floor_table(Mx)


def robin_float_log(L, m, fl, Mx):
    """log N_m(t) for t = 0..L (counting measure), renormalised float DP (zone states counted twice)"""
    f = np.zeros(L + 3)
    f[0] = 1.0
    logs = 0.0
    out = [0.0]
    for t in range(L):
        ft = 0 if t == 0 else int(fl[Mx + t])
        thr = max((m - 1) if t == 0 else (m - 1 + ft + 1), 0)
        nb = f.copy()
        nb[thr:] = 0.0
        bb = f.copy()
        bb[:thr] = 0.0
        g = nb.copy()
        g[1:] += nb[:-1]
        g += 2 * bb
        g[:int(fl[Mx + t + 1]) + 1] = 0.0
        s = g.sum()
        if s == 0:
            out.append(-math.inf)
            f = g
            continue
        logs += math.log(s)
        f = g / s
        out.append(logs)
    return out


def hard_float_log(L, Mw, fl, Mx):
    """log A_Mw(t) for t = 0..L, renormalised float DP (Mw a positive integer)"""
    f = np.zeros(L + 3)
    f[0] = 1.0
    logs = 0.0
    out = [0.0]
    for t in range(L):
        ub = (Mw - 1) if t == 0 else (Mw + int(fl[Mx + t]))
        if ub + 1 < len(f):
            f[max(ub + 1, 0):] = 0.0
        g = f.copy()
        g[1:] += f[:-1]
        g[:int(fl[Mx + t + 1]) + 1] = 0.0
        s = g.sum()
        logs += math.log(s)
        f = g / s
        out.append(logs)
    return out


def bridge_confined_fraction(x, xp, n, lo, hi):
    """exact fraction of the arrangements of the bridge x -> xp (n steps) whose path stays in (lo, hi) at steps
    1..n; returns (fraction, number of ones u) or (None, None) if xp is not reachable"""
    u = xp - x + C * n
    ui = round(u)
    if abs(u - ui) > 1e-9 or ui < 0 or ui > n:
        return None, None
    # count lattice paths with U ups among n steps, positions x + U_k - c k in (lo, hi) for k = 1..n
    f = {0: 1}
    for k in range(1, n + 1):
        g = {}
        for U, v in f.items():
            for dU in (0, 1):
                V = U + dU
                if V > ui or (k - V) > (n - ui):
                    continue
                p = x + V - C * k
                if lo < p < hi:
                    g[V] = g.get(V, 0) + v
        f = g
    good = f.get(ui, 0)
    return good / math.comb(n, ui), ui


def q_bridge_min(M, n, ys):
    """numerical min over the given starting points y of Q_y(S_t in (0,M) for 1<=t<=2n, S_2n in (0,1]) (float DP)"""
    best = 1.0
    for y in ys:
        f = np.zeros(2 * n + 2)
        f[0] = 1.0
        for k in range(1, 2 * n + 1):
            g = (1 - C) * f
            g[1:] += C * f[:-1]
            U = np.arange(len(g))
            pos = y + U - C * k
            g[(pos <= 0) | (pos >= M)] = 0.0
            f = g
        U = np.arange(len(f))
        pos = y + U - C * 2 * n
        val = f[(pos > 0) & (pos <= 1)].sum()
        best = min(best, val)
    return best


# ---------------------------------------------------------------------------------------------------
# 5. RM reduction: hard-wall counts V_k(x) and zone-avoiding counts F_k(x) from a real start x
# ---------------------------------------------------------------------------------------------------
def counts_from(x, hi, kmax):
    """list over k = 0..kmax of #{w in {0,1}^k : x + S_j in (0, hi) for 1 <= j <= k-1, x + S_k > 0}
    (the start x itself is not constrained; exact integers)"""
    f = {0: 1}
    out = []
    for j in range(kmax + 1):
        out.append(sum(v for U, v in f.items() if (j == 0 or x + U - C * j > 0)))
        g = {}
        for U, v in f.items():
            p = x + U - C * j
            if j == 0 or (0 < p < hi):
                g[U] = g.get(U, 0) + v
                g[U + 1] = g.get(U + 1, 0) + v
        f = {U: v for U, v in g.items() if x + U - C * (j + 1) > 0}
    return out
