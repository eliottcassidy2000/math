#!/usr/bin/env python3
"""procgen_pairpeak_20260926_lib -- machinery for the pairing-peak lane (session collatz-procgen-20260922).

Setting.  T(n) = n/2 (n even), (3n+1)/2 (n odd).  Pairing family: pairs {2i-1, 2i}; flipping pair i sends
2i-1 -> i-1 and 2i -> 3i.  c = log_3 2, sigma^2 = c(1-c), lambda = (1-c)/c, 1-H = 1 - H_2(c).

Main objects.
* The coupling chain of a single flip (Lemma 1 of the note): for odd v, x_1 = (v-1)/2 (flipped image) and
  y_1 = (3v+1)/2 (Collatz image) satisfy y_t = 3^(a_t) x_t + b_t along the two Collatz orbits.
* The reflected barrier process Pi^(m) (section 2 of the note): flip every odd point whose slope level
  S'_t = U_t - c t is >= m-1 (U_t = number of unflipped odd steps so far).  The modified parity word is a
  bijective image of n mod 2^L (Lemma 2), so the proportion of classes on which the process survives
  (S'_t > 0 for 1 <= t <= L) is N_m(L)/2^L, computed exactly by a DP over U.
* Hard-wall counts A_M(L) = #{u : S_t > 0 (1<=t<=L), S_t < M (0<=t<=L-1)} (bad words with peak < 3^M).

All floors floor(t c) are exact (integer comparison 3^n <= 2^t), via the peak lane's floor_table.
"""
import math
import random
import sys
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import procgen_peak_20260926_lib as PK  # read-only reuse: exact floors, rho_L, rho^peak_L

C = math.log(2) / math.log(3)
SIG2 = C * (1 - C)
LAM = (1 - C) / C
ETA = 1 - PK.H2(C)                      # 1 - H(log_3 2) = 0.050044...
RHO_B = 1 - (1 - LAM ** 2) / 14.0       # per-block factor of Theorem A


def check(cond, msg):
    """every printed claim goes through here; a failed claim aborts the run"""
    if not cond:
        raise SystemExit("CHECK FAILED: " + msg)
    print("  ok:", msg, flush=True)


def T(x):
    return x >> 1 if not (x & 1) else (3 * x + 1) >> 1


def floor_table(M):
    return PK.floor_table(3, M)


# ---------------------------------------------------------------------------------------------------
# 1. coupling chain of a single flip
# ---------------------------------------------------------------------------------------------------
def chain_step(a, b, py, px):
    """one step of (a, b) given the parities of y_t (py) and x_t (px); requires a >= 1 when py=0, px=1"""
    if py == 0 and px == 0:
        assert b % 2 == 0
        return a, b // 2
    if py == 1 and px == 1:
        assert b % 2 == 0
        n = 3 * b + 1 - 3 ** a
        assert n % 2 == 0
        return a, n // 2
    if py == 1 and px == 0:
        assert b % 2 == 1
        return a + 1, (3 * b + 1) // 2
    assert a >= 1 and b % 2 == 1
    n = b - 3 ** (a - 1)
    assert n % 2 == 0
    return a - 1, n // 2


def chain_verify(v, steps):
    """follow x_t, y_t and (a_t, b_t) from x_1 = (v-1)/2, y_1 = (3v+1)/2; verify y_t = 3^a x_t + b_t while a >= 1;
    return (merge_time or None, first time a == 0 or None, number of verified steps)"""
    x, y = (v - 1) // 2, (3 * v + 1) // 2
    a, b = 1, 2
    t = 1
    while t <= steps:
        assert y == 3 ** a * x + b, (v, t, a, b)
        if x == y:
            return t, t, t
        if a == 0:
            return None, t, t
        px, py = x & 1, y & 1
        assert px == (py + b) % 2
        a, b = chain_step(a, b, py, px)
        x, y = T(x), T(y)
        t += 1
    return None, None, t - 1


def first_opportunity_merge(v):
    """True iff the Collatz word of y_1 = (3v+1)/2 begins 1^r 0 0 (r >= 0)"""
    y = (3 * v + 1) // 2
    while y & 1:
        y = T(y)
    y = T(y)
    return (y & 1) == 0


# ---------------------------------------------------------------------------------------------------
# 2. reflected barrier process: exact survival counts N_m(L), hard-wall counts A_M(L)
# ---------------------------------------------------------------------------------------------------
def refl_count(L, m, fl, M):
    """exact N_m(L): number of words u'' in {0,1}^L whose barrier-(m-1) reflected slope walk survives.
    State: U = number of unflipped odd steps; level S' = U - c t.  At a barrier state (S' >= m-1) both letters
    give a down step (odd -> flipped, even -> halved), so the count doubles there."""
    f = {0: 1}
    for t in range(L):
        ft = 0 if t == 0 else int(fl[M + t])
        thr = (m - 1) if t == 0 else (m - 1 + ft + 1)          # S'_t >= m-1  <=>  U >= thr
        nf = {}
        for U, v in f.items():
            if U >= thr:
                nf[U] = nf.get(U, 0) + 2 * v
            else:
                nf[U + 1] = nf.get(U + 1, 0) + v
                nf[U] = nf.get(U, 0) + v
        lo = int(fl[M + t + 1]) + 1                             # S'_(t+1) > 0  <=>  U >= floor(c(t+1)) + 1
        f = {U: v for U, v in nf.items() if U >= lo}
    return sum(f.values())


def hard_count(L, Mx, fl, M):
    """exact A_Mx(L) = #{u : S_t > 0 for 1<=t<=L and S_t < Mx for 0<=t<=L-1}"""
    f = {0: 1}
    for t in range(L):
        ub = (Mx - 1) if t == 0 else (Mx + int(fl[M + t]))       # S_t < Mx  <=>  U <= ub
        f = {U: v for U, v in f.items() if U <= ub}
        nf = {}
        for U, v in f.items():
            nf[U + 1] = nf.get(U + 1, 0) + v
            nf[U] = nf.get(U, 0) + v
        lo = int(fl[M + t + 1]) + 1
        f = {U: v for U, v in nf.items() if U >= lo}
    return sum(f.values())


def refl_float(L, m, fl, M):
    """float N_m(L)/2^L and E[J] (expected number of flips of Pi^(m) over all words, stopped at death)"""
    f = np.zeros(L + 2)
    f[0] = 1.0
    EJ = 0.0
    for t in range(L):
        ft = 0 if t == 0 else int(fl[M + t])
        thr = max((m - 1) if t == 0 else (m - 1 + ft + 1), 0)
        nb = f.copy()
        nb[thr:] = 0.0
        bb = f.copy()
        bb[:thr] = 0.0
        EJ += 0.5 * bb.sum()
        g = 0.5 * nb
        g[1:] += 0.5 * nb[:-1]
        g += bb
        g[:int(fl[M + t + 1]) + 1] = 0.0
        f = g
    return float(f.sum()), EJ


def hard_float(L, Mx, fl, M):
    f = np.zeros(L + 2)
    f[0] = 1.0
    for t in range(L):
        ub = (Mx - 1) if t == 0 else (Mx + int(fl[M + t]))
        if ub + 1 < len(f):
            f[max(ub + 1, 0):] = 0.0
        g = 0.5 * f
        g[1:] += 0.5 * f[:-1]
        g[:int(fl[M + t + 1]) + 1] = 0.0
        f = g
    return float(f.sum())


def block_bound(L, m):
    """Theorem A: N_m(L)/2^L <= 2^(-(1-H)L) rho_b^floor(L/n_m), n_m = ceil(2(m+1)^2/sigma^2)"""
    n = math.ceil(2 * (m + 1) ** 2 / SIG2)
    return 2.0 ** (-ETA * L) * RHO_B ** (L // n)


def block_bound_log(L, m):
    n = math.ceil(2 * (m + 1) ** 2 / SIG2)
    return -ETA * L * math.log(2) + (L // n) * math.log(RHO_B)


# ---------------------------------------------------------------------------------------------------
# 3. the barrier process on actual integers
# ---------------------------------------------------------------------------------------------------
def barrier_run(n, L, m, fl, M):
    """run Pi^(m) on the integer n for L steps (slope-level rule); return (word, survived_by_slope, flips,
    descended_actually, cost, nflips_before_descent) where cost = sum of n/i over the flips made BEFORE the first
    descent (i = pair index of the flipped point; this is the cost of n's private certificate)"""
    x = n
    U = 0
    word = []
    flips = []
    surv = True
    desc = False
    cost = 0.0
    nfd = 0
    for t in range(L):
        ft = 0 if t == 0 else int(fl[M + t])
        thr = (m - 1) if t == 0 else (m - 1 + ft + 1)
        b = x & 1
        word.append(b)
        if b and U >= thr:
            flips.append((t, x))
            if not desc:
                cost += n / ((x + 1) // 2)
                nfd += 1
            x = (x - 1) >> 1
        elif b:
            x = (3 * x + 1) >> 1
            U += 1
        else:
            x >>= 1
        if not (U >= int(fl[M + t + 1]) + 1):
            surv = False
        if x < n:
            desc = True
    return word, surv, flips, desc, cost, nfd


def n_from_word(w, rng, extra_bits=64):
    """an integer n (with random high bits) whose Collatz parity word begins with w (Terras inverse)"""
    y, mod = 0, 1
    for j in range(len(w)):
        for cand in (y, y + mod):
            v = cand
            for _ in range(j):
                v = T(v)
            if (v & 1) == w[j]:
                y = cand
                break
        mod <<= 1
    return y + mod * ((rng.getrandbits(extra_bits) | 1) << 20)


def bad_word_sampler(L, fl, M):
    """exact backward counts g[t][e] of bad continuations; returns (sampler, |Bad_L|)"""
    g = [None] * (L + 1)
    g[L] = {e: 1 for e in range(L + 1)}
    for t in range(L - 1, -1, -1):
        d = {}
        lo = int(fl[M + t + 1]) + 1
        for e in range(t + 1):
            s = 0
            for bb in (0, 1):
                if e + bb >= lo:
                    s += g[t + 1].get(e + bb, 0)
            if s:
                d[e] = s
        g[t] = d
    total = g[0][0]

    def sample(rng):
        e = 0
        w = []
        for t in range(L):
            lo = int(fl[M + t + 1]) + 1
            c0 = g[t + 1].get(e, 0) if e >= lo else 0
            c1 = g[t + 1].get(e + 1, 0) if e + 1 >= lo else 0
            r = rng.randrange(c0 + c1)
            bb = 0 if r < c0 else 1
            w.append(bb)
            e += bb
        return w
    return sample, total


def one_fd_violation(v, L):
    """first j <= L with T^(j-1)((v-1)/2) > T^j(v) (one-flip domination fails), else None"""
    x, y = (v - 1) >> 1, (3 * v + 1) >> 1
    for j in range(1, L + 1):
        if x > y:
            return j
        x, y = T(x), T(y)
    return None


def refl_counts_all_L(Lmax, m, fl, M):
    """exact N_m(L) for L = 0..Lmax in one pass (list indexed by L)"""
    f = {0: 1}
    out = [1]
    for t in range(Lmax):
        ft = 0 if t == 0 else int(fl[M + t])
        thr = (m - 1) if t == 0 else (m - 1 + ft + 1)
        nf = {}
        for U, v in f.items():
            if U >= thr:
                nf[U] = nf.get(U, 0) + 2 * v
            else:
                nf[U + 1] = nf.get(U + 1, 0) + v
                nf[U] = nf.get(U, 0) + v
        lo = int(fl[M + t + 1]) + 1
        f = {U: v for U, v in nf.items() if U >= lo}
        out.append(sum(f.values()))
    return out


def hard_counts_all_L(Lmax, Mx, fl, M):
    """exact A_Mx(L) for L = 0..Lmax in one pass: the wall S_t < Mx is imposed at t = 0..L-1 only"""
    f = {0: 1}
    out = [1]
    for t in range(Lmax):
        ub = (Mx - 1) if t == 0 else (Mx + int(fl[M + t]))
        f = {U: v for U, v in f.items() if U <= ub}
        nf = {}
        for U, v in f.items():
            nf[U + 1] = nf.get(U + 1, 0) + v
            nf[U] = nf.get(U, 0) + v
        lo = int(fl[M + t + 1]) + 1
        f = {U: v for U, v in nf.items() if U >= lo}
        out.append(sum(f.values()))
    return out
