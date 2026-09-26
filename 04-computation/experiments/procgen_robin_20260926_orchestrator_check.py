#!/usr/bin/env python3
"""Orchestrator audit of lane `robin`, written from the note's statements; the
lane's scripts were not read.

  1. Interval arithmetic (mpmath.iv, my own bisection), the three
     computer-assisted steps:
     (a) ln kappa(theta) < -0.1 theta^2 on [0.15, pi/2];
     (b) Lemma Z: F(theta) = kappa sin 2theta - 2(1-c) e^(-beta c) sin((2+c)theta) > 0 on [0.15, pi/4];
     (c) end lemma grid: n(M) >= n_min(M + 1/4) on the grid M = 110, 110.25, ..., 3000.
  2. Lemma E(i): c f(s+1-c) + (1-c) f(s-c) = kappa f(s) for f = e^(beta s) sin(theta s + phi).
  3. Exact counts: N_m(L) (reflected barrier) against Theorem 1, for m = 2..12, L <= 200;
     A_M(L) (hard wall) against Theorem 2's crude constant, for M = 4..16, L <= 200;
     and N_m(L) <= A_(m+2)(L) * (m+2)^17.
"""
import math, random
import mpmath
from mpmath import iv, mpf

mpmath.mp.dps = 40
iv.dps = 30


def check(cond, msg):
    if not cond:
        raise SystemExit("CHECK FAILED: " + msg)
    print("  ok:", msg)


# ---------------------------------------------------------------- interval functions
Civ = iv.log(2) / iv.log(3)


def beta_iv(th):
    return iv.log((1 - Civ) * iv.sin(Civ * th) / (Civ * iv.sin((1 - Civ) * th)))


def kappa_iv(th):
    b = beta_iv(th)
    return Civ * iv.exp(b * (1 - Civ)) * iv.cos((1 - Civ) * th) + (1 - Civ) * iv.exp(-b * Civ) * iv.cos(Civ * th)


def bisect_prove(fun_upper_negative, lo, hi, max_boxes=200000):
    """prove fun(box) < 0 on [lo, hi] (interval endpoints given as intervals: use outer bounds)."""
    stack = [(lo, hi)]
    boxes = 0
    while stack:
        a, b = stack.pop()
        box = iv.mpf([a, b])
        val = fun_upper_negative(box)
        if val.b < 0:
            boxes += 1
            continue
        if b - a < mpf("1e-12"):
            return False, boxes
        mid = (a + b) / 2
        stack.append((a, mid)); stack.append((mid, b))
        if boxes > max_boxes:
            return False, boxes
    return True, boxes


pi_iv = iv.pi
lo = mpf("0.15")
hi_half = (pi_iv / 2).b      # outward enclosure of pi/2
hi_quarter = (pi_iv / 4).b   # outward enclosure of pi/4
ok, nb = bisect_prove(lambda th: iv.log(kappa_iv(th)) + iv.mpf("0.1") * th * th, lo, hi_half)
check(ok, f"(a) ln kappa(theta) < -0.1 theta^2 on [0.15, pi/2] (independent interval bisection, {nb} boxes)")


def negF(th):
    b = beta_iv(th)
    k = kappa_iv(th)
    F = k * iv.sin(2 * th) - 2 * (1 - Civ) * iv.exp(-b * Civ) * iv.sin((2 + Civ) * th)
    return -F


ok, nb = bisect_prove(negF, lo, hi_quarter)
check(ok, f"(b) Lemma Z: F(theta) > 0 on [0.15, pi/4] (independent interval bisection, {nb} boxes)")

# (c) end lemma grid
c = mpmath.log(2) / mpmath.log(3)


def n_of(Mv):
    # largest n with 4 n ln(2n) <= (M-2)^2 (monotone in n)
    target = (Mv - 2) ** 2
    lo_n, hi_n = 1, 1
    while 4 * hi_n * mpmath.log(2 * hi_n) <= target:
        hi_n *= 2
    while hi_n - lo_n > 1:
        mid = (lo_n + hi_n) // 2
        if 4 * mid * mpmath.log(2 * mid) <= target:
            lo_n = mid
        else:
            hi_n = mid
    return lo_n


def nmin_of(Mv):
    return int(mpmath.ceil((Mv / 2 + 1) / (c - mpf(1) / 2)))


bad = 0
Mv = mpf(110)
while Mv <= 3000:
    if n_of(Mv) < nmin_of(Mv + mpf(1) / 4):
        bad += 1
    Mv += mpf(1) / 4
check(bad == 0, "(c) end lemma: n(M_a) >= n_min(M_a + 1/4) on the grid 110..3000 (step 1/4)")

# 2. Lemma E(i)
random.seed(1)
cf = mpmath.log(2) / mpmath.log(3)
worst = 0
for _ in range(200):
    th = mpf(random.uniform(0.01, 3.0 / float(cf) * 0.99))
    phi = mpf(random.uniform(0, 6.28))
    s = mpf(random.uniform(-5, 5))
    be = mpmath.log((1 - cf) * mpmath.sin(cf * th) / (cf * mpmath.sin((1 - cf) * th)))
    ka = cf * mpmath.exp(be * (1 - cf)) * mpmath.cos((1 - cf) * th) + (1 - cf) * mpmath.exp(-be * cf) * mpmath.cos(cf * th)
    f = lambda x: mpmath.exp(be * x) * mpmath.sin(th * x + phi)
    lhs = cf * f(s + 1 - cf) + (1 - cf) * f(s - cf)
    worst = max(worst, abs(lhs - ka * f(s)))
check(worst < mpf("1e-30"), "Lemma E(i): the eigen-identity holds to 1e-30 at 200 random (theta, phi, s)")

# 3. exact counts
def N_barrier(m, L):
    """reflected barrier: state U (up count), time t; S' = U - c t; zone [m-1, m-c)."""
    cur = {0: 1}
    out = [1]
    for t in range(0, L):
        nxt = {}
        for U, cnt in cur.items():
            in_zone = (3 ** (U - m + 1) >= 2 ** t if U - m + 1 >= 0 else False) and \
                      ((U - m < 0) or 3 ** (U - m) < 2 ** (t - 1) if t >= 1 else U - m < 0)
            if in_zone:
                cands = [(U, 2 * cnt)]              # both letters: S' -> S' - c
            else:
                cands = [(U + 1, cnt), (U, cnt)]
            for U2, w in cands:
                if 3 ** U2 > 2 ** (t + 1):          # S'_{t+1} > 0
                    nxt[U2] = nxt.get(U2, 0) + w
        cur = nxt
        out.append(sum(cur.values()))
    return out


def A_hard(M, L):
    """words with S_t > 0 (1<=t<=L) and S_t < M (0<=t<=L-1), S_t = e_t - c t."""
    cur = {0: 1}
    out = [1]
    for t in range(0, L):
        nxt = {}
        for e, cnt in cur.items():
            # S_t < M required at time t (t <= L-1)
            if not (e - M < 0 or 3 ** (e - M) < 2 ** t):
                continue
            for b in (0, 1):
                e2 = e + b
                if 3 ** e2 > 2 ** (t + 1):
                    nxt[e2] = nxt.get(e2, 0) + cnt
        cur = nxt
        out.append(sum(cur.values()))
    return out


Hbits = -(cf * mpmath.log(cf, 2) + (1 - cf) * mpmath.log(1 - cf, 2))
lam = (1 - cf) / cf


def beta_f(th):
    return mpmath.log((1 - cf) * mpmath.sin(cf * th) / (cf * mpmath.sin((1 - cf) * th)))


def kappa_f(th):
    be = beta_f(th)
    return cf * mpmath.exp(be * (1 - cf)) * mpmath.cos((1 - cf) * th) + (1 - cf) * mpmath.exp(-be * cf) * mpmath.cos(cf * th)


worst1 = -mpf(10) ** 9
for m in range(2, 13):
    th = mpmath.pi / (m + 2)
    Ns = N_barrier(m, 200)
    for L in range(0, 201):
        rhs = mpmath.exp(abs(beta_f(th)) * (m - cf)) * mpmath.power(2, Hbits * L) * kappa_f(th) ** L
        assert Ns[L] <= rhs, (m, L, Ns[L], rhs)
        if Ns[L] > 0:
            worst1 = max(worst1, mpmath.log(Ns[L] / rhs))
check(True, f"Theorem 1 against exact N_m(L), m = 2..12, L <= 200 (max ln(N/bound) = {float(worst1):.3f})")

worst2 = mpf(10) ** 9
for M in range(4, 17):
    th = mpmath.pi / M
    Pc = M * mpmath.exp(abs(beta_f(th))) / (2 * cf * (1 - cf) * lam ** M)
    As = A_hard(M, 200)
    for L in range(1, 201):
        lhs_bound = mpmath.power(2, Hbits * L) * kappa_f(th) ** L / Pc
        assert As[L] >= lhs_bound, (M, L, As[L], lhs_bound)
        worst2 = min(worst2, mpmath.log(As[L] / lhs_bound))
check(True, f"Theorem 2 (crude constant) against exact A_M(L), M = 4..16, L <= 200 (min log-margin {float(worst2):.2f})")

maxratio = 0
for m in range(2, 11):
    Ns = N_barrier(m, 160)
    As = A_hard(m + 2, 160)
    for L in range(1, 161):
        if As[L] > 0:
            maxratio = max(maxratio, Ns[L] / As[L])
            assert Ns[L] <= (m + 2) ** 17 * As[L]
check(maxratio < 1.001, f"Corollary 3: N_m(L) <= (m+2)^17 A_(m+2)(L); observed max N_m/A_(m+2) = {maxratio:.7f} (m <= 10, L <= 160)")
