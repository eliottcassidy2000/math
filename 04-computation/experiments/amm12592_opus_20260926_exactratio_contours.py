#!/usr/bin/env python3
"""EXACT-RATIO VARIANT (opus, collatz-exponent-atlas-20260926) of the THM-4468 contour certificate.

Differences from amm12592_procgen_20260923_hyp9128_contours.py (which is otherwise copied verbatim):
  * c is a command-line parameter (--cnum/--cden, default 158/100);
  * the bottom-regime constant theta_1 is the EXACT binomial ratio bound. Lemma L0(b) gives
    |[x^r] A| <= C(A0 + r - 1, r) with A0 = K + 2m, and the original certificate bounded
    C(A0+r-1, r)/C(R_0, r) by ((A0 + r - 1)/R_0)^r (largest factor to the r-th power), which
    forces c > 1.586. The exact ratio is P_r = prod_{j<r} (A0 + j)/(R_0 - j). With the rationals
    a = 2 - c + 1/8 >= A0/N and b = c - 1 - 1/N_A <= R_0/N (N >= N_A), P_r(N) <= G_r(N) :=
    prod_{j<r} (aN + j)/(bN - j); each MAJORANT factor is increasing in j and decreasing in N
    (d/dN = -j(a+b)/(bN-j)^2), so P_r(N) <= G_r(N_A) for r <= r_A = 3N_A/128, and
    log G_r(N) <= N I(r/N), I(tau) = int_0^tau log((a+s)/(b-s)) ds (left Riemann sum of an
    increasing integrand). I is convex with I(0) = 0, so once I(3/128) < 0 (certified by interval
    arithmetic) the chord gives P_r(N) <= exp(r_A I(3/128)/(3/128)) for r_A < r < 3N/128. The script
    uses theta_1 = max(max_{r<=r_A} G_r(N_A), exp(r_A I(3/128)/(3/128))), and max_r G_r(N_A) = G_1 = a/b.
    CORRECTION (independent audit, 2026-09-26): the first version used the EXACT products at N_A,
    claiming the exact factors (A0+j)/(R0-j) decrease in N; they do not (the residue ceil(cN) - cN
    varies along N = 16*4^k: at c = 197/125, P_1(16384) = 0.953057 > P_1(4096) = 0.952946, with
    limit a/(c-1) = 0.953125). The majorant factors are monotone; the bound is now uniform in N.
  Everything else (levels, eps_D, eps_T, middle regime) is the original interval certificate,
  re-run at the new c. A PASS means THM-4468's construction (kappa = 1/16, B = 4) works with
  deadlines t_i = min(4N, ceil(c(N+i))) for all N >= N_A, given the finite certificates below N_A.

HYP-9128, analytic part: interval-arithmetic certificates for the fold contour bounds.

Construction (see the note amm12592_procgen_20260923_hyp9128_proof.md):
  c = 159/100, super-block [N,4N) (N a power of two, N >= N_A), deadlines t_i = min(4N, ceil(c(N+i))),
  handoff state S_N(w) = (1+w)^m (1 + w^m + ... + w^{15m}), m = N/16, target
  F_0(u) = ((1+u)/2)^{N-1} S_N((1-u^2)/4),  Long's fold  F_{i+1} = (F_i - W_i)/y, y = (1-u)/2.

Rate function (N-uniform, rigorous):   w = (1-u^2)/4,
  log|F_0(u)| <= N*[a + b/16 + (15/16) max(0,o)] - a + log 16,
  a = log|(1+u)/2|, b = log|1+w|, o = log|w|       (uses |sum_{j<16} w^{mj}| <= 16 max(1,|w|)^{15m}).

Certified items (all by adaptive interval subdivision of circles, mpmath.iv, prec 60):
  V1  levels 1 <= i <= h: coefficient c^{(i)}_{R_{i-1}} and tail tau_i bounds  ->  M_i <= 1/2;
  V2D level-0 bottom-regime tail term eps_D;   V2T level-0 top-regime term eps_T;
  V2M level-0 middle regime: coefficient bounds B_t times the pairing factor g1^{min(t,R_0-t)}.
The script prints every cell (circle, certified exponent, constant) and the final inequalities at N = N_A,
which are monotone in N beyond N_A (checked symbolically in the printout).
Usage: python3 amm12592_procgen_20260923_hyp9128_contours.py [--NA 4096] [--quick]
"""
from __future__ import annotations
import argparse, math, sys, time
from fractions import Fraction
import numpy as np
import mpmath as mp

iv = mp.iv
iv.prec = 60
C = Fraction(159, 100)   # overwritten from the command line in main()
CF = float(C)
LOG16 = math.log(16.0)


# ----------------------------------------------------------------------------- float helpers (for choosing circles)
def f_parts(u):
    eps = 1e-300
    w = (1 - u ** 2) / 4
    a = np.log(np.abs((1 + u) / 2) + eps)
    b = np.log(np.abs(1 + w) + eps)
    o = np.log(np.abs(w) + eps)
    lam = np.log(np.abs(u) + eps)
    mu = np.log(np.abs((u - 1) / 2) + eps)
    return a, b, o, lam, mu


def f_rate(u):
    a, b, o, lam, mu = f_parts(u)
    return a + b / 16 + (15 / 16) * np.maximum(0, o), lam, mu


TH = np.linspace(0, math.pi, 1201)


def float_max(x0, rho, vbox, sbox, tie=False):
    """tie=True: sigma = (c-1)(1+v) + eta with eta in sbox (the level relation), else sigma in sbox."""
    u = x0 + rho * np.exp(1j * TH)
    G, lam, mu = f_rate(u)
    best = -np.inf
    for v in vbox:
        for s in sbox:
            sig = (CF - 1) * (1 + v) + s if tie else s
            best = max(best, np.max(G - sig * lam - v * mu))
    return best


def choose_circle(vbox, sbox, need01=True, need_m1=False, gap=0.0, tie=False):
    best = (np.inf, None)
    for x0 in np.linspace(-0.7, 1.7, 49):
        for rho in np.exp(np.linspace(np.log(0.05), np.log(9.0), 90)):
            if abs(x0) >= 0.97 * rho:
                continue
            if need01 and (abs(1 - x0) >= 0.9 * rho or rho - abs(1 - x0) < gap):
                continue
            if need_m1 and (abs(-1 - x0) >= 0.9 * rho or rho - abs(-1 - x0) < gap):
                continue
            val = float_max(x0, rho, vbox, sbox, tie)
            if val < best[0]:
                best = (val, (x0, rho))
    # local refinement
    val0, (x0, rho) = best
    step_x, step_r = 0.05, 0.05
    for _ in range(40):
        improved = False
        for dx, dr in [(step_x, 0), (-step_x, 0), (0, step_r), (0, -step_r)]:
            x1, r1 = x0 + dx, rho * math.exp(dr)
            if abs(x1) >= 0.97 * r1 or (need01 and (abs(1 - x1) >= 0.9 * r1 or r1 - abs(1 - x1) < gap)) \
                    or (need_m1 and (abs(-1 - x1) >= 0.9 * r1 or r1 - abs(-1 - x1) < gap)):
                continue
            v1 = float_max(x1, r1, vbox, sbox, tie)
            if v1 < val0:
                val0, x0, rho, improved = v1, x1, r1, True
        if not improved:
            step_x /= 2
            step_r /= 2
    return val0, (round(x0, 6), round(rho, 6))


# ----------------------------------------------------------------------------- interval evaluation on a circle arc
def iv_parts(x0, rho, T):
    X = iv.mpf(x0) + iv.mpf(rho) * iv.cos(T)
    Y = iv.mpf(rho) * iv.sin(T)
    X2, Y2 = X ** 2, Y ** 2
    XY2 = (2 * X * Y) ** 2
    a = iv.log(((1 + X) ** 2 + Y2) / 4) / 2
    b = iv.log(((5 - X2 + Y2) ** 2 + XY2) / 16) / 2
    o = iv.log(((1 - X2 + Y2) ** 2 + XY2) / 16) / 2
    lam = iv.log(X2 + Y2) / 2
    um1 = iv.log((X - 1) ** 2 + Y2) / 2          # log|u-1|
    up1 = iv.log((X + 1) ** 2 + Y2) / 2          # log|u+1|
    return a, b, o, lam, um1, up1


def pos(I):
    return iv.mpf([max(mp.mpf(0), mp.mpf(I.a)), max(mp.mpf(0), mp.mpf(I.b))])


def certify_circle(x0, rho, exponent_fn, extra_fn, target, NA, maxdepth=28, extra_cap=40.0):
    """Adaptive bisection of theta in [0, pi] (the integrands are conjugation symmetric).
    exponent_fn(parts) -> interval of the N-rate Psi; extra_fn(parts) -> interval of the O(1) correction E.
    On every final piece: Psi.b <= target (< 0) and E.b <= extra_cap.
    Returns (max_p Psi_p, max_p (NA*Psi_p + E_p), #pieces).  For N >= NA the arc bound
    max_p (N Psi_p + E_p) <= max_p (NA Psi_p + E_p) + (N - NA) max_p Psi_p  is therefore certified."""
    todo = [(mp.mpf(0), mp.mpf(mp.pi) + mp.mpf(2) ** -50, 0)]
    maxrate, maxQ, pieces = -mp.inf, -mp.inf, 0
    while todo:
        lo, hi, d = todo.pop()
        T = iv.mpf([lo, hi])
        P = iv_parts(x0, rho, T)
        R = exponent_fn(P)
        E = extra_fn(P)
        bad = (R.b > target) or (E.b > extra_cap)
        if bad and d < maxdepth:
            mid = (lo + hi) / 2
            todo += [(lo, mid, d + 1), (mid, hi, d + 1)]
            continue
        if bad:
            raise AssertionError(("cannot certify", x0, rho, float(lo), float(hi), float(R.b), float(E.b), target))
        maxrate = max(maxrate, mp.mpf(R.b))
        maxQ = max(maxQ, NA * mp.mpf(R.b) + mp.mpf(E.b))
        pieces += 1
    return float(maxrate), float(maxQ), pieces


def main():
    global C, CF
    ap = argparse.ArgumentParser()
    ap.add_argument("--NA", type=int, default=4096)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--cnum", type=int, default=158)
    ap.add_argument("--cden", type=int, default=100)
    args = ap.parse_args()
    C = Fraction(args.cnum, args.cden)
    CF = float(C)
    NA = args.NA
    t0 = time.time()
    c1 = CF - 1.0
    vmax = 4 / CF - 1 + 1.0 / NA
    print(f"=== HYP-9128 contour certificates: c = {C} = {CF}, kappa = 1/16, B = 4, N_A = {NA} ===")
    print(f"level range v = i/N in [0, {vmax:.6f}], sigma = (c-1)(1+v) + eta, eta in [-2/N_A, 1/N_A]")

    # ------------------------------------------------------------------ V1
    dv = 0.02 if args.quick else 0.01
    cells = []
    v = 0.0
    while v < vmax - 1e-12:
        va, vb = v, min(v + dv, vmax)
        cells.append((va, vb, -2.0 / NA, 1.0 / NA))
        v = vb
    worst_rate, worst_logM = -1e9, -1e9
    print(f"\n[V1] {len(cells)} v-cells (width {dv}); per cell: circle (x0, rho), certified max rate Psi, log prefactor")
    for (va, vb, ea, eb) in cells:
        fv, (x0, rho) = choose_circle([va, vb], [ea, eb], tie=True)
        target = fv + 0.004
        V = iv.mpf([va, vb])
        Et = iv.mpf([ea, eb])
        c1i = iv.mpf(C.numerator) / C.denominator - 1
        # sigma = (c-1)(1+v) + eta:  Psi = G - (c-1) lam - v((c-1) lam + mu) - eta lam
        rate_fn = lambda P: (P[0] + P[1] / 16 + iv.mpf(15) / 16 * pos(P[2]) - c1i * P[3]
                             - V * (c1i * P[3] + P[4] - iv.log(iv.mpf(2))) - Et * P[3])
        # correction for coefficient: -a - lam + log16 ; for tail: -a - log|u-1| + log16 ; take max of both
        def extra_fn(P):
            e1, e2 = (-P[0] - P[3]), (-P[0] - P[4])
            hi = max(mp.mpf(e1.b), mp.mpf(e2.b))
            return iv.mpf([hi, hi]) + iv.log(iv.mpf(16))
        assert target < 0
        rate, Q, pieces = certify_circle(x0, rho, rate_fn, extra_fn, target, NA)
        # M_i <= rho*(max|coef integrand| + max|tail integrand|) <= 2 rho exp(max_p (N Psi_p + E_p))
        logM = math.log(2 * rho) + Q
        extra = Q - NA * rate
        worst_rate = max(worst_rate, rate)
        worst_logM = max(worst_logM, logM)
        print(f"  v in [{va:.2f},{vb:.4f}]: circle ({x0:+.4f}, {rho:.4f}) rate <= {rate:+.5f} (float {fv:+.5f}) "
              f"extra <= {extra:.3f} pieces {pieces:4d}  log M_i bound at N_A: {logM:+.2f}", flush=True)
    print(f"  => worst certified rate {worst_rate:+.5f};  max_i M_i <= exp({worst_logM:.2f}) = {math.exp(worst_logM):.3e} at N = N_A "
          f"(bound decreasing in N since every rate < 0)")
    assert math.exp(worst_logM) <= 0.5

    # ------------------------------------------------------------------ level-0 constants
    r1frac = 3 / 128
    qbar = r1frac / (c1 - r1frac)                        # q_1 = r1/(R0 - r1 + 1) <= qbar
    pbar = r1frac / c1                                    # p = r/R0 >= pbar on the middle range
    g1 = 1 - 2 * pbar * (1 - pbar)
    # --- exact-ratio bottom bound (opus variant) ---
    a_ex = iv.mpf(2) - iv.mpf(C.numerator) / C.denominator + iv.mpf(1) / 8      # >= (K + 2m)/N
    b_ex = iv.mpf(C.numerator) / C.denominator - 1 - iv.mpf(1) / NA             # <= R_0/N
    tau1 = iv.mpf(3) / 128
    # I(tau) = int_0^tau log((a+s)/(b-s)) ds = [(a+s)log(a+s) - (a+s)] - [a log a - a] + [(b-s)log(b-s) - (b-s)]_0^tau
    def prim(x):
        return x * iv.log(x) - x
    I_tau1 = (prim(a_ex + tau1) - prim(a_ex)) + (prim(b_ex - tau1) - prim(b_ex))
    assert I_tau1.b < 0, f"I(3/128) not certified negative: {I_tau1}"
    # majorant running products at N = N_A for 1 <= r <= r_A = 3 N_A/128: with a >= A0/N and b <= R_0/N
    # for N >= N_A, P_r(N) <= G_r(N) := prod_{j<r} (aN + j)/(bN - j), and each majorant factor is
    # decreasing in N (d/dN = -j(a+b)/(bN-j)^2 <= 0), so P_r(N) <= G_r(N_A).  (The EXACT factors
    # (A0+j)/(R0-j) are NOT monotone in N because ceil(cN) - cN varies with N: error found by the
    # independent audit of 2026-09-26; the first version of this script used them.)  For
    # r_A < r < 3N/128 the convexity chord gives P_r(N) <= exp(N I(r/N)) <= exp(r I(tau1)/tau1)
    # <= exp(r_A I(tau1)/tau1).
    cN = -((-C.numerator * NA) // C.denominator)          # ceil(c N_A)
    K_A = 2 * NA - cN
    A0 = Fraction(K_A + 2 * (NA // 16))
    R0_A = Fraction(cN - NA - 1)
    rA = (3 * NA) // 128
    a_fr = 2 - C + Fraction(1, 8)
    b_fr = C - 1 - Fraction(1, NA)
    assert A0 <= a_fr * NA and R0_A >= b_fr * NA, "majorant constants must dominate at N_A"
    prod = Fraction(1)
    theta_A = Fraction(0)
    for j in range(rA):
        prod *= (a_fr * NA + j) / (b_fr * NA - j)
        theta_A = max(theta_A, prod)
    theta_chord = float(iv.exp(rA * I_tau1 / tau1).b)
    theta_one = float(theta_A) * (1 + 2 ** -40)            # rounded up
    theta_r1 = float(iv.exp(NA * I_tau1).b)
    theta1 = max(theta_one, theta_chord, theta_r1)
    print(f"[exact ratio] N_A = {NA}: K = {K_A}, A0 = K + 2m = {A0}, R_0 = {R0_A}, r_A = {rA}; a = {a_fr} >= A0/N_A, b = {b_fr} <= R_0/N_A; "
          f"max_r<=r_A majorant product G_r(N_A) = {theta_one:.6f} (G_1 = a/b = {float(a_fr / b_fr):.6f}; the exact P_1(N_A) = {float(A0 / R0_A):.6f} is not uniform in N); "
          f"G_(r_A) = {float(prod):.3e}; chord bound for r > r_A: {theta_chord:.3e}")
    theta1_crude = (2 - CF + 1 / 8 + r1frac) / (c1 - r1frac - 1.0 / NA)
    s0a, s0b = c1 - 1.0 / NA, c1
    print(f"\n[level 0] r1 = 3N/128;  qbar = {qbar:.5f}, pbar = {pbar:.5f}, g1 = {g1:.6f} (log g1 = {math.log(g1):.5f}), "
          f"theta1 majorant-ratio(N_A) = {theta1:.5f} (r<=r_A value {theta_one:.5f}, r=r1 value {theta_r1:.3e}; I(3/128) in [{float(I_tau1.a):+.6f},{float(I_tau1.b):+.6f}]; "
          f"crude theta1 would be {theta1_crude:.5f}), sigma0 in [{s0a:.6f}, {s0b:.6f}]")

    S0 = iv.mpf([s0a, s0b])
    qb = iv.mpf(qbar)

    # V2D: eps_D contour (encloses 0, 1; stays >= 0.25 away from 1 so that |lambda| qbar < 1)
    fv, (x0, rho) = choose_circle([0.0], [s0a, s0b], need01=True, gap=0.25)
    rate_fn = lambda P: P[0] + P[1] / 16 + iv.mpf(15) / 16 * pos(P[2]) - S0 * P[3]

    def extra_D(P):
        # -a + log16 + log( 2 qbar / (|z-1|^2 (1 - |lambda| qbar)) ),  |lambda| = |z+1|/|z-1|
        lamabs = iv.exp(P[5] - P[4])
        den = 1 - lamabs * qb
        if not den.a > 0:
            return iv.mpf([mp.inf, mp.inf])
        return -P[0] + iv.log(iv.mpf(16)) + iv.log(2 * qb) - 2 * P[4] - iv.log(den)
    rD, QD, pD = certify_circle(x0, rho, rate_fn, extra_D, fv + 0.004, NA)
    eD = QD - NA * rD
    logepsD = math.log(rho) + QD
    print(f"[V2D] D-contour ({x0:+.4f},{rho:.4f}): rate <= {rD:+.5f}, extra <= {eD:.3f}, eps_D <= exp({logepsD:.2f})")

    # V2T: eps_T contour (encloses -1, 0, 1)
    fv, (x0, rho) = choose_circle([0.0], [s0a, s0b], need01=True, need_m1=True, gap=0.25)

    def extra_T(P):
        lamabs = iv.exp(P[4] - P[5])         # |z-1|/|z+1|
        den = 1 - lamabs * qb
        if not den.a > 0:
            return iv.mpf([mp.inf, mp.inf])
        return -P[0] + iv.log(iv.mpf(16)) + iv.log(iv.mpf(2)) - P[4] - P[5] - iv.log(den)
    rT, QT, pT = certify_circle(x0, rho, rate_fn, extra_T, fv + 0.004, NA)
    eT = QT - NA * rT
    logepsT = math.log(rho) + QT
    print(f"[V2T] top-contour ({x0:+.4f},{rho:.4f}): rate <= {rT:+.5f}, extra <= {eT:.3f}, eps_T <= exp({logepsT:.2f})")

    # bottom regime margin
    R0min = c1 * NA - 1
    marg_bottom = R0min * (1 - theta1 - math.exp(logepsD))
    print(f"[bottom] margin >= ((c-1)N_A - 1)(1 - theta1 - eps_D) = {marg_bottom:.3f} (need >= 5; increasing in N)")
    assert marg_bottom >= 5
    marg_top = R0min * (1 - math.exp(logepsT))
    print(f"[top]    margin >= ((c-1)N_A - 1)(1 - eps_T) = {marg_top:.3f} (need >= 5); top boundary |u_(0,R0)| <= eps_T = {math.exp(logepsT):.2e} <= 1")
    assert marg_top >= 5 and math.exp(logepsT) <= 1

    # V2M: middle regime.  For r in [r1, R0 - r1]:
    #   |e_{0,r}| <= sum_t |w_t| |K_t(r)|,   |K_t(r)| <= sqrt(2 R0) g1^{min(t,R0-t)} binom(R0,r)   (pairing + binomial lower bound)
    #   sum_t |w_t| g1^{min} <= sum_cells (N dt + 1) max_{t in cell} B_t g1^{min}  (+ tail sum tau_0)
    dt = 0.01 if args.quick else 0.005
    logg1 = math.log(g1)
    worst_mid_rate = -1e9
    tau = 0.0
    tcells = []
    while tau < s0b - 1e-12:
        tcells.append((tau, min(tau + dt, s0b)))
        tau += dt
    print(f"\n[V2M] {len(tcells)} tau-cells of width {dt} for B_t (t = tau N), plus the tail sum tau_0 (t = R_0)")
    cell_terms = []      # (log value at N_A of (N dt+1) rho exp(...), exponent per N)
    for (ta, tb) in tcells:
        pair = max(min(ta, s0a - tb), 0.0)
        fv, (x0, rho) = choose_circle([0.0], [ta, tb], need01=False)
        Tb = iv.mpf([ta, tb])
        rate_fn_t = lambda P: P[0] + P[1] / 16 + iv.mpf(15) / 16 * pos(P[2]) - Tb * P[3]
        extra_fn_t = lambda P: -P[0] - P[3] + iv.log(iv.mpf(16))
        r_, Q_, p_ = certify_circle(x0, rho, rate_fn_t, extra_fn_t, fv + 0.004, NA)
        expo = r_ + pair * logg1
        assert expo < 0, "middle-regime cell exponent not negative"
        worst_mid_rate = max(worst_mid_rate, expo)
        val = math.log(NA * (tb - ta) + 1) + math.log(rho) + Q_ + NA * pair * logg1
        cell_terms.append(val)
        print(f"  tau in [{ta:.3f},{tb:.3f}]: circle ({x0:+.4f},{rho:.4f}) rate <= {r_:+.5f}, pairing {pair:.4f}*log g1, "
              f"exponent {expo:+.5f}; log[(N dt+1) B g1^min] <= {val:+.2f}", flush=True)
    fv, (x0, rho) = choose_circle([0.0], [s0a, s0b], need01=True)
    extra_tail = lambda P: -P[0] - P[4] + iv.log(iv.mpf(16))
    r_, Q_, p_ = certify_circle(x0, rho, rate_fn, extra_tail, fv + 0.004, NA)
    worst_mid_rate = max(worst_mid_rate, r_)
    cell_terms.append(math.log(rho) + Q_)
    print(f"  tail sum tau_0: circle ({x0:+.4f},{rho:.4f}) rate <= {r_:+.5f}, log bound {math.log(rho) + Q_:+.2f}")
    mx = max(cell_terms)
    logsum = mx + math.log(sum(math.exp(t - mx) for t in cell_terms))
    R0max = c1 * NA
    mid_total = 0.5 * math.log(2 * R0max) + logsum
    print(f"  => sqrt(2 R0) * sum <= exp({mid_total:.2f}) = {math.exp(mid_total):.3e} (need <= 1/2, then margin >= binom(R0,r)/2 >= 5)")
    assert math.exp(mid_total) <= 0.5
    print(f"  middle-regime worst exponent {worst_mid_rate:+.5f}; sqrt(N)(N dt+1) exp(N*exponent) decreasing for N >= {1.5/abs(worst_mid_rate):.0f} (<= N_A: {1.5/abs(worst_mid_rate) <= NA})")
    assert 1.5 / abs(worst_mid_rate) <= NA
    print(f"\nAll analytic inequalities certified at N = N_A = {NA}; each bound is of the form poly(N) exp(-delta N) with "
          f"delta >= {min(-worst_rate, -rD, -rT):.4f} (levels, tails) and the middle-regime exponent < 0, hence holds for all N >= N_A.")
    print(f"[{time.time() - t0:.0f}s]")
    return 0


if __name__ == "__main__":
    sys.exit(main())
