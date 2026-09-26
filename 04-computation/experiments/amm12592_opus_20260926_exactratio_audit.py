#!/usr/bin/env python3
"""Independent adversarial audit of the exact-ratio lemma of THM-4494 (AMM 12592, level-0 bottom regime of THM-4468).

Target: 01-canon/theorems/THM-4494-amm12592-exact-ratio-bottom-regime-c-below-log2-3.md and the theta step of
04-computation/experiments/amm12592_opus_20260926_exactratio_contours.py.

Notation (THM-4468, Lemma L0):  N = 16*4^k,  m = N/16,  K = 2N - ceil(cN),  A_0 = K + 2m,  R_0 = ceil(cN) - N - 1,
bottom regime 1 <= r < r_1 = 3N/128,  P_r = C(A_0 + r - 1, r)/C(R_0, r) = prod_{j<r} (A_0 + j)/(R_0 - j).
Lemma (THM-4494): with a = 2 - c + 1/8 >= A_0/N, b = c - 1 - 1/N_A <= R_0/N (N >= N_A) and
I(tau) = int_0^tau log((a+s)/(b-s)) ds:  (1) each factor increases in j and decreases in N; (2) log P_r <= N I(r/N);
(3) I convex, I(0) = 0, so I(tau) <= (tau/tau_1) I(tau_1) on (0, tau_1]; (4) P_r(N) <= P_r(N_A) for r <= r_A = 3N_A/128
and P_r(N) <= exp(r_A I(tau_1)/tau_1) for r_A < r < 3N/128.  The script takes
theta = max( max_{r<=r_A} P_r(N_A) [exact rationals at N_A], exp(r_A I(tau_1)/tau_1), exp(N_A I(tau_1)) ).

What this audit computes, for c in {158/100, 197/125, 203/128, 159/100} and N_A = 4096, all in exact rationals unless
stated: K, m, A_0, R_0, r_A, every P_r for r <= r_A and its maximum, the crude theta_1 of THM-4468, I(3/128) with
mpmath at 50 digits (primitive AND quadrature); then the two monotonicity claims (in j at N_A; in N for
N in {4096, 8192, 16384} and, along the admissible sequence, 65536 ... 2^26), the Riemann-sum direction, the convexity
chord, the uniform validity of the script's theta over N >= N_A, and a corrected theta built from the majorant
factors (aN + j)/(bN - j) (which ARE monotone in N).  Nothing here is imported from the audited scripts.
"""
from __future__ import annotations
import math
from fractions import Fraction
from math import comb
import mpmath as mp

mp.mp.dps = 50
NA = 4096
TAU1 = Fraction(3, 128)
CS = [(Fraction(158, 100), "158/100"), (Fraction(197, 125), "197/125"), (Fraction(203, 128), "203/128"),
      (Fraction(159, 100), "159/100 (= 1.59)")]
# eps_D certified by the derived contour script at N_A (from the .out files); negligible, quoted for the margins
EPS_D = {Fraction(158, 100): math.exp(-120.74), Fraction(197, 125): math.exp(-114.08)}


def ceil_frac(x: Fraction) -> int:
    return -((-x.numerator) // x.denominator)


def params(c: Fraction, N: int):
    cN = ceil_frac(c * N)
    K = 2 * N - cN
    m = N // 16
    return dict(cN=cN, K=K, m=m, A0=K + 2 * m, R0=cN - N - 1, delta=Fraction(cN) - c * N)


def factor(c: Fraction, N: int, j: int) -> Fraction:
    p = params(c, N)
    return Fraction(p["A0"] + j, p["R0"] - j)


def products(c: Fraction, N: int, rmax: int):
    """[P_1, ..., P_rmax] as exact Fractions."""
    out, prod = [], Fraction(1)
    for j in range(rmax):
        prod *= factor(c, N, j)
        out.append(prod)
    return out


def ab(c: Fraction, N_A: int = NA):
    return 2 - c + Fraction(1, 8), c - 1 - Fraction(1, N_A)


def majorant_factor(c: Fraction, N, j) -> Fraction:
    a, b = ab(c)
    return (a * N + j) / (b * N - j)


def majorant_products(c: Fraction, N: int, rmax: int):
    out, prod = [], Fraction(1)
    for j in range(rmax):
        prod *= majorant_factor(c, N, j)
        out.append(prod)
    return out


def mpf(x: Fraction):
    return mp.mpf(x.numerator) / x.denominator


def I_prim(c: Fraction, tau):
    a, b = ab(c)
    a, b, tau = mpf(a), mpf(b), mpf(tau)
    prim = lambda x: x * mp.log(x) - x
    return (prim(a + tau) - prim(a)) + (prim(b - tau) - prim(b))


def I_quad(c: Fraction, tau):
    a, b = ab(c)
    a, b, tau = mpf(a), mpf(b), mpf(tau)
    return mp.quad(lambda s: mp.log((a + s) / (b - s)), [0, tau])


def crude_theta1(c: Fraction) -> Fraction:
    return (2 - c + Fraction(1, 8) + TAU1) / (c - 1 - TAU1 - Fraction(1, NA))


def logP(P: Fraction):
    return mp.log(mpf(P))


def main():
    print("=== THM-4494 exact-ratio lemma: independent audit (exact rationals + mpmath 50 digits) ===")
    print(f"N_A = {NA}, tau_1 = 3/128, r_A = 3 N_A/128 = {3 * NA // 128}")
    print(f"crude threshold of THM-4468: theta_1 < 1 iff c > 203/128 = {203 / 128:.7f}; log_2 3 = {math.log2(3):.7f}")
    rA = 3 * NA // 128

    # ------------------------------------------------------------ 0. identity check P_r = C(A0+r-1, r)/C(R0, r)
    c = Fraction(197, 125)
    p = params(c, NA)
    P = products(c, NA, rA)
    for r in (1, 2, 3, 50, rA):
        exact = Fraction(comb(p["A0"] + r - 1, r), comb(p["R0"], r))
        assert exact == P[r - 1], ("product identity", r)
    print("\n[0] identity C(A0+r-1,r)/C(R0,r) == prod_{j<r}(A0+j)/(R0-j) checked at r = 1,2,3,50,r_A (c = 197/125): OK")

    # ------------------------------------------------------------ 1. per-c table at N_A
    print("\n[1] per-c quantities at N = N_A (exact; floats shown to 6-7 digits)")
    summary = {}
    for c, name in CS:
        p = params(c, NA)
        a, b = ab(c)
        P = products(c, NA, rA)
        Pmax = max(P)
        rmax = P.index(Pmax) + 1
        Pmin = min(P)
        rmin = P.index(Pmin) + 1
        th_crude = crude_theta1(c)
        I1 = I_prim(c, TAU1)
        I1q = I_quad(c, TAU1)
        assert abs(I1 - I1q) < mp.mpf(10) ** -40, "primitive vs quadrature disagree"
        chord = mp.exp(rA * I1 / mpf(TAU1))
        theta_r1 = mp.exp(NA * I1)
        theta_script = max(mpf(Pmax), chord, theta_r1)      # what the derived script uses (exact P_max -> 50 digits)
        # majorant products at N_A (valid uniform bound, see [4]); a/b is the r = 1 value
        G = majorant_products(c, NA, rA)
        Gmax = max(G)
        rGmax = G.index(Gmax) + 1
        theta_corr = max(mpf(Gmax), chord)
        R0min = (c - 1) * NA - 1
        epsD = EPS_D.get(c, 0.0)
        marg_script = float(R0min) * (1 - float(theta_script) - epsD)
        marg_corr = float(R0min) * (1 - float(theta_corr) - epsD)
        summary[c] = dict(theta_script=theta_script, theta_corr=theta_corr, Gmax=Gmax)
        print(f"\n  c = {name} = {float(c):.6f}")
        print(f"    ceil(cN_A) = {p['cN']} (delta = ceil(cN) - cN = {float(p['delta']):.6f}), K = {p['K']}, m = {p['m']}, "
              f"A_0 = K + 2m = {p['A0']}, R_0 = {p['R0']}, r_A = {rA}")
        print(f"    a = 2 - c + 1/8 = {a} = {float(a):.6f} (A_0/N_A = {float(Fraction(p['A0'], NA)):.6f} <= a: {Fraction(p['A0'], NA) <= a}); "
              f"b = c - 1 - 1/N_A = {b} = {float(b):.6f} (R_0/N_A = {float(Fraction(p['R0'], NA)):.6f} >= b: {Fraction(p['R0'], NA) >= b})")
        print(f"    exact P_r at N_A: P_1 = {p['A0']}/{p['R0']} = {float(P[0]):.7f}; max_{{r<=r_A}} P_r = {float(Pmax):.7f} at r = {rmax}; "
              f"min = {float(Pmin):.4e} at r = {rmin}; P_(r_A) = {float(P[-1]):.4e}")
        print(f"    crude theta_1 (THM-4468) = {float(th_crude):.6f}  [{'< 1' if th_crude < 1 else '>= 1: crude certificate fails'}]")
        print(f"    I(3/128) = {mp.nstr(I1, 12)} (primitive) / {mp.nstr(I1q, 12)} (quadrature)  [{'< 0' if I1 < 0 else '>= 0: lemma inapplicable'}]")
        print(f"    chord bound exp(r_A I/tau_1) = {mp.nstr(chord, 8)}; exp(N_A I) = {mp.nstr(theta_r1, 8)}")
        print(f"    script theta = max(exact P_max, chord, exp(N_A I)) = {mp.nstr(theta_script, 8)}")
        print(f"    majorant products at N_A: G_1 = a/b = {float(a / b):.7f}; max_{{r<=r_A}} G_r = {float(Gmax):.7f} at r = {rGmax}; G_(r_A) = {float(G[-1]):.4e}")
        print(f"    corrected uniform theta = max(G_max, chord) = {mp.nstr(theta_corr, 8)}")
        print(f"    bottom margin ((c-1)N_A - 1)(1 - theta - eps_D), eps_D = {epsD:.2e}: script theta -> {marg_script:.3f}; corrected theta -> {marg_corr:.3f} (need >= 5)")

    # ------------------------------------------------------------ 2. monotonicity in j (at N_A)
    print("\n[2] monotonicity of the exact factor (A_0+j)/(R_0-j) in j at N = N_A, j = 0..r_A-1")
    for c, name in CS:
        fs = [factor(c, NA, j) for j in range(rA)]
        inc = all(fs[j + 1] > fs[j] for j in range(rA - 1))
        p = params(c, NA)
        # crossing of 1: first j with factor >= 1
        cross = next((j for j in range(rA) if fs[j] >= 1), None)
        print(f"  c = {name}: strictly increasing in j: {inc}; factor_0 = {float(fs[0]):.6f}, factor_(r_A-1) = {float(fs[-1]):.6f}; "
              f"first j with factor >= 1: {cross} (analytic: j >= (R_0 - A_0)/2 = {(p['R0'] - p['A0']) / 2:.1f}); "
              f"denominators positive (r_A < R_0 = {p['R0']}): {rA < p['R0']}")
    print("  (general: numerator A_0 + j increases, denominator R_0 - j decreases and stays positive for j < r_1 < R_0 -> increasing; CONFIRMED)")

    # ------------------------------------------------------------ 3. monotonicity in N: exact vs majorant
    print("\n[3] monotonicity in N of the exact factor f_j(N) = (A_0(N)+j)/(R_0(N)-j) and of the majorant g_j(N) = (aN+j)/(bN-j)")
    print("    (claim 1 of the lemma says f_j decreases in N; the script uses it as P_r(N) <= P_r(N_A) for r <= r_A)")
    Ns_task = [4096, 8192, 16384]
    Ns_adm = [16 * 4 ** k for k in range(4, 12)]          # 4096 ... 2^26, the admissible sizes N = 16*4^k >= N_A
    for c, name in CS:
        print(f"\n  c = {name}")
        print(f"    {'N':>9} {'admissible':>10} {'delta':>8} {'A_0':>10} {'R_0':>10} {'f_0=P_1':>11} {'f_1':>11} {'f_95':>11} {'P_2':>11} {'P_96':>11} {'g_0=a/b':>11} {'g_95':>11}")
        base = None
        viol_exact, viol_maj = [], []
        for N in sorted(set(Ns_task + Ns_adm)):
            p = params(c, N)
            adm = (N % 16 == 0) and (N // 16 & (N // 16 - 1) == 0) and (int(round(math.log(N // 16, 4))) == math.log(N // 16, 4))
            f0, f1, f95 = factor(c, N, 0), factor(c, N, 1), factor(c, N, 95)
            Pn = products(c, N, 96)
            g0, g95 = majorant_factor(c, N, 0), majorant_factor(c, N, 95)
            if N == NA:
                base = dict(f=[factor(c, N, j) for j in range(rA)], g=[majorant_factor(c, N, j) for j in range(rA)], P=Pn)
            else:
                for j in range(rA):
                    if factor(c, N, j) > base["f"][j]:
                        viol_exact.append((N, j))
                    if majorant_factor(c, N, j) > base["g"][j]:
                        viol_maj.append((N, j))
            print(f"    {N:9d} {str(adm):>10} {float(p['delta']):8.4f} {p['A0']:10d} {p['R0']:10d} {float(f0):11.7f} {float(f1):11.7f} {float(f95):11.7f} "
                  f"{float(Pn[1]):11.7f} {float(Pn[95]):11.4e} {float(g0):11.7f} {float(g95):11.7f}")
        bad = sorted(set(N for N, _ in viol_exact))
        print(f"    exact factor: violations of f_j(N) <= f_j(N_A) at N in {bad} ({len(viol_exact)} (N,j) pairs)"
              f"{'; NONE' if not viol_exact else ''}")
        if viol_exact:
            worst = max(((N, j, factor(c, N, j) - base['f'][j]) for N, j in viol_exact), key=lambda t: t[2])
            print(f"      largest excess: N = {worst[0]}, j = {worst[1]}: f_j(N) - f_j(N_A) = {float(worst[2]):.3e}")
            v1 = [(N, float(products(c, N, 1)[0])) for N in bad if products(c, N, 1)[0] > base['P'][0]]
            print(f"      P_1(N) > P_1(N_A) = {float(base['P'][0]):.7f} at: {[(N, round(v, 7)) for N, v in v1]}")
        print(f"    majorant factor: violations of g_j(N) <= g_j(N_A): {len(viol_maj)} (expected 0: d/dN g_j = -(a+b) j/(bN-j)^2 <= 0, = 0 for j = 0)")
        # sup of P_1 over the admissible sequence
        sup_P1 = max((products(c, N, 1)[0], N) for N in Ns_adm)
        a, b = ab(c)
        print(f"    sup_(admissible N <= 2^26) P_1(N) = {float(sup_P1[0]):.7f} at N = {sup_P1[1]}; limit a/(c-1) = {float((2 - c + Fraction(1, 8)) / (c - 1)):.7f}; "
              f"uniform bound a/b = {float(a / b):.7f}; script's theta = {mp.nstr(summary[c]['theta_script'], 8)}")
        ts = summary[c]["theta_script"]
        exceeded = (mpf(sup_P1[0]) > ts) or (mpf((2 - c + Fraction(1, 8)) / (c - 1)) > ts)
        print(f"    => script theta is a valid uniform bound for P_1 over N >= N_A: "
              f"{'NO (exceeded by the exact P_1 at a larger admissible N and/or by the limit a/(c-1))' if exceeded else 'not contradicted'}; "
              f"corrected theta = {mp.nstr(summary[c]['theta_corr'], 8)} >= a/b >= P_1(N) for all N >= N_A: {summary[c]['theta_corr'] >= mpf(a / b)}")

    # ------------------------------------------------------------ 4. Riemann-sum direction: log P_r(N) <= N I(r/N)
    print("\n[4] Riemann-sum inequality log P_r(N) <= N I(r/N) (left sum of the increasing integrand log((a+s)/(b-s)))")
    for c, name in CS:
        if I_prim(c, TAU1) >= 0 and c != Fraction(203, 128):
            pass
        worst = None
        ok = True
        for N in (4096, 16384, 65536):
            Pn = products(c, N, rA)
            for r in range(1, rA + 1):
                lhs = logP(Pn[r - 1])
                rhs = N * I_prim(c, Fraction(r, N))
                gap = rhs - lhs
                if gap < 0:
                    ok = False
                if worst is None or gap < worst[0]:
                    worst = (gap, N, r)
        # also the pure majorant sum (no ceilings): sum_j log g_j(N) <= N I(r/N)
        okm = True
        for N in (4096, 16384):
            for r in range(1, rA + 1):
                s = sum(mp.log(mpf(majorant_factor(c, N, j))) for j in range(r))
                if s > N * I_prim(c, Fraction(r, N)):
                    okm = False
        print(f"  c = {name}: holds for all 1 <= r <= r_A at N in (4096, 16384, 65536): {ok}; smallest slack N I(r/N) - log P_r = "
              f"{mp.nstr(worst[0], 6)} at N = {worst[1]}, r = {worst[2]}; majorant sum <= N I(r/N): {okm}")
    print("  (direction: phi(s) = log((a+s)/(b-s)) is increasing, so (1/N) sum_{j<r} phi(j/N) <= int_0^{r/N} phi; CONFIRMED)")

    # ------------------------------------------------------------ 5. convexity and chord
    print("\n[5] convexity of I and the chord bound I(tau) <= (tau/tau_1) I(tau_1) on (0, tau_1]")
    for c, name in CS:
        a, b = ab(c)
        I1 = I_prim(c, TAU1)
        # I''(tau) = 1/(a+tau) + 1/(b-tau) > 0 on [0, tau_1] since a > 0, b - tau_1 > 0
        conv = (a > 0) and (b - TAU1 > 0)
        chord_ok = True
        worst = None
        for k in range(1, 97):
            tau = Fraction(k, 4096)
            lhs = I_prim(c, tau)
            rhs = mpf(tau) / mpf(TAU1) * I1
            if lhs > rhs + mp.mpf(10) ** -45:
                chord_ok = False
            if worst is None or rhs - lhs < worst[0]:
                worst = (rhs - lhs, tau)
        # per-r uniform consequence: for every r in [1, r_1): P_r <= exp(r I(tau_1)/tau_1) <= exp(I(tau_1)/tau_1)
        print(f"  c = {name}: I'' = 1/(a+tau) + 1/(b-tau) > 0 on [0, 3/128]: {conv} (a = {float(a):.4f}, b - tau_1 = {float(b - TAU1):.4f}); "
              f"chord holds at tau = k/4096, k = 1..96: {chord_ok} (min slack {mp.nstr(worst[0], 4)} at tau = {worst[1]}); "
              f"I(3/128) = {mp.nstr(I1, 8)}; exp(I(tau_1)/tau_1) [chord at r = 1, valid for ALL r >= 1 and ALL N >= N_A] = {mp.nstr(mp.exp(I1 / mpf(TAU1)), 8)}")

    # ------------------------------------------------------------ 6. c-threshold of the lemma: I(3/128; c) = 0
    print("\n[6] where the exact-ratio lemma stops: root of I(3/128; c) = 0 (N_A = 4096) by bisection")
    lo, hi = mp.mpf("1.5"), mp.mpf("1.7")
    f = lambda cc: (lambda a, b: ((a + t) * mp.log(a + t) - (a + t) - (a * mp.log(a) - a)) + ((b - t) * mp.log(b - t) - (b - t) - (b * mp.log(b) - b)))(
        2 - cc + mp.mpf(1) / 8, cc - 1 - mp.mpf(1) / NA)
    t = mpf(TAU1)
    assert f(lo) > 0 and f(hi) < 0
    for _ in range(200):
        mid = (lo + hi) / 2
        if f(mid) > 0:
            lo = mid
        else:
            hi = mid
    print(f"  I(3/128; c) < 0 iff c > {mp.nstr(hi, 12)}   (197/125 = 1.576 is above: {197 / 125 > float(hi)}; 63/40 = 1.575 is below: {63 / 40 < float(hi)})")
    print(f"  sanity: at c = 203/128 (crude threshold) I(3/128) = {mp.nstr(I_prim(Fraction(203, 128), TAU1), 6)} < 0 and the crude theta_1 = {float(crude_theta1(Fraction(203, 128))):.6f}")

    # ------------------------------------------------------------ 7. reproduction of the .out numbers
    print("\n[7] reproduction of the derived script's printed theta values")
    for c, name, printed in [(Fraction(197, 125), "197/125", 0.952946), (Fraction(158, 100), "158/100", 0.939789)]:
        ts = summary[c]["theta_script"]
        print(f"  c = {name}: audit recomputation of the script's theta = {mp.nstr(ts, 7)}; printed in .out: {printed}; agree to 6 decimals: {abs(float(ts) - printed) < 1e-6}")
    print("\n=== end of audit ===")


if __name__ == "__main__":
    main()
