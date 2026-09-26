#!/usr/bin/env python3
"""collatz_dipspectrum_20260926_orders_audit.py -- independent adversarial audit of THM-4498 (Theorem 5 and
section 1.2b of collatz_dipspectrum_20260926_entropy_curve.md) and of the post-audit addendum of the thin-divergence
note (Lemma 1.4b and section 1.6b of collatz_thin_20260925_thin_divergent_orbits.md, the UPDATE line of THM-4476).
Auditor: separate agent (Fable), 2026-09-26, worktree collatz-exponent-atlas-20260926.  Nothing here is a proof;
the sections re-derive the constants of the written proofs, test every inequality that can be tested exactly, and
re-implement the two controls of collatz_dipspectrum_20260926_orders.py with EXACT integer arithmetic.

Sections of the output (numbered as the audit claims):
 (0) constants: alpha, rho*, h*, lambda*, a*, theta_0, theta_1, c_0.
 (1) section 1.2b: the geometric ratio bound, its threshold t_0(rho), the constant and its blow-up at log_4 3,
     Stirling's upper bound, and an exact check of the block bound sum_(o >= o_0) C(t,o) <= bound(t) for t <= 3000.
 (2) Theorem 5(a): the block-end inequality of step (vii) (checked at n = 2^t), the margin with the stated K,
     the Hoeffding exponent, exact hypergeometric tails against the cited bound, the constants table (3) of the
     original script, Stirling's lower bound, |h'| <= 1 on [rho, rho + 1/t'], the exact fraction of V with
     final rise F(v) <= delta, and an integer check of the sufficient condition S_i(w) >= -ct + 2 (both sheets).
 (3) Theorem 5(c): the delta' >= 29 sum (exact and factorised), K >= 40, binomial symmetry in place of the CLT,
     and an integer check of the multiplicative carry at gamma = 1/2 < log_2(3/2) with tight margins (both sheets).
 (4) Remark (i): D/X at the 0.7925 column of the brute-force table.
 (5) controls: EXACT integer DP for W_t(gamma) (thresholds 3^(100 o) > 2^(100 j - (100-p) t) with gamma = p/100;
     exact ties at o = 0 excluded by the strict inequality), the float DP of the original script mirrored, every
     tie flip listed, the normalised table (1) recomputed and compared with the filed .out, local exponents,
     extension to t = 10000, exact binomial tail sums, the "13.8" formula; the brute-force table (2) recomputed.
 (6) Lemma 1.4b: every constant, and sum_(o >= o_0) C(k,o) <= 12.8 2^(k h(rho)) k^(-1/2) on a grid of theta, k.
 (7) section 1.6b: the concavity constants, inequality (iv), monotonicity, well-foundedness, the size of X_3.
 (8) provenance: sha256 of collatz_dipspectrum_20260926_orders.py / .out (raw bytes, CR count).
Usage: python3 collatz_dipspectrum_20260926_orders_audit.py [OUTPATH] [--quick]
"""
import hashlib
import math
import os
import random
import re
import sys
import time

from mpmath import mp, mpf, ceil as mpceil, floor as mpfloor, log as mplog

mp.dps = 50
ALPHA_MP = mplog(3) / mplog(2)
ALPHA = math.log2(3)
L32 = math.log2(1.5)          # 0.585 of the notes
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
QUICK = '--quick' in sys.argv
ARGS = [a for a in sys.argv[1:] if not a.startswith('--')]
OUTPATH = ARGS[0] if ARGS else os.path.join(REPO, '05-knowledge', 'results', 'collatz_dipspectrum_20260926_orders_audit.out')
OUT = open(OUTPATH, 'w', encoding='utf-8', newline='\n')
T0 = time.time()


def say(s=''):
    print(s)
    OUT.write(s + '\n')
    OUT.flush()


def h(p):
    if p <= 0 or p >= 1:
        return 0.0
    return -(p * math.log2(p) + (1 - p) * math.log2(1 - p))


def hprime(p):
    return math.log2((1 - p) / p)


def E(gamma):
    return h(max(0.5, gamma / ALPHA))


def log2int(n):
    return math.log2(n) if n > 0 else float('-inf')


# ---------------------------------------------------------------- (0) constants
rho_star = 1 / ALPHA
h_star = h(rho_star)
lam_star = math.log2(rho_star / (1 - rho_star)) / ALPHA
lam_star_alt = math.log(1 / L32) / math.log(3)
a_star = lam_star / h_star - 0.5
theta0 = 1 - ALPHA / 2
theta1 = theta0 / 2
c0 = theta0
rho_min = (1 - theta1) / ALPHA

say("== (0) constants ==")
say("   alpha = log_2 3 = %.12f   log_2(3/2) = %.9f   log_4 3 = %.9f" % (ALPHA, L32, ALPHA / 2))
say("   rho* = log_3 2 = %.9f   h* = h(rho*) = %.9f   (note: 0.9499555)" % (rho_star, h_star))
say("   lambda* = log_2(rho*/(1-rho*))/alpha = %.9f = log_3(1/log_2(3/2)) = %.9f = -h'(rho*)/alpha = %.9f   (note: 0.488077)"
    % (lam_star, lam_star_alt, -hprime(rho_star) / ALPHA))
say("   a* = lambda*/h* - 1/2 = %.9f   (note: 0.013789);  lambda*/h* = %.6f (note bounds it by 0.52)" % (a_star, lam_star / h_star))
say("   theta_0 = 1 - alpha/2 = c_0 = %.9f   theta_1 = theta_0/2 = %.9f   rho_min = (1-theta_1)/alpha = %.9f   (note says rho >= 0.5655: %s)"
    % (theta0, theta1, rho_min, "TRUE" if rho_min >= 0.5655 else "FALSE, it is %.6f >= 0.5654 only" % rho_min))

# ---------------------------------------------------------------- (1) section 1.2b
say("")
say("== (1) section 1.2b: geometric tail above o_0 = ceil(rho t - 2), ratio (t-o)/(o+1) <= (1 - rho + 2/t)/(rho - 2/t) ==")
say("   gamma     rho      r=(1-rho)/rho  1/(1-r)=rho/(2rho-1)   t_0 (ratio<1 needs t > 4/(2rho-1))   asymptotic block constant (rho/(1-rho))^2/((1-r) sqrt(2 pi rho(1-rho)))")
for g in (0.7930, 0.80, 0.82, 0.85, 0.88, 0.91, 0.94, 0.97, 1.0):
    rho = g / ALPHA
    r = (1 - rho) / rho
    t0 = math.floor(4 / (2 * rho - 1)) + 1
    Cb = (rho / (1 - rho)) ** 2 / ((1 - r) * math.sqrt(2 * math.pi * rho * (1 - rho)))
    say("   %.4f   %.6f   %.6f      %10.3f          %6d            %10.3f" % (g, rho, r, 1 / (1 - r), t0, Cb))
say("   -> the constant is uniform on compact subsets of (log_4 3, 1] (t_0 and 1/(1-r) are monotone in rho) and blows up like rho/(2 rho - 1) as gamma -> log_4 3.")
say("   exact check of the block bound: S(t) = sum_(o >= o_0) C(t, o) versus B(t) = (rho/(1-rho))^2 2^(t h(rho)) / (delta_t sqrt(2 pi t p (1-p))), p = o_0/t, delta_t = 1 - (1-rho+2/t)/(rho-2/t)")
say("   also the chain's pieces: max_(o>=o_0) (t-o)/(o+1) <= (1-rho+2/t)/(rho-2/t) [ratio ok?], C(t,o_0) <= 2^(t h(p))/sqrt(2 pi t p(1-p)) [Stirling ok?], 2^(t h(p)) <= (rho/(1-rho))^2 2^(t h(rho)) [h' ok?]")
for g in (0.82, 0.91, 0.97, 1.0):
    rho = g / ALPHA
    for t in ((100, 300, 1000, 3000) if not QUICK else (100, 300)):
        o0 = int(mpceil(mpf(g) / ALPHA_MP * t - 2))
        S = sum(math.comb(t, o) for o in range(o0, t + 1))
        p = o0 / t
        rt = (1 - rho + 2 / t) / (rho - 2 / t)
        delta_t = 1 - rt
        if delta_t <= 0:
            say("   gamma=%.2f t=%5d o_0=%5d: t < t_0(rho) = %d, the ratio bound is %.4f >= 1 and 1.2b gives NO bound for this block (the block is covered by the finite n_0 term)"
                % (g, t, o0, math.floor(4 / (2 * rho - 1)) + 1, rt))
            continue
        logB = 2 * math.log2(rho / (1 - rho)) + t * h(rho) - math.log2(delta_t) - 0.5 * math.log2(2 * math.pi * t * p * (1 - p))
        ratio_ok = (t - o0) / (o0 + 1) <= rt
        stirling_ok = log2int(math.comb(t, o0)) <= t * h(p) - 0.5 * math.log2(2 * math.pi * t * p * (1 - p))
        hp_ok = t * h(p) <= 2 * math.log2(rho / (1 - rho)) + t * h(rho)
        say("   gamma=%.2f t=%5d o_0=%5d p=%.5f  S/B = %.4f  ratio ok: %s  Stirling ok: %s  h' ok: %s" % (g, t, o0, p, 2 ** (log2int(S) - logB), ratio_ok, stirling_ok, hp_ok))

# Stirling upper and lower bounds, all k, for a range of n
say("   Stirling bounds used: 2^(n h(p))/sqrt(8 n p(1-p)) <= C(n,k) <= 2^(n h(p))/sqrt(2 pi n p(1-p)), p = k/n, 1 <= k <= n-1:")
worst_up = -1e9
worst_lo = 1e9
ns = list(range(2, 401)) + [500, 1000, 2000, 3000]
for n in ns:
    for k in range(1, n):
        p = k / n
        lc = log2int(math.comb(n, k))
        up = n * h(p) - 0.5 * math.log2(2 * math.pi * n * p * (1 - p))
        lo = n * h(p) - 0.5 * math.log2(8 * n * p * (1 - p))
        worst_up = max(worst_up, lc - up)
        worst_lo = min(worst_lo, lc - lo)
say("   max over n in [2,400] u {500,1000,2000,3000}, all k, of log2 C(n,k) - log2(upper) = %.3e (must be <= 0: %s); min of log2 C - log2(lower) = %.3e (must be >= 0: %s)"
    % (worst_up, worst_up <= 1e-9, worst_lo, worst_lo >= -1e-9))
say("   [time %.1fs]" % (time.time() - T0))

# ---------------------------------------------------------------- (2) Theorem 5(a)
say("")
say("== (2) Theorem 5(a) lower bound ==")
say("   (i) drift: rho alpha - 1 = gamma - 1 = -c: identity.")
say("   (ii) E Y_i = i j/t' with j = ceil(rho t') <= rho t' + 1, so E Y_i <= rho i + i/t' <= rho i + 1: identity.")
say("   (iii) Hoeffding (1963, Thm 4 -> Thm 2, values in [0,1]): P(mean dev >= x) <= exp(-2 i x^2); x = c/(2 alpha) gives exp(-c^2 i/(2 alpha^2)): matches the note.")
say("   (vii) block end: the note claims 4 * 2^((gamma-1)(t+1)) >= 4 n^(gamma-1) for n in [2^t, 2^(t+1)). With gamma - 1 = -c < 0, n^(gamma-1) is DEcreasing in n,")
say("        so its maximum over the block is at n = 2^t, namely 2^(-ct) > 2^(-c(t+1)). The claimed inequality is FALSE for every n in the block:")
for g in (0.82, 0.85, 0.88, 0.91, 0.94, 0.97):
    c = 1 - g
    t = 100
    lhs = 4 * 2 ** (-c * (t + 1))
    rhs_at_2t = 4 * 2 ** (-c * t)
    say("        gamma=%.2f t=%d: 4*2^((gamma-1)(t+1)) = %.6e < 4 n^(gamma-1) at n = 2^t = %.6e  (ratio 2^(-c) = %.4f)" % (g, t, lhs, rhs_at_2t, lhs / rhs_at_2t))
say("        What is needed for M_i n >= 4 n^gamma on the whole block is S_i(w) >= -ct + 2 (n = 2^t is the worst case).")
say("        The margin actually available from the stated K, (c + 0.585) K >= delta_0 + c + 2, is S_(K+l)(w) >= (c+0.585)K - ct - delta_0 >= -ct + c + 2 > -ct + 2:")
say("        the '+c' in the K condition is exactly the slack that repairs the step (or, without it, the weaker M_i n >= 4 * 2^(-c) n^gamma >= 3.5 n^gamma still beats the carry).")
say("   (viii) carry: 4 n^gamma - (3/2)^t >= n^gamma iff (3/2)^t <= 3 n^gamma; with n >= 2^t this holds for EVERY t when gamma >= log_2(3/2) (no t_0 needed).")
say("   (vi) margin algebra: 0.585 K - c t' - delta_0 = (c + 0.585) K - c t - delta_0; the note's target -c(t+1) + 2 needs (c+0.585)K >= delta_0 - c + 2, weaker than the stated condition.")
say("   constants table (3) of the original script recomputed (i_0 = min{i : q^i/(1-q) <= 1/2}, q = exp(-c^2/(2 alpha^2))):")
K_of = {}
for g in (0.82, 0.85, 0.88, 0.91, 0.94, 0.97):
    c = 1 - g
    q = math.exp(-c * c / (2 * ALPHA * ALPHA))
    i0 = 1
    while q ** i0 / (1 - q) > 0.5:
        i0 += 1
    d0 = max(ALPHA, L32 * i0)
    K0 = math.ceil((d0 + c + 2) / (c + L32))
    Kneed = math.ceil((d0 + 2) / (c + L32))
    K_of[g] = (i0, d0, K0)
    margin = (c + L32) * K0 - d0   # S_(K+l)(w) + c t >= this
    say("      gamma=%.2f c=%.2f i_0=%6d delta_0=%9.2f K_0=%6d (K needed for -ct+2: %6d); with K_0: S_(K+l) + ct >= %.3f >= 2 + c = %.3f: %s; cost 2^(-K_0 E) = %.2e"
        % (g, c, i0, d0, K0, Kneed, margin, 2 + c, margin >= 2 + c, 2 ** (-K0 * E(g))))
say("   Hoeffding without replacement, exact hypergeometric tails: P(Y_i - E Y_i >= c i/(2 alpha)) / exp(-c^2 i/(2 alpha^2)), Y_i ~ Hyp(t', j, i), j = ceil(rho t'):")
for g in (0.82, 0.91):
    c = 1 - g
    tp = 400 if not QUICK else 120
    j = int(mpceil(mpf(g) / ALPHA_MP * tp))
    worst = 0.0
    worst_i = 0
    denom_cache = {}
    for i in range(1, tp + 1):
        mean = i * j / tp
        thr = mean + c * i / (2 * ALPHA)
        y0 = math.floor(thr) + 1  # smallest integer > thr; ties (thr integer) counted as >=: use ceil
        if abs(thr - round(thr)) < 1e-12:
            y0 = int(round(thr))
        num = sum(math.comb(j, y) * math.comb(tp - j, i - y) for y in range(max(y0, 0), min(i, j) + 1))
        P = num / math.comb(tp, i)
        bound = math.exp(-c * c * i / (2 * ALPHA * ALPHA))
        if P / bound > worst:
            worst, worst_i = P / bound, i
    say("      gamma=%.2f t'=%d j=%d: max_i P/bound = %.4f at i = %d (must be <= 1: %s)" % (g, tp, j, worst, worst_i, worst <= 1))
say("   |h'(p)| = |log_2((1-p)/p)| on [rho, rho + 1/t']: at rho* + 1/28 = %.6f it is %.5f (<= 1: %s); at rho* + 1/27 it is %.5f (> 1): the '|h'| <= 1' step needs t' >= 28."
    % (rho_star + 1 / 28, abs(hprime(rho_star + 1 / 28)), abs(hprime(rho_star + 1 / 28)) <= 1, abs(hprime(rho_star + 1 / 27))))


def frac_final_rise(tp, j, delta):
    """exact count of words of length tp with j ones whose final rise max_i (alpha Y_i - i) is <= delta
    (Y_i = ones among the last i letters): DP on the reversed word with the barrier Y_i <= (i + delta)/alpha.
    delta is an mpf; the 1e-30 nudge only matters at exact ties (i + delta)/alpha in Z, which the <= admits."""
    cap = [int(mpfloor((i + delta) / ALPHA_MP + mpf('1e-30'))) for i in range(tp + 1)]
    row = [1]  # index y = ones among the first i reversed letters
    for i in range(1, tp + 1):
        top = min(i, cap[i], j)
        ext = row + [0]
        new = [ext[y] + (ext[y - 1] if y >= 1 else 0) for y in range(top + 1)]
        row = new
    cnt = row[j] if j < len(row) else 0
    return cnt, math.comb(tp, j)


say("   exact fraction of V (words of length t' with j = ceil(rho t') ones) whose final rise F(v) is <= delta:")
for g in (0.82, 0.91):
    tp = 1000 if not QUICK else 300
    j = int(mpceil(mpf(g) / ALPHA_MP * tp))
    i0 = K_of[g][0]
    d0 = K_of[g][1]
    parts = []
    for delta in (ALPHA_MP, 2 * ALPHA_MP, mpf(5), mpf(10), mpf(20), mpf(50), mpf(i0) * (ALPHA_MP - 1)):
        cnt, tot = frac_final_rise(tp, j, delta)
        parts.append("F<=%.2f: %.4f" % (float(delta), cnt / tot))
    say("      gamma=%.2f t'=%d j=%d:  " % (g, tp, j) + "  ".join(parts) + "   (the note needs >= 1/2 at delta_0 = %.1f: trivially true since 0.585 t' < delta_0)" % d0)


def residue_from_word(word, b):
    x = 0
    for k in range(len(word)):
        for cand in (x, x + (1 << k)):
            y = cand
            ok = True
            for i in range(k + 1):
                if (y & 1) != word[i]:
                    ok = False
                    break
                y = y // 2 if y % 2 == 0 else (3 * y + b) // 2
            if ok:
                x = cand
                break
        else:
            raise RuntimeError("Terras bijection failed")
    return x


def orbit_nodip(n, t, b, gp, gq):
    """T_b^i(n)^gq >= n^gp for 0 <= i <= t, i.e. T^i(n) >= n^(gp/gq); returns (ok, first_dip_index)."""
    y = n
    npow = n ** gp
    for i in range(t + 1):
        if y ** gq < npow:
            return False, i
        y = y // 2 if y % 2 == 0 else (3 * y + b) // 2
    return True, -1


def sampled_words(t, lo, hi, nsamp, seed):
    rnd = random.Random(seed)
    res = []
    for _ in range(nsamp):
        w = [rnd.randint(0, 1) for _ in range(t)]
        o = 0
        mn = 1e9
        for i, x in enumerate(w, 1):
            o += x
            mn = min(mn, o * ALPHA - i)
        if lo <= mn <= hi:
            res.append(w)
    return res


say("   integer check of the sufficient condition (a): words w of length t with min_i S_i(w) in a tight window, S_i = o_i alpha - i; n = the representative in [2^t, 2^(t+1)); no dip = T_b^i(n) >= n^gamma for i <= t:")
t = 40
g = 0.82
c = 1 - g
gp, gq = 41, 50
nsamp = 60000 if not QUICK else 10000
for (lo, hi, label) in ((-c * t + 2, -c * t + 2.6, "corrected margin  [-ct+2, -ct+2.6]"), (-c * (t + 1) + 2, -c * t + 2 - 1e-9, "the note's margin [-c(t+1)+2, -ct+2)"), (-c * t - 3, -c * t - 2, "below the barrier [-ct-3, -ct-2] (M_i n <= n^gamma/4: must dip)")):
    ws = sampled_words(t, lo, hi, nsamp, 7)
    res = {}
    for b in (1, -1):
        ok_cnt = 0
        for w in ws:
            x = residue_from_word(w, b)
            n = x + (1 << t)
            ok, _ = orbit_nodip(n, t, b, gp, gq)
            ok_cnt += ok
        res[b] = ok_cnt
    say("      gamma=0.82 t=40 %s: %d words; no-dip on plus sheet %d/%d, minus sheet %d/%d" % (label, len(ws), res[1], len(ws), res[-1], len(ws)))
say("   [time %.1fs]" % (time.time() - T0))

# ---------------------------------------------------------------- (3) Theorem 5(c)
say("")
say("== (3) Theorem 5(c) ==")
say("   (i) uniform drift: alpha/2 - 1 = %.6f = -c_0 = -(1 - log_4 3): identity." % (ALPHA / 2 - 1))
say("   (ii) reduction: S_(K+l)(w) = 0.585 K + S_l(v) > 0.585 K - c_0 t' - delta' = 0.585K - c_0 t' - (c_0+0.585)K + 2 = -c_0 (t' + K) + 2 = -c_0 t + 2 >= -c t + 2 because c >= c_0 and t > 0: correct direction, and the target -ct + 2 is the RIGHT block end (unlike (a)).")
say("   (iii) Hoeffding, independent letters: P(alpha Y_i - i >= delta'/2) = P(Y_i/i - 1/2 >= (c_0 i + delta'/2)/(alpha i)) <= exp(-2 (c_0 i + delta'/2)^2/(alpha^2 i)) <= e^(-2 c_0 delta'/alpha^2) e^(-2 c_0^2 i/alpha^2).")
q = math.exp(-2 * c0 * c0 / ALPHA / ALPHA)


def hoeff_sum_exact(dp):
    s = 0.0
    i = 1
    while True:
        term = math.exp(-2 * (c0 * i + dp / 2) ** 2 / (ALPHA * ALPHA * i))
        s += term
        if term < 1e-30 and i > 100:
            break
        i += 1
    return s


def hoeff_sum_fact(dp):
    return math.exp(-2 * c0 * dp / ALPHA / ALPHA) * q / (1 - q)


say("      q = e^(-2 c_0^2/alpha^2) = %.6f, sum_(i>=1) q^i = q/(1-q) = %.4f, e^(-2 c_0 delta'/alpha^2) at delta' = 29: %.6f" % (q, q / (1 - q), math.exp(-2 * c0 * 29 / ALPHA / ALPHA)))
say("      factorised bound at delta' = 29: %.5f (<= 1/4: %s); exact sum at delta' = 29: %.5f" % (hoeff_sum_fact(29), hoeff_sum_fact(29) <= 0.25, hoeff_sum_exact(29)))
dmin_f = 1
while hoeff_sum_fact(dmin_f) > 0.25:
    dmin_f += 1
dmin_e = 1
while hoeff_sum_exact(dmin_e) > 0.25:
    dmin_e += 1
say("      smallest integer delta' with the factorised bound <= 1/4: %d; with the exact sum <= 1/4: %d" % (dmin_f, dmin_e))
Kmin = math.ceil((29 + 2) / (c0 + L32))
say("      delta' = (c_0 + 0.585) K - 2 >= 29 iff K >= %.3f, i.e. K >= %d (note: K >= 40: %s)" % (31 / (c0 + L32), Kmin, Kmin <= 40))
say("   (iv) endpoint: P(S_(t')(v) >= E S_(t')(v)) = P(Bin(t',1/2) >= t'/2) >= 1/2 EXACTLY for every t' (binomial symmetry), so the CLT is not needed and the good fraction is >= 1 - 1/4 - 1/2 = 1/4 for all t' >= 1:")
worst = 1.0
for tp in range(1, 401):
    P = sum(math.comb(tp, y) for y in range(math.ceil(tp / 2), tp + 1)) / 2 ** tp
    worst = min(worst, P)
say("      min over t' <= 400 of P(Bin(t',1/2) >= t'/2) = %.6f (>= 1/2: %s)" % (worst, worst >= 0.5))
say("   (v) multiplicative carry below log_2(3/2): integer check at gamma = 1/2, t = 30 (additive bound useless: (3/2)^30 = %.2e > n^(1/2) = 2^15 = %d), words with min_i S_i in [-ct+2, -ct+2.6] = [-13, -12.4]:" % (1.5 ** 30, 2 ** 15))
t = 30
g = 0.5
c = 0.5
ws = sampled_words(t, -c * t + 2, -c * t + 2.6, 60000 if not QUICK else 10000, 11)
for b in (1, -1):
    ok_cnt = 0
    first_fail = None
    for w in ws:
        x = residue_from_word(w, b)
        n = x + (1 << t)
        ok, idx = orbit_nodip(n, t, b, 1, 2)
        ok_cnt += ok
        if not ok and first_fail is None:
            first_fail = (n, idx)
    say("      sheet b=%+d: %d/%d representatives have no dip below n^(1/2) within 30 steps%s" % (b, ok_cnt, len(ws), "" if first_fail is None else "; FIRST FAILURE n=%d at i=%d" % first_fail))
ws2 = sampled_words(t, -c * t - 3, -c * t - 2, 60000 if not QUICK else 10000, 12)
for b in (1, -1):
    ok_cnt = sum(orbit_nodip(residue_from_word(w, b) + (1 << t), t, b, 1, 2)[0] for w in ws2)
    say("      control, min S_i in [-ct-3, -ct-2] (M_i n <= n^gamma/4), sheet b=%+d: %d/%d have no dip (expected 0)" % (b, ok_cnt, len(ws2)))
say("   [time %.1fs]" % (time.time() - T0))

# ---------------------------------------------------------------- (4) Remark (i)
say("")
say("== (4) Remark (i): at Korec's endpoint (c) with c = c_0 gives Theta(X); the brute-force column '0.792' (gamma = 0.7925, plus sheet) normalised by X ==")
bf_path = os.path.join(REPO, '05-knowledge', 'results', 'collatz_dipspectrum_20260926.out')
bf_txt = open(bf_path, encoding='utf-8', errors='replace').read().split('== sheet 3n-1')[0]
bf_gs = [float(x) for x in next(l for l in bf_txt.split('\n') if l.strip().startswith('gamma:')).split(':')[1].split()]
bf_rows = [(int(m.group(1)), [int(x) for x in m.group(2).split()]) for m in re.finditer(r"X=2\^(\d+) D:\s+([0-9 ]+)", bf_txt)]
say("   T:    " + "  ".join("%7d" % T for T, _ in bf_rows))
say("   D/X:  " + "  ".join("%7.4f" % (D[0] / 2 ** T) for T, D in bf_rows) + "   (the script does not count n = 1; gamma = 0.7925 is 0.00002 above log_4 3, where h(rho) - 1 = %.2e)" % (h(0.7925 / ALPHA) - 1))

# ---------------------------------------------------------------- (5) controls
say("")
say("== (5) exact DP for W_t(gamma) = #{w in {0,1}^t : o_j alpha - j > (gamma - 1) t, 1 <= j <= t}, gamma = p/100 ==")
say("   exact criterion: o alpha - j > (p/100 - 1) t  <=>  3^(100 o) > 2^E, E = 100 j - (100 - p) t; for E > 0 no tie is possible (alpha irrational),")
say("   E < 0 admits every o >= 0, and E = 0 is an EXACT TIE at o = 0 (1 > 1 is false, so o = 0 is excluded by the strict inequality).")
GAMMAS = [0.82, 0.85, 0.88, 0.91, 0.94, 0.97, 1.0]
PS = [82, 85, 88, 91, 94, 97, 100]
TS_ALL = [25, 50, 100, 150, 200, 300, 400, 500, 600]
TS_BIG = [800, 1000, 1500, 2000, 3000]
TS_EXT = [4000, 6000, 8000, 10000]
if QUICK:
    TS_BIG = [800]
    TS_EXT = []
TMAX = max(TS_ALL + TS_BIG + TS_EXT)


def make_B(tmax):
    """B[o] = floor(100 o alpha) exactly = bit_length(3^(100 o)) - 1."""
    B = [0] * (tmax + 1)
    x = 1
    m = 3 ** 100
    for o in range(1, tmax + 1):
        x *= m
        B[o] = x.bit_length() - 1
    return B


B_TAB = make_B(TMAX)


def omin_exact(t, p):
    res = []
    ptr = 0
    for j in range(1, t + 1):
        Ej = 100 * j - (100 - p) * t
        if Ej < 0:
            om = 0
        elif Ej == 0:
            om = 1
        else:
            while B_TAB[ptr] < Ej:
                ptr += 1
            om = ptr
        res.append(om)
    return res


def omin_float(t, gamma):
    """mirrors the original script's comparison o*ALPHA - j > (gamma - 1.0)*t in double precision."""
    thr = (gamma - 1.0) * t
    res = []
    o = 0
    for j in range(1, t + 1):
        while not (o * ALPHA - j > thr):
            o += 1
        res.append(o)
    return res


def W_from_thresholds(t, om):
    row = [1]
    lo = 0
    for j in range(1, t + 1):
        m = om[j - 1]
        if m > j or not row:
            return 0
        A = row[m - lo:] + [0]
        Bv = ([0] if m == lo else []) + row[(m - 1 - lo if m > lo else 0):]
        row = [a + b for a, b in zip(A, Bv)]
        lo = m
    return sum(row)


# reproduce the original script's float DP for a cross-check of the mirror
import importlib.util
spec = importlib.util.spec_from_file_location("orders_orig", os.path.join(HERE, 'collatz_dipspectrum_20260926_orders.py'))
orig = importlib.util.module_from_spec(spec)
spec.loader.exec_module(orig)

# parse the filed table (1)
orders_out = open(os.path.join(HERE, 'collatz_dipspectrum_20260926_orders.out'), encoding='utf-8').read()
filed = {}
for line in orders_out.split('\n'):
    m = re.match(r"\s+t=\s*(\d+):\s+(.*)$", line)
    if m:
        t = int(m.group(1))
        nums = [float(x) for x in re.findall(r"[0-9]+\.[0-9]+", m.group(2).split('(')[0])]
        filed[t] = nums

say("   tie flips: (t, gamma) pairs where j = (100-p) t/100 is an integer and the float threshold admits o = 0 at that j:")
flips = []
exactW = {}
floatW = {}
for t in TS_ALL + TS_BIG:
    gl = GAMMAS if t <= 600 else [0.82, 0.91, 1.0]
    for g in gl:
        p = PS[GAMMAS.index(g)]
        oe = omin_exact(t, p)
        of = omin_float(t, g)
        diff = [(j + 1, oe[j], of[j]) for j in range(t) if oe[j] != of[j]]
        We = W_from_thresholds(t, oe)
        Wf = W_from_thresholds(t, of)
        exactW[(t, g)] = We
        floatW[(t, g)] = Wf
        if t <= 600:
            Wo = orig.W_dp(t, g)
            if Wo != Wf:
                say("      MIRROR MISMATCH: the original script's W_dp(%d, %s) = %d differs from the threshold mirror %d" % (t, g, Wo, Wf))
        ties = [j + 1 for j in range(t) if 100 * (j + 1) - (100 - p) * t == 0]
        for (j, a, b) in diff:
            flips.append((t, g, j, a, b))
        if diff:
            say("      t=%5d gamma=%.2f: tie at j=%s; threshold (j, exact, float) = %s; W_float - W_exact = 10^%.3f (relative %.2e)"
                % (t, g, ties, diff, math.log10(Wf - We), (Wf - We) / We))
        elif ties:
            say("      t=%5d gamma=%.2f: tie at j=%s but the float comparison happened to reject o = 0 there (no flip)" % (t, g, ties))
if not flips:
    say("      none")
say("   non-tie safety: min over 1 <= o <= %d of the distance of 100 o alpha to the nearest integer:" % TMAX)
mind = 1.0
mino = 0
for o in range(1, TMAX + 1):
    v = 100 * o * ALPHA_MP
    d = float(abs(v - mpfloor(v + mpf('0.5'))))
    if d < mind:
        mind, mino = d, o
say("      %.3e at o = %d, i.e. |o alpha - (rational with denominator 100)| >= %.3e, far above the double-precision error (~1e-12): non-tie comparisons cannot flip." % (mind, mino, mind / 100))

say("   table (1) recomputed with the EXACT DP: W_t(gamma) t^(1/2) 2^(-tE) [gamma = 1: t^(3/2)]; second number = float DP (= filed); third = filed value from the .out; flag if exact differs from filed by > 0.00005")
say("   gamma:        " + "  ".join("%7.3f" % g for g in GAMMAS))
say("   E(gamma):     " + "  ".join("%7.4f" % E(g) for g in GAMMAS))
norm_exact = {}
for t in TS_ALL + TS_BIG:
    gl = GAMMAS if t <= 600 else [0.82, 0.91, 1.0]
    cells = []
    for gi, g in enumerate(gl):
        pw = 1.5 if g == 1.0 else 0.5
        ne = 2 ** (log2int(exactW[(t, g)]) + pw * math.log2(t) - t * E(g))
        nf = 2 ** (log2int(floatW[(t, g)]) + pw * math.log2(t) - t * E(g))
        norm_exact[(t, g)] = ne
        fv = filed.get(t, [None] * 7)[gi] if t in filed and gi < len(filed[t]) else None
        flag = ""
        if fv is not None and abs(ne - fv) > 0.00005:
            flag = "!"
        cells.append("%.4f/%.4f/%s%s" % (ne, nf, ("%.4f" % fv) if fv is not None else "n/a", flag))
    say("   t=%5d: " % t + "  ".join(cells))
say("   local exponent beta(t1,t2) = -[log2(W_(t2) 2^(-t2 E)) - log2(W_(t1) 2^(-t1 E))]/log2(t2/t1)  (Theta(t^(-1/2)) means beta -> 0.5; the -3/2 ballot order means 1.5):")
seq = TS_ALL + TS_BIG
for g in GAMMAS:
    ts = [t for t in seq if (t, g) in exactW]
    parts = []
    for a, b in zip(ts, ts[1:]):
        beta = -((log2int(exactW[(b, g)]) - b * E(g)) - (log2int(exactW[(a, g)]) - a * E(g))) / math.log2(b / a)
        parts.append("%d->%d: %.3f" % (a, b, beta))
    say("      gamma=%.2f: " % g + "  ".join(parts))

if TS_EXT:
    say("   extension of the exact DP to t = %s for gamma in {0.82, 0.91, 0.94, 0.97}: normalised W_t t^(1/2) 2^(-tE), local exponent beta from the previous t, and the exact tail Sigma_t = sum_(o > rho t) C(t,o) normalised the same way, with W_t/Sigma_t:" % TS_EXT)
    for g in (0.82, 0.91, 0.94, 0.97):
        p = PS[GAMMAS.index(g)]
        rho = g / ALPHA
        r = (1 - rho) / rho
        bound_W = 1 / ((1 - r) * math.sqrt(2 * math.pi * rho * (1 - rho)))
        bound_note = (rho / (1 - rho)) ** 2 * bound_W
        say("      gamma=%.2f rho=%.5f r=%.5f: asymptotic bound for the normalised tail 1/((1-r) sqrt(2 pi rho(1-rho))) = %.3f; with the extra (rho/(1-rho))^2 of the D_b count (the '13.8' formula) = %.3f"
            % (g, rho, r, bound_W, bound_note))
        prev = None
        for t in [600, 1000, 2000, 3000] + TS_EXT:
            if (t, g) not in exactW:
                exactW[(t, g)] = W_from_thresholds(t, omin_exact(t, p))
            W = exactW[(t, g)]
            o0 = int(mpfloor(mpf(g) / ALPHA_MP * t)) + 1
            Sig = sum(math.comb(t, o) for o in range(o0, t + 1))
            rt = (t - o0) / (o0 + 1)
            pt = o0 / t
            bound_t = 1 / ((1 - rt) * math.sqrt(2 * math.pi * pt * (1 - pt))) * 2 ** (t * (h(pt) - h(rho)))
            nW = 2 ** (log2int(W) + 0.5 * math.log2(t) - t * h(rho))
            nS = 2 ** (log2int(Sig) + 0.5 * math.log2(t) - t * h(rho))
            beta = "" if prev is None else "beta=%.3f" % (-((log2int(W) - t * h(rho)) - (log2int(prev[1]) - prev[0] * h(rho))) / math.log2(t / prev[0]))
            say("         t=%6d: W_t norm = %8.4f  %-12s  Sigma_t norm = %8.4f (finite-t bound %.3f)  W_t/Sigma_t = %.4f   [%.0fs]" % (t, nW, beta, nS, bound_t, W / Sig, time.time() - T0))
            prev = (t, W)

say("   table (2), THM-4487's brute-force D_+(2^T, gamma) normalised: D T^(1/2)/2^(T E) recomputed from collatz_dipspectrum_20260926.out (plus sheet):")
say("    gamma:   " + "  ".join("%7.3f" % g for g in bf_gs))
filed2 = {}
for line in orders_out.split('\n'):
    m = re.match(r"\s+T=(\d+):\s+(.*)$", line)
    if m:
        filed2[int(m.group(1))] = [float(x) for x in m.group(2).split()]
maxdiff = 0.0
for T, D in bf_rows:
    vals = [2 ** (math.log2(d) + 0.5 * math.log2(T) - T * E(g)) for g, d in zip(bf_gs, D)]
    if T in filed2:
        maxdiff = max(maxdiff, max(abs(a - b) for a, b in zip(vals, filed2[T])))
    say("    T=%2d:    " % T + "  ".join("%7.4f" % v for v in vals))
say("    max |recomputed - filed| over table (2) = %.1e" % maxdiff)
say("    local exponent beta(T, T+2) of D 2^(-TE) for the brute-force counts (window-length sawtooth included):")
for gi, g in enumerate(bf_gs):
    parts = []
    for (T1, D1), (T2, D2) in zip(bf_rows, bf_rows[1:]):
        beta = -((math.log2(D2[gi]) - T2 * E(g)) - (math.log2(D1[gi]) - T1 * E(g))) / math.log2(T2 / T1)
        parts.append("%d->%d: %.2f" % (T1, T2, beta))
    say("      gamma=%.3f: " % g + "  ".join(parts))
say("   [time %.1fs]" % (time.time() - T0))

# ---------------------------------------------------------------- (6) Lemma 1.4b
say("")
say("== (6) Lemma 1.4b: theta in (0, theta_1], rho = (1-theta)/alpha in [%.6f, %.6f), k >= 200, o_0 = ceil(rho k - 2) ==" % (rho_min, rho_star))
say("   (rho - 1/2) k >= 13 for k >= 200: worst (rho_min - 1/2) * 200 = %.4f (>= 13: %s)" % ((rho_min - 0.5) * 200, (rho_min - 0.5) * 200 >= 13))
worst_ratio = (1 - rho_min + 2 / 200) / (rho_min - 2 / 200)
say("   ratio bound (1 - rho + 2/k)/(rho - 2/k) <= 0.81: worst (rho_min, k=200) = %.5f (%s); 1/(1-0.81) = %.4f <= 5.3 (%s)" % (worst_ratio, worst_ratio <= 0.81, 1 / 0.19, 1 / 0.19 <= 5.3))
say("   p(1-p) >= rho*(1-rho*) = %.6f (note: 0.2329; 1/sqrt(2 pi 0.232858) = %.5f, note uses 0.827: safe %s)" % (rho_star * (1 - rho_star), 1 / math.sqrt(2 * math.pi * rho_star * (1 - rho_star)), 0.827 >= 1 / math.sqrt(2 * math.pi * rho_star * (1 - rho_star))))
sup_hh = (rho_star / (1 - rho_star)) ** 2
say("   (rho/(1-rho))^2 < (rho*/(1-rho*))^2 = %.5f  (note says <= 2.92: %s -- the supremum is 2.9224, so the rounded constant should read 2.93)" % (sup_hh, sup_hh <= 2.92))
say("   product of the note's rounded constants 5.3 * 2.92 * 0.827 = %.4f; product of the true suprema (1/0.19) * 2.92243 * 0.82673 = %.4f; both <= 12.8: %s"
    % (5.3 * 2.92 * 0.827, (1 / 0.19) * sup_hh / math.sqrt(2 * math.pi * rho_star * (1 - rho_star)), (1 / 0.19) * sup_hh / math.sqrt(2 * math.pi * rho_star * (1 - rho_star)) <= 12.8))
say("   representatives: 2^k > X/2 so a class mod 2^k meets [1, X] in at most 2 integers: factor 2.  k >= log_2 X - 1 >= 0.995 log_2 X iff log_2 X >= 200: given.")
say("   final constant: 2 * 12.8 / sqrt(0.995) = %.4f <= 27: %s" % (2 * 12.8 / math.sqrt(0.995), 2 * 12.8 / math.sqrt(0.995) <= 27))
say("   numerical verification: R(k, theta) = sum_(o >= o_0) C(k, o) * k^(1/2) / 2^(k h(rho)) on theta in {1e-9, 1e-6, 1e-4, 1e-3, 0.005, 0.01, ..., 0.10, theta_1} and every k in [200, 3000] (QUICK: step 25):")
thetas = [1e-9, 1e-6, 1e-4, 1e-3, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, theta1]
rhos = [(1 - mpf(th)) / ALPHA_MP for th in thetas]
hr = [h(float(rr)) for rr in rhos]
maxR = 0.0
argmax = None
maxR_theta = {th: 0.0 for th in thetas}
ks = range(200, 3001, 25 if QUICK else 1)
for k in ks:
    o_lo = int(0.565 * k) - 3
    row = [math.comb(k, o_lo)]
    for o in range(o_lo, k):
        row.append(row[-1] * (k - o) // (o + 1))
    suf = [0] * (len(row) + 1)
    for i in range(len(row) - 1, -1, -1):
        suf[i] = suf[i + 1] + row[i]
    for ti, th in enumerate(thetas):
        o0 = int(mpceil(rhos[ti] * k - 2))
        S = suf[o0 - o_lo]
        R = 2 ** (log2int(S) + 0.5 * math.log2(k) - k * hr[ti])
        if R > maxR_theta[th]:
            maxR_theta[th] = R
        if R > maxR:
            maxR, argmax = R, (k, th)
say("      max R over the grid = %.4f at (k, theta) = %s  (<= 12.8: %s; <= 13.5 for the two-representative form 2 S <= 27 2^(k h) k^(-1/2): %s)" % (maxR, argmax, maxR <= 12.8, maxR <= 13.5))
say("      max R per theta: " + "  ".join("%.0e:%.3f" % (th, maxR_theta[th]) if th < 1e-3 else "%.4g:%.3f" % (th, maxR_theta[th]) for th in thetas))
say("   [time %.1fs]" % (time.time() - T0))

# ---------------------------------------------------------------- (7) section 1.6b
say("")
say("== (7) section 1.6b ==")
say("   (i) concavity: h(rho) <= h* + h'(rho*)(rho - rho*) with rho - rho* = -theta/alpha and h'(rho*) = log_2((1-rho*)/rho*) = %.6f = -alpha lambda* (%.6f): so h(rho) <= h* + lambda* theta."
    % (hprime(rho_star), -ALPHA * lam_star))
worst_gap = -1.0
for i in range(1, 10001):
    th = theta1 * i / 10000
    rho = (1 - th) / ALPHA
    worst_gap = max(worst_gap, h(rho) - (h_star + lam_star * th))
say("      max over theta in (0, theta_1] (10^4 grid) of h(rho) - h* - lambda* theta = %.3e (must be <= 0: %s)" % (worst_gap, worst_gap <= 1e-12))
say("   (ii) (R'): 1.5 gives N(X) <= k + k N(X^(1-theta)) + #F_b(X, theta) with (E) <= k <= L = log_2 X; 1.4b bounds #F_b for theta <= theta_1 and X >= 2^200 uniformly in theta: (R') follows.")
say("   (iii) X^(1 - theta_X) = X 2^(-theta_X L) = X 2^(-c_1 log_2 L) = X L^(-c_1); X^(h(rho)) <= X^(h* + lambda* theta_X) = X^(h*) L^(lambda* c_1): identities.")
say("   (iv) lambda* c_1 - 1/2 = a* + eta lambda*/h* = a* + %.5f eta; need <= a - eta = a* + 3 eta (eta = (a-a*)/4): true since %.5f <= 3. Checked on a grid of a:" % (lam_star / h_star, lam_star / h_star))
for a in (0.0138, 0.014, 0.02, 0.05, 0.1, 0.5, 1.0):
    if a <= a_star:
        say("      a=%.4f <= a*: outside the theorem" % a)
        continue
    eta = (a - a_star) / 4
    c1 = (1 + eta) / h_star
    say("      a=%.4f eta=%.6f c_1=%.6f  lambda* c_1 - 1/2 = %.6f  a - eta = %.6f  ok: %s" % (a, eta, c1, lam_star * c1 - 0.5, a - eta, lam_star * c1 - 0.5 <= a - eta))
say("   (v) the step uses only log_2 Y <= L for Y = X L^(-c_1) (no monotonicity of Y -> Y^(h*) (log_2 Y)^a is needed; it is monotone anyway: d/dY[h* ln Y + a ln log_2 Y] = h*/Y + a/(Y ln Y) > 0 for Y > 1);")
say("       2 L^(-eta) + L^(-a) <= 1 and 27 <= K close the induction; K = max(X_3, 27) covers the base 2 <= X < X_3 since N(X) <= X.")
say("   (vi) well-foundedness: N(X) depends only on floor(X), and Y = X L^(-c_1) <= X / 200^(c_1) for X >= 2^200: at eta -> 0, 200^(1/h*) = %.1f, so floor(Y) < floor(X) and the induction on floor(X) is well-founded;"
    % (200 ** (1 / h_star)))
say("       the note says 'strong induction' over real X without this remark: a presentational gap, not an error.")
say("   size of X_3 (2 L^(-eta) + L^(-a) <= 1): log_2 log_2 X_3 for a = 0.02, 0.05, 0.1:")
for a in (0.02, 0.05, 0.1):
    eta = (a - a_star) / 4
    lo, hi = 1.0, 1e6
    # solve in log2 L: f(L) = 2 L^-eta + L^-a - 1 = 0; increase log2 L until <= 0
    l2 = 1.0
    while 2 * (2 ** l2) ** (-eta) + (2 ** l2) ** (-a) > 1:
        l2 *= 1.01
    say("      a=%.3f eta=%.5f: log_2 L >= %.1f, i.e. X_3 ~ 2^(2^%.0f)" % (a, eta, l2, l2))
say("   (vii) Corollary 4: the recursion in 1.7 is N_A(X) <= k N_A(X^(1-theta)) + #F_b + #F_(-b) + (k+1)(|b|/3+1) with base N_A <= 2X: (R) with a doubled no-dip term and an extra O(k) term, so 27 -> 54 and K adjusts; the polylog argument transfers, but '(R) verbatim' is not literally accurate.")

# ---------------------------------------------------------------- (8) provenance
say("")
say("== (8) provenance: sha256 of the THM-4498 script and output (raw bytes) ==")
for fn, claimed in (('collatz_dipspectrum_20260926_orders.py', '87f31934e76b6c958038ce46b4d26cd1de111e602c9ec08418d97275d8109225'),
                    ('collatz_dipspectrum_20260926_orders.out', '270c93827ae4e7e2da91bcc49e43505093664d6d850540b6fd9bee3727b2873f')):
    data = open(os.path.join(HERE, fn), 'rb').read()
    raw = hashlib.sha256(data).hexdigest()
    lf = hashlib.sha256(data.replace(b'\r\n', b'\n')).hexdigest()
    say("   %s: raw sha256 %s (header: %s, match %s); CR bytes %d; LF-normalised %s" % (fn, raw, claimed, raw == claimed, data.count(b'\r'), lf))
say("")
say("total time %.1fs" % (time.time() - T0))
OUT.close()
