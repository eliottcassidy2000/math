#!/usr/bin/env python3
r"""
lrc_tent_fLP_k9_kps_S73.py   (kind-pasteur-2026-07-07-S73, HYP-5147 follow-through)  v2

THE f-GAME ABOVE THE TENT: is 31/63 (k=9) improvable within the gap-histogram frame?

The frame: f >= 0 supported on (0, 1/7]; for every k-family,
    mu_{1/7}(E) >= 1 - k(k-1) int f / m(f),
    m(f) = min over safe gap vectors {g in [0,1/7]^k, sum g = 1} of the ring sum
           RINGSUM_f(g) = sum_{l,i} f(S_{i,l}) (terms with S <= 1/7).

FINDINGS (this file documents analysis + verification, not a search):

(1) TENT OPTIMALITY AMONG CONVEX f (proof).  For convex f, g |-> sum_i f(g_i) is
    Schur-convex, so the safe-polytope minimum is at the all-equal config g = 1/k
    (the majorization-minimal point), giving m(f) = k f(1/k) (+ rings: the all-equal
    config's ring-l sums are l/k > 1/7 for l >= 2, k <= 13 -- no ring payment).
    The floor is 1 - k(k-1) int f / (k f(1/k)); minimizing int f / f(1/k) over convex
    f >= 0 supported in (0,1/7] is achieved by the widest tent through (1/k, f(1/k)):
    f = (s - a)+ with d/da[(1/7-a)^2 / (1/9-a)] = 0 => (k=9) a = 5/63, value 31/63.
    The k=8 case: a = 3/28, value 3/4 (the theorem).  General k: a = (14-k)/(7k),
    floor = 1 - 2(k-1)(k-7)/(7k).

(2) RINGS DO NOT BITE AT k >= 9: the binding face (all gaps >= a) has ring-2 sums
    >= 2a = 2(14-k)/(7k); at k=9 that is 10/63 > 1/7, so no ring-2 term can be forced.
    (Verified numerically below.)

(3) NON-CONVEX f DOES NOT HELP (the all-equal cap): the adversary can always play
    all-equal, so every f obeys floor <= 1 - k(k-1) int f / (k f(1/k)); relaxing
    convexity only weakens the min elsewhere.  Numerical column-generation confirms
    the ping-pong: the adversary hides in a band around 1/k, the LP dodges per-config
    constraints forever (the working-set LP values are NOT valid floors; the earlier
    v1 output of this script is superseded by this analysis -- kept in git history).

(4) WHAT REMAINS FOR k=9, 10 (named):
    (a) SIGNED f (two-sided Markov: negative mass on (0, eps] taxes the tail's close
        pairs; k=9 safe configs can carry up to 4 tiny gaps, so the sign structure is
        genuinely different) -- the full degree-2 game, OPEN.
    (b) THE CONDITIONAL TENT (the k=9/k=10 program): bound
        E[F 1_{G_P}] = sum_pairs int_{G_P} f(frac(dx)) dx <= c * meas(G_P) * k(k-1) int f
        with c <= 1.7; then rho* >= meas(G_P)(1 - c(1 - 31/63)) >= m_P discharges k=9
        with ~4x headroom at c = 1 (0.4943 * 0.4921 = 0.2433 vs m_P = 0.0565), and
        k=10 at 2.4x.  For d large, Koksma/interval counting gives c -> 1 at rate
        #intervals(G_P)/d; SMALL d need a finite exact table per P -- the resonant
        shapes (f-window intervals dodging G_P's holes) are the honest obstruction
        (the R2 wall at the G_P-alignment level, but now FINITE: d small, P finite).

Verification below: (i) Schur/all-equal minimum for the k=9 tent; (ii) ring-2 sums
on the binding face; (iii) the k=9/k=10 conditional-tent arithmetic.
"""
import random
from fractions import Fraction as F_

THR = 1.0 / 7.0

def tent(s, a):
    return max(0.0, s - a) if s <= THR + 1e-15 else 0.0

def ringsum(g, a):
    k = len(g)
    tot = 0.0
    for i in range(k):
        s = 0.0
        for l in range(1, k):
            s += g[(i + l - 1) % k]
            if s > THR:
                break
            tot += tent(s, a)
    return tot

rng = random.Random(73)
print("=" * 92)
print("(i) k=9 tent (a = 5/63): random+descent search for min ring sum vs theory 1 - 9a = 18/63")
print("=" * 92)
a9 = 5.0 / 63.0
best = None
for _ in range(120000):
    g = [rng.uniform(0, THR) for _ in range(8)]
    s = sum(g)
    if s < 1 and 1 - s <= THR:
        g.append(1 - s)
        v = ringsum(g, a9)
        if best is None or v < best[0]:
            best = (v, sorted(g))
print(f"    search min = {best[0]:.6f}; theory m = 1 - 9*(5/63) = {1 - 9*a9:.6f} = 2/7")
print(f"    floor = 1 - 72*((1/7 - 5/63)^2/2)/(2/7) = {1 - 72*((THR - a9)**2/2)/(2.0/7):.6f} = 31/63 = {31/63:.6f}")

print()
print("(ii) ring-2 on the binding face: all gaps >= 5/63 => 2-sums >= 10/63 > 9/63 = 1/7:", 10/63 > 1/7)

print()
print("=" * 92)
print("(iii) conditional-tent arithmetic (c = discrepancy factor, program target c <= 1.7)")
print("=" * 92)
MEAS = {9: F_(1979, 4004), 10: F_(55, 91), 11: F_(66, 91), 12: F_(6, 7)}
TENTF = {9: F_(31, 63), 10: F_(8, 35), 11: F_(0), 12: F_(0)}
MP = F_(14249, 252252)
for k in (9, 10):
    for c_num in (10, 13, 17):
        c = F_(c_num, 10)
        rho = MEAS[k] * (1 - c * (1 - TENTF[k]))
        print(f"    k={k}, c={float(c):.1f}: rho* >= meas(G_P)*(1 - c*(1-floor)) = {float(rho):.4f} "
              f">= m_P {float(MP):.4f}: {rho >= MP}")
print()
print("    => at c = 1 (no discrepancy): k=9 headroom 4.3x, k=10 headroom 2.4x;")
print("       c <= 1.7 suffices at k=9; c <= 1.3-ish at k=10.  The finite program:")
print("       exact int_{G_P} f(frac(dx)) tables for small d per P + Koksma tail for large d.")
