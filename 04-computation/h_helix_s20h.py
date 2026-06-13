#!/usr/bin/env python3
"""H as product of distances — helix perspective. kind-pasteur-2026-03-22-S20h"""
import sys
from math import log, sqrt, pi
sys.stdout.reconfigure(line_buffering=True)

print("H AS FREE ENERGY: log H = sum of log-distances from x=2 to each root")
print()

# Energy levels: E_i = -log(2 - r_i) = -log(2 + |r_i|)
# For degree-1 polynomial with root r = -1/alpha_1:
# E = -log(2 + 1/alpha_1) = -log((2*alpha_1 + 1)/alpha_1)

print("ENERGY LEVELS (one root per polynomial, n=5):")
for a in [1, 2, 3, 4, 5, 6, 10, 20, 100]:
    r = -1.0 / a
    E = -log(2 - r)
    H = 2*a + 1
    print(f"  alpha_1={a:>3d}: root={r:>8.4f}, E={E:>8.4f}, H={H:>5d}, logH={log(H):>8.4f}")

print()
print(f"ENERGY BANDWIDTH: [{-log(3):.4f}, {-log(2):.4f}] = width {log(3/2):.6f} = log(3/2)")
print(f"  log(3/2) = {log(3/2):.6f} nats = {log(3/2)/log(2):.6f} bits")
print(f"  3/2 = reciprocal of 2/3 = the modular group generator ratio")
print()

# The RIEMANN SURFACE perspective
# Root at r = -1/a has log(r) = -log(a) + i*pi (on the principal branch)
# Point x=2 has log(2) = 0.693
# Riemann distance = |log(2) - log(r)| = sqrt((log(2)+log(a))^2 + pi^2)

print("RIEMANN SURFACE DISTANCES (from log(2) to log(root)):")
print(f"  {'alpha':>6s} {'root':>8s} {'|2-r|':>8s} {'Riem':>8s} {'amp':>8s}")
for a in [1, 2, 4, 5, 6, 10, 50, 100]:
    r = -1.0/a
    std = 2 + 1.0/a
    riem = sqrt((log(2) + log(a))**2 + pi**2)
    amp = (a * riem) / (2*a + 1) if a > 0 else 0
    print(f"  {a:>6d} {r:>8.4f} {std:>8.4f} {riem:>8.4f} {amp:>8.4f}")

print()
print(f"AMPLIFICATION LIMIT: a*Riem / (2a+1) -> pi/2 = {pi/2:.6f} as a->inf")
print()

# Check the limit
# Riem ~ sqrt(log(2a)^2 + pi^2) ~ sqrt(log(a)^2 + pi^2)
# a * Riem / (2a+1) ~ a * sqrt(log(a)^2 + pi^2) / (2a) = sqrt(log(a)^2 + pi^2)/2
# This -> infinity, not pi/2. Let me recheck.
# Actually H_Riem = a * Riem where Riem = distance in LOG space.
# But H_std = |2-r| * a = (2+1/a)*a = 2a+1.
# Riem = sqrt((log(2)+log(a))^2 + pi^2) = sqrt(log(2a)^2 + pi^2).
# a * Riem = a * sqrt(log(2a)^2 + pi^2). This grows as a*log(a).
# H_std = 2a+1 ~ 2a. Ratio ~ sqrt(log(2a)^2 + pi^2)/2 -> infinity.
# So the amplification grows WITHOUT BOUND. My earlier claim of pi/2 was wrong.

print("CORRECTION: amplification grows as sqrt(log(a)^2 + pi^2)/2, not pi/2.")
print("The helix amplification is UNBOUNDED as alpha_1 -> inf.")
print()

# But the pi^2 term in the Riemann distance IS interesting.
# For small alpha_1 (few cycles), log(a)^2 is small and pi^2 dominates:
# Riem ~ pi (the distance across the branch cut).
# For large alpha_1, log(a)^2 dominates and pi^2 becomes negligible.

print("TWO REGIMES:")
print("  Small alpha_1: Riemann distance ~ pi (dominated by the branch cut)")
print("  Large alpha_1: Riemann distance ~ log(alpha_1) (dominated by magnitude)")
print(f"  Crossover: when log(alpha_1) ~ pi, i.e., alpha_1 ~ e^pi = {pow(2.71828, pi):.1f}")
print(f"  e^pi = {pow(2.71828, pi):.4f} ~ 23. The crossover is at alpha_1 ~ 23.")
print(f"  23 = |BT| - 1 = the moonshine number!")
print()

# Actually e^pi = 23.14... and 23 IS a significant number:
# 23 is the number of umbral moonshine cases.
# |BT| - 1 = 24 - 1 = 23.
# At alpha_1 ~ 23: the branch-cut distance (pi) equals the magnitude distance (log 23).
# Below 23 cycles: the topology of the complex plane (the branch cut at pi) dominates.
# Above 23 cycles: the size of the roots (log alpha_1) dominates.

print("THE MOONSHINE CROSSOVER:")
print(f"  At alpha_1 = 23: log(23) = {log(23):.4f}, pi = {pi:.4f}")
print(f"  log(23)/pi = {log(23)/pi:.4f} ~ 1. The crossover point.")
print(f"  23 = |BT| - 1 = umbral moonshine count.")
print()
print("  Below 23 cycles: the TOPOLOGY of the complex plane (branch cut)")
print("  controls the Riemann distance. This is the modular regime.")
print("  Above 23 cycles: the MAGNITUDE of the roots controls.")
print("  This is the classical (non-modular) regime.")
print()

# A more practical observation:
# log H = log(alpha_d) + sum log(2 - r_i)
# For all-real-negative roots: log H = log(alpha_d) + sum log(2 + |r_i|)
# Each term is positive. log H is a SUM of positive contributions.
# This means: log H is a CONCAVE function of the root positions.
# Moving roots closer together INCREASES log H (by Jensen's inequality
# on the convex function -log).

# Actually: is log H concave in the root positions?
# d(log H)/dr_i = -1/(2 - r_i) < 0 (since 2-r_i > 0).
# d^2(log H)/dr_i^2 = -1/(2-r_i)^2 < 0. YES, concave.
# So log H is concave in each root position.
# MOVING A ROOT CLOSER TO x=2 (less negative) INCREASES H.
# The tournament that MAXIMIZES H has roots as close to 2 as possible.

print("CONCAVITY OF log H:")
print("  d^2(log H)/dr_i^2 = -1/(2-r_i)^2 < 0 for all r_i < 2.")
print("  log H is CONCAVE in each root position.")
print("  CONSEQUENCE: among all polynomials with the same degree and")
print("  leading coefficient, H is maximized when roots are as close")
print("  to x=2 as possible (least negative).")
print()
print("  For tournament Omega: roots r_i = -1/m_i where m_i are related")
print("  to the sizes of maximal cliques in Omega.")
print("  Smaller cliques (smaller m_i) -> roots closer to 0 -> farther from 2")
print("  -> LOWER H. Larger cliques -> roots closer to 0... hmm.")
print()

# Actually: for degree-1, root = -1/alpha_1.
# Larger alpha_1 -> root closer to 0 -> farther from -infinity, closer to 0.
# |2 - r| = 2 + 1/alpha_1, which DECREASES as alpha_1 increases.
# So H = alpha_1 * (2 + 1/alpha_1) = 2*alpha_1 + 1.
# H INCREASES with alpha_1 even though |2-r| decreases.
# Because the leading coefficient alpha_1 grows faster than the distance shrinks.

# For degree d: H = alpha_d * prod(2 - r_i).
# If we ADD a root (increase degree by 1), we multiply H by (2 - r_new) * (alpha_{d+1}/alpha_d).
# Adding a root with r_new < 0: factor (2 - r_new) > 2. So H roughly doubles per added root.
# This is the EXPONENTIAL GROWTH of H with the number of independent cycle packings.

print("EXPONENTIAL GROWTH:")
print("  Each new root (new alpha_k coefficient) multiplies H by > 2.")
print("  d roots -> H > alpha_d * 2^d.")
print("  For floor(n/3) roots: H > alpha_d * 2^{floor(n/3)}.")
print("  This gives a LOWER BOUND on H from the polynomial degree alone.")
print()

# Compute the lower bound
for n_val in [3, 4, 5, 6, 7, 8, 9, 10]:
    d = n_val // 3
    lower = 2**d  # alpha_d >= 1, each distance >= 2
    print(f"  n={n_val}: degree={d}, lower bound 2^{d} = {2**d}")
