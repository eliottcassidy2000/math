#!/usr/bin/env python3
"""
odd_even_ratio_limit.py — opus-2026-03-13-S67i

Investigate the limit log(F_odd)/log(F_even) → 6 discovered in
fibonacci_word_ergodic.py. This is an UNEXPECTED integer limit.

Also reconcile the Paley identification issue: at p=11, det(I+A)
gives different values than expected from THM-162. Need to check
whether the circulant parameterization matches the QR-based Paley.

ANALYTICAL DERIVATION of the ratio:
  log(F_odd) = sum_{k odd} log(1 + Q_k)
  log(F_even) = sum_{k even} log(1 + Q_k)

  For large p:
  - Odd modes: Q_k = 1/(4sin^2(k*pi/(2p))) ~ p^2/(pi^2*k^2) for small k
    So log(1+Q_k) ~ log(Q_k) ~ 2*log(p) - 2*log(k) - log(pi^2)
    But there are ~m/2 odd modes, each contributing differently

  - Even modes: Q_k = 1/(4cos^2(k*pi/(2p))) -> 1/4 for all k
    So log(1+Q_k) -> log(5/4) for each of ~m/2 even modes
    Total: log(F_even) ~ (m/2)*log(5/4)

  - For odd modes: the dominant contribution comes from the integral
    sum_{k odd, 1..m} log(1+Q_k) ~ integral of log(1 + 1/(4sin^2(u))) du
    from 0 to pi/2
    = integral of log((4sin^2(u) + 1)/(4sin^2(u))) du
    = integral of log(4sin^2(u) + 1) du - integral of log(4sin^2(u)) du

    The second integral: int_0^{pi/2} log(4sin^2(u)) du
    = pi*log(2)/2 + 2*int_0^{pi/2} log(sin(u)) du = pi*log(2)/2 - pi*log(2) = -pi*log(2)/2
    Wait, int_0^{pi/2} log(sin u) du = -pi*log(2)/2
    So int_0^{pi/2} log(4sin^2 u) du = pi*log(4)/2 + 2*(-pi*log(2)/2) = pi*log(4)/2 - pi*log(2) = 0

    And int_0^{pi/2} log(4sin^2 u + 1) du = ?
    Let's compute this numerically and see if the ratio works out to 6.
"""

import math
import numpy as np

print("=" * 70)
print("ODD/EVEN RATIO LIMIT ANALYSIS")
print("=" * 70)

phi = (1 + math.sqrt(5)) / 2

# First: verify the ratio numerically at large p
print("\nlog(F_odd)/log(F_even) for large p:")
primes_3mod4 = []
for p in range(7, 5000, 2):
    if p % 4 != 3:
        continue
    if all(p % d != 0 for d in range(3, int(p**0.5)+1, 2)):
        primes_3mod4.append(p)

for p in primes_3mod4[:40]:
    m = (p - 1) // 2
    log_odd = 0
    log_even = 0
    for k in range(1, m + 1):
        if k % 2 == 1:
            Q = 1.0 / (4 * math.sin(k * math.pi / (2 * p))**2)
            log_odd += math.log(1 + Q)
        else:
            Q = 1.0 / (4 * math.cos(k * math.pi / (2 * p))**2)
            log_even += math.log(1 + Q)
    ratio = log_odd / log_even if log_even > 0 else 0
    if p <= 200 or p > 4000:
        print(f"  p={p:5d}: ratio = {ratio:.8f}")

# Analytical computation of the ratio in the limit
print("\n\n=== ANALYTICAL DERIVATION ===")

# The key integrals (over [0, pi/2], normalized by m ~ p/2):
# log(F_even)/m ~ (1/2) * log(5/4) = log(5/4)/2

# For odd modes: log(F_odd) ~ (m/pi) * int_0^{pi/2} log(1 + 1/(4sin^2 u)) du
# where we substituted k*pi/(2p) -> u, dk -> (2p/pi) du, and divided by 2 (odd k only)
# Actually: sum over odd k from 1 to m of f(k*pi/(2p))
# ~ (p/pi) * int_0^{pi/2} f(u) du [sum over all k]
# But only odd k contributes, so factor of 1/2

# Wait, let me be more careful.
# For PALEY Fibonacci product, Q_k for k=1..m:
# k odd: Q_k = 1/(4sin^2(k*pi/(2p)))
# k even: Q_k = 1/(4cos^2(k*pi/(2p)))

# Contribution of odd modes:
# S_odd = sum_{k=1,3,5,...,m} log(1 + 1/(4sin^2(k*pi/(2p))))
# The odd k from 1 to m: there are roughly m/2 terms
# With u_k = k*pi/(2p), spacing Delta u = pi/p (since k increments by 2)
# S_odd ~ (p/pi) * int_0^{pi/2} log(1 + 1/(4sin^2 u)) du / 2

# Wait: k goes from 1 to m by 2. So k = 1,3,...,m (if m odd) or m-1 (if m even)
# u_k = k*pi/(2p), step in u is 2*pi/(2p) = pi/p
# S_odd ~ (1/(pi/p)) * int = (p/pi) * int_0^{pi/2} log(1 + 1/(4sin^2 u)) du
# But this counts both odd and even terms! For odd only, multiply by 1/2:
# Actually no. If we sum over k=1,3,5,...  with step 2, and u_k = k*pi/(2p):
# du = (step in k) * pi/(2p) = 2 * pi/(2p) = pi/p
# Number of terms ~ m/2
# S_odd ~ (m/2) * (2/m) * (p/pi) * int * ...
#
# Let me just compute numerically:

# scipy not available, using manual quadrature

# Integral of log(1 + 1/(4sin^2 u)) from 0 to pi/2
# Note: this diverges logarithmically at u=0!

def integrand_odd(u):
    if u < 1e-15:
        return 0  # regularize
    return math.log(1 + 1/(4*math.sin(u)**2))

def integrand_even(u):
    return math.log(1 + 1/(4*math.cos(u)**2))

# Numerical integration
N_quad = 10000
du = (math.pi / 2) / N_quad
I_odd = sum(integrand_odd(du/2 + i*du) * du for i in range(N_quad))
I_even = sum(integrand_even(du/2 + i*du) * du for i in range(N_quad))

print(f"\nNumerical integrals over [0, pi/2]:")
print(f"  I_odd  = int log(1 + 1/(4sin^2 u)) du = {I_odd:.8f}")
print(f"  I_even = int log(1 + 1/(4cos^2 u)) du = {I_even:.8f}")
print(f"  I_odd / I_even = {I_odd/I_even:.8f}")

# But wait — the odd integral diverges! It's log(Q) ~ log(1/u^2) near u=0
# The actual sum S_odd has a finite number of terms starting at u_1 = pi/(2p)
# So S_odd ~ (p/(2pi)) * [int + correction from cutoff]

# Let me compute the per-mode averages more carefully
print("\nPer-mode average contributions:")
for p in [997, 1999, 4999]:
    m = (p - 1) // 2
    n_odd = (m + 1) // 2  # number of odd k in 1..m
    n_even = m // 2        # number of even k in 1..m

    log_odd = 0
    log_even = 0
    for k in range(1, m + 1):
        if k % 2 == 1:
            Q = 1.0 / (4 * math.sin(k * math.pi / (2 * p))**2)
            log_odd += math.log(1 + Q)
        else:
            Q = 1.0 / (4 * math.cos(k * math.pi / (2 * p))**2)
            log_even += math.log(1 + Q)

    avg_odd = log_odd / n_odd
    avg_even = log_even / n_even

    print(f"  p={p}: n_odd={n_odd}, n_even={n_even}")
    print(f"    avg log(1+Q) per odd mode:  {avg_odd:.6f}")
    print(f"    avg log(1+Q) per even mode: {avg_even:.6f} (log(5/4)={math.log(5/4):.6f})")
    print(f"    ratio of averages: {avg_odd/avg_even:.6f}")
    print(f"    ratio of totals:   {log_odd/log_even:.6f}")
    print(f"    n_odd/n_even = {n_odd/n_even:.6f}")
    print(f"    total ratio = (avg_odd/avg_even) * (n_odd/n_even) = "
          f"{(avg_odd/avg_even) * (n_odd/n_even):.6f}")

# The total ratio = (n_odd/n_even) * (avg_odd/avg_even)
# n_odd/n_even -> 1 as m -> infinity
# So the ratio -> avg_odd/avg_even
# avg_even -> log(5/4)
# avg_odd -> diverges (because of the k=1 term which grows as log(p))

# AHA! The ratio is NOT converging to 6!
# Let me check larger p more carefully

print("\n\nRatio at very large p:")
for p in primes_3mod4:
    if p < 500:
        continue
    m = (p - 1) // 2
    log_odd = 0
    log_even = 0
    for k in range(1, m + 1):
        if k % 2 == 1:
            Q = 1.0 / (4 * math.sin(k * math.pi / (2 * p))**2)
            log_odd += math.log(1 + Q)
        else:
            Q = 1.0 / (4 * math.cos(k * math.pi / (2 * p))**2)
            log_even += math.log(1 + Q)
    ratio = log_odd / log_even
    print(f"  p={p:5d}: ratio = {ratio:.8f}, "
          f"log_odd={log_odd:.4f}, log_even={log_even:.4f}")
    if p > 3000:
        break

print("\n\n=== THE RATIO DIVERGES! ===")
print("The ratio log(F_odd)/log(F_even) is NOT converging to 6.")
print("It's slowly increasing because log(F_odd) grows as m*log(phi)*2")
print("while log(F_even) grows as m*log(5/4)/2.")
print(f"The asymptotic ratio should be: 2*log(phi) / (log(5/4)/2)")
print(f"  = 4*log(phi)/log(5/4) = {4*math.log(phi)/math.log(5/4):.6f}")

# Wait, that's not right either. Let me think...
# log(F_p) = sum_all log(1+Q_k) ~ m * 2*log(phi) [from kappa_analytical]
# Actually: log(F_p) ~ 2m*log(phi)
# log(F_even) ~ m/2 * log(5/4) [each even mode contributes ~log(5/4)]
# log(F_odd) = log(F_p) - log(F_even) ~ 2m*log(phi) - m*log(5/4)/2
# ratio = (2m*log(phi) - m*log(5/4)/2) / (m*log(5/4)/2)
#        = 4*log(phi)/log(5/4) - 1

asymptotic_ratio = 4 * math.log(phi) / math.log(5/4) - 1
print(f"\n  Corrected asymptotic ratio = 4*log(phi)/log(5/4) - 1 = {asymptotic_ratio:.6f}")
print(f"  Compare with observed: ~6.0")

# Hmm, 4*log(phi)/log(5/4) - 1
val = 4 * math.log(phi) / math.log(5/4)
print(f"\n  4*log(phi)/log(5/4) = {val:.6f}")
print(f"  4*log(phi)/log(5/4) - 1 = {val - 1:.6f}")

# So the ratio -> 4*log(phi)/log(5/4) - 1 ≈ ?
# log(phi) = 0.48121, log(5/4) = 0.22314
# 4*0.48121/0.22314 = 8.627..., minus 1 = 7.627
# That's not 6 either!

# The issue is that log(F_p)/m is NOT exactly 2*log(phi).
# From kappa_analytical: log(F_p)/m = 2*log(phi) asymptotically,
# but the convergence is slow.

# Let me compute the ACTUAL log(F)/m and log(F_even)/m
print("\n\nActual growth rates:")
for p in primes_3mod4:
    if p < 100:
        continue
    m = (p - 1) // 2
    log_odd = 0
    log_even = 0
    for k in range(1, m + 1):
        if k % 2 == 1:
            Q = 1.0 / (4 * math.sin(k * math.pi / (2 * p))**2)
            log_odd += math.log(1 + Q)
        else:
            Q = 1.0 / (4 * math.cos(k * math.pi / (2 * p))**2)
            log_even += math.log(1 + Q)

    log_total = log_odd + log_even
    n_even = m // 2

    print(f"  p={p:5d}, m={m:4d}: log(F)/m={log_total/m:.6f}, "
          f"log(F_odd)/m={log_odd/m:.6f}, "
          f"log(F_even)/m={log_even/m:.6f}, "
          f"log(F_even)/(m/2)={log_even/n_even:.6f}")
    if p > 3000:
        break

print(f"\n  Expected log(F)/m -> 2*log(phi) = {2*math.log(phi):.6f}")
print(f"  Expected log(F_even)/(m/2) -> log(5/4) = {math.log(5/4):.6f}")

# So log(F_even)/m -> log(5/4)/2 = 0.11157
# log(F_odd)/m -> 2*log(phi) - log(5/4)/2 = 0.96242 - 0.11157 = 0.85085
# Ratio = log(F_odd)/log(F_even) = (0.85085*m) / (0.11157*m) = 7.627
#
# But we observed ~6.0 at p~200. So the convergence to ~7.6 is VERY slow
# because of the log(p) correction in log(F_odd).
#
# Actually wait — the log(p) term in the k=1 mode:
# log(1+Q_1) ~ log(p^2/(pi^2)) ~ 2*log(p) for large p
# This is O(log p), while the rest of the sum is O(m) = O(p).
# So the correction is O(log(p)/p) which is slow but vanishes.
# The ratio at finite p is BELOW the asymptotic value.

print(f"\n  Asymptotic ratio = (2*log(phi) - log(5/4)/2) / (log(5/4)/2)")
print(f"                   = {(2*math.log(phi) - math.log(5/4)/2) / (math.log(5/4)/2):.6f}")
print(f"                   = 4*log(phi)/log(5/4) - 1")
print(f"                   = {4*math.log(phi)/math.log(5/4) - 1:.6f}")

# Now: is there a clean form?
# 4*log(phi)/log(5/4) = 4*log(phi)/log(5/4)
# log(phi) = log((1+sqrt(5))/2)
# log(5/4) = log(5) - log(4) = log(5) - 2*log(2)
# Not obvious. Let's check:
print(f"\n  Checking for clean form:")
print(f"  4*log(phi)/log(5/4) = {4*math.log(phi)/math.log(5/4):.10f}")
print(f"  Is it close to an integer? {abs(val - round(val)):.6f}")
print(f"  Nearest integer: {round(val)}")
print(f"  Ratio/phi = {val/phi:.6f}")
print(f"  Ratio - phi = {val - phi:.6f}")
print(f"  Ratio*log(2) = {val*math.log(2):.6f}")

# The actual ratio 4*log(phi)/log(5/4) ≈ 8.627...
# The observed slow convergence from below ~6 is expected

print("\n\n=== RECONCILE WITH PALEY det(I+A) ===")
print("The circulant parameterization needs checking.")

# At p=7, Paley tournament has QR set = {1, 2, 4} mod 7
# The quadratic residues mod 7: 1^2=1, 2^2=4, 3^2=2 -> QR = {1,2,4}
# In our parameterization, S = set of k in {1,..,m} such that k is a QR
# For p=7, m=3: QR ∩ {1,2,3} = {1,2} (since 3 is NQR mod 7)
# So Paley should be S={1,2}, NOT S={1,2,3}!

print("\nPaley orientation set (QRs in {1,..,m}):")
for p in [7, 11, 19, 23, 31, 43]:
    m = (p - 1) // 2
    QR = set()
    for k in range(1, p):
        QR.add((k*k) % p)
    # QR ∩ {1,..,m}
    paley_S = QR & set(range(1, m+1))
    print(f"  p={p}: QR mod p = {sorted(QR)}, Paley S = {sorted(paley_S)}")

# Now recompute det(I+A) with CORRECT Paley identification
print("\ndet(I+A) with CORRECT Paley identification:")
for p in [7, 11]:
    m = (p - 1) // 2
    QR = set()
    for k in range(1, p):
        QR.add((k*k) % p)
    paley_S = QR & set(range(1, m+1))

    print(f"\n  p={p}, Paley S = {sorted(paley_S)}:")

    results = []
    for bits in range(1 << m):
        S_set = set()
        for k in range(m):
            if bits & (1 << k):
                S_set.add(k + 1)

        A = np.zeros((p, p), dtype=int)
        for i in range(p):
            for k in range(1, m+1):
                if k in S_set:
                    j = (i + k) % p
                else:
                    j = (i - k) % p
                A[i][j] = 1

        eigenvals = np.linalg.eigvals(A.astype(float))
        det_val = abs(np.prod(1 + eigenvals))
        sum4 = sum(abs(e)**4 for e in eigenvals)

        # Spectral variance (non-zero eigenvalues only)
        mags = sorted([abs(e) for e in eigenvals if abs(e) > 0.01])
        var_mag = np.var(mags) if len(mags) > 0 else 0

        is_paley = (S_set == paley_S)
        results.append((S_set, det_val, sum4, var_mag, is_paley))

    results.sort(key=lambda x: -x[1])
    for S_set, det_val, sum4, var_mag, is_paley in results[:5]:
        tag = " <-- PALEY" if is_paley else ""
        print(f"    S={str(sorted(S_set)):20s}: det={det_val:10.2f}, "
              f"sum4={sum4:8.1f}, var={var_mag:.4f}{tag}")
    print(f"    ...")
    for S_set, det_val, sum4, var_mag, is_paley in results[-3:]:
        tag = " <-- PALEY" if is_paley else ""
        print(f"    S={str(sorted(S_set)):20s}: det={det_val:10.2f}, "
              f"sum4={sum4:8.1f}, var={var_mag:.4f}{tag}")

    # Find Paley
    paley_result = [r for r in results if r[4]]
    if paley_result:
        pr = paley_result[0]
        rank_det = sorted(results, key=lambda x: -x[1]).index(pr) + 1
        rank_sum4 = sorted(results, key=lambda x: x[2]).index(pr) + 1
        rank_var = sorted(results, key=lambda x: x[3]).index(pr) + 1
        print(f"\n    Paley rank by det(I+A):    {rank_det}/{len(results)}")
        print(f"    Paley rank by sum|lam|^4:  {rank_sum4}/{len(results)} (1=min)")
        print(f"    Paley rank by spec_var:    {rank_var}/{len(results)} (1=min)")

print("\n\nDONE — odd_even_ratio_limit.py complete")
