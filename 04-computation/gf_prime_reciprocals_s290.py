#!/usr/bin/env python3
"""
gf_prime_reciprocals_s290.py -- opus-2026-03-24-S290

GENERATING FUNCTIONS INVOLVING PRIME RECIPROCALS

From the OCF: H(T) = I(Ω(T), 2) = finite Euler product over odd cycles.
From the Burnside formula: V_n = (1/n!) Σ_{odd σ} 2^{c(g_E)}

The PRIMES enter through:
1. The cycle index: only ODD cycle lengths contribute
2. The OCF: odd cycles are the "primes" of tournament theory
3. The arc neutrality: numerators 1, 3, 11, 79 involve tournament primes
4. The residual: factors 7×11×19, 79×199 at n=8,9

GENERATING FUNCTION CANDIDATES:

GF_1: The tournament zeta function
  ζ_T(s) = Π_{odd k ≥ 3} (1 - k^{-s})^{-1}
  Evaluated at s values related to n.

GF_2: The prime reciprocal sum
  S(n) = Σ_{p prime, p ≤ 2n-3} 1/p
  Connected to the Mertens function M(n) ≈ ln ln n + M.

GF_3: The Euler product for V_n
  V_n / (2^m / n!) → 1 as n → ∞.
  The correction: V_n × n! / 2^m = 1 + Σ corrections from non-identity σ.
  These corrections involve 2^{c(g_E) - m} for each non-identity odd σ.

Author: opus-2026-03-24-S290
"""

import sys
from math import comb, factorial, gcd, log, log2
from fractions import Fraction
from functools import reduce
from collections import Counter
sys.stdout.reconfigure(line_buffering=True)

def banner(title):
    print("\n" + "="*72)
    print("  " + title)
    print("="*72 + "\n")

def pair_orbits(ct):
    within = sum((c - 1) // 2 for c in ct)
    across = sum(gcd(ct[i], ct[j]) for i in range(len(ct)) for j in range(i+1, len(ct)))
    return within + across

def multinomial(n, ct):
    ct_counter = Counter(ct)
    denom = 1
    for c, a in ct_counter.items():
        denom *= (c ** a) * factorial(a)
    return factorial(n) // denom

def gen_odd_partitions(n):
    result = []
    def _gen(remaining, max_part, current):
        if remaining == 0:
            result.append(tuple(sorted(current)))
            return
        start = min(remaining, max_part)
        if start % 2 == 0: start -= 1
        for p in range(start, 0, -2):
            _gen(remaining - p, p, current + [p])
    _gen(n, n, [])
    return list(set(result))

# ============================================================
banner("THE BURNSIDE CORRECTION AS PRIME PRODUCT")
# ============================================================

# V_n × n! / 2^m = 1 + correction
# correction = Σ_{non-identity odd σ} mult(σ) × 2^{c(g_E) - m} / n!

# For the identity: mult = 1, c(g_E) = m. Contribution = 2^m / n! × n! / 2^m = 1.
# For non-identity: mult(σ) × 2^{c(g_E)} / n!. So correction/1 = this sum.

print("  V_n × n! / 2^m = 1 + correction:\n")

for n in range(3, 16):
    m = comb(n, 2)
    nfact = factorial(n)
    parts = gen_odd_partitions(n)

    total = 0
    identity_term = 0
    for ct in parts:
        mult = multinomial(n, ct)
        po = pair_orbits(ct)
        contrib = mult * pow(2, po)
        total += contrib
        if ct == tuple([1]*n):
            identity_term = contrib

    V = total // nfact
    ratio = Fraction(total, nfact * pow(2, m)) if m < 200 else float(total) / (nfact * pow(2, m))
    correction = ratio - 1 if isinstance(ratio, Fraction) else ratio - 1.0

    if n <= 10:
        print("  n=%2d: V=%d, V×n!/2^m = %s, correction = %s" %
              (n, V, ratio, correction))
    else:
        print("  n=%2d: V=%.4e, correction ≈ %.8f" %
              (n, float(V), float(correction)))

# ============================================================
banner("THE CORRECTION AS A PRODUCT OVER ODD PRIMES")
# ============================================================

# The correction involves terms from 3-cycles, 5-cycles, 7-cycles, etc.
# Each k-cycle contributes: (number of k-cycles in S_n) × 2^{Δ(pair_orbits)}

# For a SINGLE k-cycle with n-k fixed points:
# pair_orbits = C(n-k, 2) + (k-1)/2 + (n-k) × 1
# Δ = pair_orbits - m = C(n-k,2) + (k-1)/2 + (n-k) - C(n,2)
#   = -C(k,2) + (k-1)/2 + (n-k) - k(n-k)
# Hmm, this is getting complex. Let me just compute.

print("  CONTRIBUTION OF EACH CYCLE LENGTH TO V_n:\n")

for n in range(3, 12):
    m = comb(n, 2)
    nfact = factorial(n)
    parts = gen_odd_partitions(n)

    by_max_cycle = {}
    for ct in parts:
        max_c = max(ct)
        mult = multinomial(n, ct)
        po = pair_orbits(ct)
        if max_c not in by_max_cycle:
            by_max_cycle[max_c] = 0
        by_max_cycle[max_c] += mult * pow(2, po)

    if n <= 8:
        print("  n=%d:" % n)
        for k in sorted(by_max_cycle.keys()):
            frac = Fraction(by_max_cycle[k], nfact * pow(2, m))
            print("    max_cycle=%d: contribution = %s = %.8f" % (k, frac, float(frac)))

# ============================================================
banner("PRIME RECIPROCAL GENERATING FUNCTION")
# ============================================================

# Consider: f(n) = log(V_n × n! / 2^m)
# This is the log of (1 + correction).
# If correction is small: f(n) ≈ correction.

# The correction at large n is dominated by the 3-cycle term:
# corr_3 ≈ C(n,3) × 2^{m - m + Δ_3} / n! where Δ_3 = pair_orbits((1,...,1,3)) - m

# At n=5: partition (1,1,3), pair_orbits = 4, m = 10.
# Δ = 4 - 10 = -6. So 2^{-6} = 1/64.
# mult = C(5,3) × 2 / 3 × 1 = 20. (Actually 20 from the multinomial.)
# Contribution: 20 × 2^4 / 120 = 20 × 16 / 120 = 320/120 = 8/3.
# Wait: pair_orbits for (1,1,3) = (3-1)/2 + 2×1 = 1 + 2 = 3... no.
# pair_orbits((1,1,3)): within = (1-1)/2 + (1-1)/2 + (3-1)/2 = 0+0+1 = 1
# across = gcd(1,1) + gcd(1,3) + gcd(1,3) = 1+1+1 = 3. Total = 4.
# At n=5: m=10. So 2^4 = 16 vs 2^10 = 1024. Ratio = 16/1024 = 1/64.
# Contribution to V: 20 × 16 / 120 = 320/120 = 8/3 ≈ 2.67.
# As fraction of identity: 8/3 × 120/(1024) = 320/1024 ≈ 0.3125.

print("""
  THE PRIME RECIPROCAL CONNECTION:

  The correction to V_n from the identity-only approximation is:
  correction(n) = Σ_{k=3,5,7,...} corr_k(n)

  where corr_k(n) comes from partitions containing a k-cycle.

  For a SINGLE k-cycle (one k-cycle + (n-k) fixed points):
  corr_k(n) = C(n, k) × (k-1)! / k × 2^{pair_orbits((1^{n-k}, k)) - m}

  The 2^{...} factor is:
  2^{(k-1)/2 + (n-k) - C(k,2)} = 2^{(k-1)/2 + n - k - k(k-1)/2}
  = 2^{(k-1)/2 × (1-k) + n - k}
  = 2^{n - k - (k-1)(k-1)/2}

  For k=3: 2^{n-3-2} = 2^{n-5}
  For k=5: 2^{n-5-8} = 2^{n-13}
  For k=7: 2^{n-7-18} = 2^{n-25}

  The k-cycle contribution DECAYS as 2^{-k(k-1)/2} relative to the identity.
  This is SUPER-EXPONENTIAL in k (the "prime").

  So: correction(n) ≈ corr_3(n) × (1 + 2^{-6} × ... + 2^{-18} × ... + ...)
  The higher "primes" (k=5,7,...) contribute NEGLIGIBLY.

  THE GENERATING FUNCTION:
  G(x) = Σ_n V_n × n! / 2^{C(n,2)} × x^n
  = Σ_n (1 + correction(n)) × x^n
  ≈ 1/(1-x) + Σ_n correction(n) × x^n

  The correction GF involves:
  Σ_n corr_3(n) x^n = Σ_n C(n,3)(2/3)2^{n-5-m+Δ} x^n = ...

  This is getting complex. Let me just compute the numerical values.
""")

# Compute correction sequence
print("  CORRECTION SEQUENCE (V_n × n! / 2^m - 1):\n")
corrections = []
for n in range(3, 16):
    m = comb(n, 2)
    nfact = factorial(n)
    parts = gen_odd_partitions(n)
    total = sum(multinomial(n, ct) * pow(2, pair_orbits(ct)) for ct in parts)
    V = total // nfact

    if m < 100:
        corr = Fraction(total - nfact * pow(2, m), nfact * pow(2, m))
        corrections.append(float(corr))
        if n <= 10:
            print("  n=%2d: correction = %s = %.10f" % (n, corr, float(corr)))
        else:
            print("  n=%2d: correction ≈ %.12f" % (n, float(corr)))
    else:
        corr = total / (nfact * 2.0**m) - 1
        corrections.append(corr)
        print("  n=%2d: correction ≈ %.12f" % (n, corr))

# Does the correction relate to prime reciprocals?
print("\n  COMPARISON WITH PRIME RECIPROCAL SUMS:")
print("  Σ 1/p for odd primes p ≤ 2n-3:\n")

def is_prime(k):
    if k < 2: return False
    for p in range(2, int(k**0.5)+1):
        if k % p == 0: return False
    return True

for n in range(3, 16):
    # Sum of 1/p for odd primes p ≤ 2n-3
    prime_sum = sum(Fraction(1, p) for p in range(3, 2*n-2, 2) if is_prime(p))
    corr = corrections[n-3]
    if abs(corr) > 0:
        ratio = float(prime_sum) / corr if corr != 0 else 0
    else:
        ratio = 0
    print("  n=%2d: Σ1/p = %s = %.6f, correction = %.8f, ratio = %.4f" %
          (n, prime_sum, float(prime_sum), corr, ratio))

# ============================================================
banner("THE MERTENS CONSTANT CONNECTION")
# ============================================================

# Mertens' theorem: Σ_{p ≤ x} 1/p ≈ ln ln x + M
# where M ≈ 0.2615 is the Meissel-Mertens constant.

# Our correction(n) involves contributions from primes 3, 5, 7, ...
# The leading term from k=3 is proportional to 2^{n-5} × C(n,3) / n!
# = 2^{n-5} / (3! × (n-3)!) × n! / n! ... hmm.

# Let me just look at log(1 + correction):
print("  log(V_n × n! / 2^m) sequence:\n")
for n in range(3, 16):
    corr = corrections[n-3]
    log_corr = log(1 + corr) if corr > -1 else float('nan')
    print("  n=%2d: log(1+corr) = %.10f" % (n, log_corr))

# Ratio of consecutive corrections
print("\n  Correction ratios:")
for i in range(1, len(corrections)):
    if corrections[i-1] != 0:
        ratio = corrections[i] / corrections[i-1]
        print("  n=%d→%d: ratio = %.6f" % (i+3, i+4, ratio))

print("\nDone. opus-2026-03-24-S290")
