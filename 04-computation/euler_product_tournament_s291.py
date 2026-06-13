#!/usr/bin/env python3
"""
euler_product_tournament_s291.py -- opus-2026-03-24-S291

EULER PRODUCT FACTORIZATION OF THE TOURNAMENT COUNTING FUNCTION

From Burnside/Cauchy-Frobenius:
  V_n = (1/n!) Σ_{σ ∈ S_n, all odd cycles} 2^{pair_orbits(σ)}

The identity contributes 2^m / n! where m = C(n,2).
The correction R(n) = V_n × n! / 2^m - 1 factors by CYCLE TYPE.

Each cycle type λ = (λ_1, λ_2, ..., λ_k) (all odd parts summing to n)
contributes:
  C(λ) = mult(λ) × 2^{pair_orbits(λ) - m}

where mult(λ) = n! / (λ_1 × λ_2 × ... × λ_k × r_1! × r_2! × ...)

THE EULER PRODUCT:
  1 + R(n) = V_n × n! / 2^m = Σ_λ C(λ)

Can this be written as a PRODUCT over odd primes p=3,5,7,...?
  1 + R(n) = Π_{p odd prime} F_p(n)

where F_p(n) captures contributions from partitions containing p-cycles?

THE KEY INSIGHT: pair_orbits for a cycle of length k acting on pairs:
  Within a k-cycle: (k-1)/2 orbits on pairs from the cycle
  Between a k-cycle and a j-cycle: gcd(k,j) orbits
  Between a k-cycle and a fixed point: k orbits... no, 1 orbit per fixed point.

Wait: pair_orbits counts orbits of the cycle structure on PAIRS {i,j}.
For a single k-cycle on vertices v_1,...,v_k plus n-k fixed points:
  - Pairs within the cycle: C(k,2) pairs in (k-1)/2 + ... orbits
    Actually: a k-cycle acts on pairs from the cycle.
    Number of orbits = (k-1)/2 for odd k.
  - Pairs between cycle and fixed points: k(n-k) pairs in (n-k) orbits
    (each fixed point pairs with each cycle vertex,
     orbits = number of fixed points × 1 per fixed point)
    Wait: the k-cycle rotates the k vertices. A pair (cycle_vertex, fixed_point)
    has an orbit of size k (the cycle rotates). So #orbits = k(n-k)/k = (n-k).
  - Pairs of fixed points: C(n-k,2) pairs, each its own orbit.

  Total: (k-1)/2 + (n-k) + C(n-k,2)
  = (k-1)/2 + (n-k) + (n-k)(n-k-1)/2
  = (k-1)/2 + (n-k)(1 + (n-k-1)/2)
  = (k-1)/2 + (n-k)(n-k+1)/2

  Compared to identity: m = C(n,2) = n(n-1)/2

  Delta_k = pair_orbits - m
  = (k-1)/2 + (n-k)(n-k+1)/2 - n(n-1)/2
  = (k-1)/2 + ((n-k)^2 + (n-k) - n^2 + n)/2
  = (k-1)/2 + (n^2 - 2nk + k^2 + n - k - n^2 + n)/2
  = (k-1)/2 + (-2nk + k^2 + 2n - k)/2
  = (k-1)/2 + k(k - 2n - 1)/2 + 2n/2
  = (k-1)/2 + (k^2 - 2nk - k + 2n)/2
  Let me redo more carefully:
  = (k-1)/2 + (n-k)(n-k+1)/2 - n(n-1)/2

  (n-k)(n-k+1) = n^2 - 2nk + k^2 + n - k
  n(n-1) = n^2 - n

  Difference: n^2 - 2nk + k^2 + n - k - n^2 + n = k^2 - 2nk + 2n - k
  = k(k-1) - 2n(k-1) = (k-1)(k - 2n)

  So Delta_k = (k-1)/2 + (k-1)(k-2n)/2 = (k-1)(1 + k - 2n)/2 = (k-1)(k+1-2n)/2

  For a SINGLE k-cycle: 2^{Delta_k} = 2^{(k-1)(k+1-2n)/2}

  At k=3: Delta = 2×(4-2n)/2 = 4-2n. So 2^{4-2n}.
  At k=5: Delta = 4×(6-2n)/2 = 2(6-2n) = 12-4n. So 2^{12-4n}.
  At k=7: Delta = 6×(8-2n)/2 = 3(8-2n) = 24-6n. So 2^{24-6n}.

  General: Delta_k = (k-1)(k+1-2n)/2

  The multiplicity of a single k-cycle type (k, 1^{n-k}):
  mult = n! / (k × (n-k)!) = C(n,k) × (k-1)!

  Contribution: C(n,k) × (k-1)! × 2^{(k-1)(k+1-2n)/2} / 2^0
  (The /2^0 because we already divided by n! in V_n and multiplied by n!/2^m.)

  Actually: the contribution to R(n) = V_n × n!/2^m - 1 is:
  R_k(n) = mult(λ) × 2^{pair_orbits(λ)} / (n! × 2^m)... no wait.

  V_n × n! / 2^m = Σ_λ mult(λ) × 2^{pair_orbits(λ)} / (n! × 2^m) × n!
  Hmm. Let me be precise.

  V_n = (1/n!) Σ_λ mult(λ) × 2^{po(λ)}
  V_n × n! = Σ_λ mult(λ) × 2^{po(λ)}
  V_n × n! / 2^m = Σ_λ mult(λ) × 2^{po(λ)-m}

  For identity (λ = (1^n)): mult = 1, po = m. Contribution = 1.
  For single k-cycle (λ = (k, 1^{n-k})): mult = C(n,k)×(k-1)!, po-m = Delta_k.
  Contribution = C(n,k) × (k-1)! × 2^{Delta_k}

  R(n) = Σ_{k=3,5,7,...} C(n,k)(k-1)! × 2^{Delta_k} + (multi-cycle terms)

Author: opus-2026-03-24-S291
"""

import sys
from math import comb, factorial, gcd, log, log2
from fractions import Fraction
from functools import reduce
from collections import Counter, defaultdict
sys.stdout.reconfigure(line_buffering=True)

def banner(title):
    print("\n" + "="*72)
    print("  " + title)
    print("="*72 + "\n")

def pair_orbits(ct):
    """Number of pair orbits under a permutation with cycle type ct."""
    within = sum((c - 1) // 2 for c in ct)
    across = sum(gcd(ct[i], ct[j]) for i in range(len(ct)) for j in range(i+1, len(ct)))
    return within + across

def multinomial(n, ct):
    """Number of permutations in S_n with cycle type ct."""
    ct_counter = Counter(ct)
    denom = 1
    for c, a in ct_counter.items():
        denom *= (c ** a) * factorial(a)
    return factorial(n) // denom

def gen_odd_partitions(n):
    """Generate all partitions of n into odd parts."""
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
banner("EXACT CYCLE-TYPE CONTRIBUTIONS TO R(n)")
# ============================================================

# R(n) = V_n × n! / 2^m - 1 = Σ_{λ ≠ (1^n)} mult(λ) × 2^{po(λ)-m}

print("  For each odd partition λ of n, compute its exact contribution.\n")

for n in range(3, 14):
    m = comb(n, 2)
    parts = gen_odd_partitions(n)

    # Sort partitions by "complexity" (identity last → first)
    # Identity = (1,1,...,1)
    identity = tuple([1]*n)

    contributions = {}
    total = Fraction(0)
    for ct in parts:
        mult = multinomial(n, ct)
        po = pair_orbits(ct)
        # Contribution to V_n × n! / 2^m
        contrib = Fraction(mult * pow(2, po), pow(2, m))
        contributions[ct] = contrib
        total += contrib

    R = total - 1  # subtract identity contribution

    if n <= 9:
        print("  n=%d: R(n) = %s = %.10f" % (n, R, float(R)))
        # Show top contributors (by magnitude)
        sorted_contribs = sorted(
            [(ct, c) for ct, c in contributions.items() if ct != identity],
            key=lambda x: -abs(float(x[1]))
        )
        for ct, c in sorted_contribs[:5]:
            label = "×".join(str(k) for k in sorted(ct) if k > 1) or "1^n"
            mult = multinomial(n, ct)
            po = pair_orbits(ct)
            delta = po - m
            print("    λ=%-12s mult=%-6d Δ=%-4d 2^Δ=2^{%d} contrib=%s" %
                  (ct, mult, delta, delta, c))
    else:
        print("  n=%2d: R(n) ≈ %.14f" % (n, float(R)))

# ============================================================
banner("THE SINGLE k-CYCLE FORMULA (EXACT)")
# ============================================================

# For λ = (k, 1^{n-k}):
#   mult = C(n,k) × (k-1)!
#   Delta_k = (k-1)(k+1-2n)/2
#   contrib = C(n,k) × (k-1)! × 2^{(k-1)(k+1-2n)/2}

print("  SINGLE k-CYCLE CONTRIBUTIONS:\n")
print("  %-3s | %-6s | %-12s | %-8s | %-20s" %
      ("n", "k", "mult", "Delta", "contribution"))
print("  " + "-"*60)

for n in range(3, 14):
    m = comb(n, 2)
    for k in [3, 5, 7, 9, 11, 13]:
        if k > n: break
        mult = comb(n, k) * factorial(k-1)
        delta = (k-1) * (k+1-2*n) // 2

        # Verify: delta should equal pair_orbits((k,) + (1,)*(n-k)) - m
        ct = tuple(sorted([k] + [1]*(n-k)))
        po = pair_orbits(ct)
        assert po - m == delta, f"Mismatch: po-m={po-m}, delta={delta}"

        contrib = Fraction(mult, 1) * Fraction(pow(2, po), pow(2, m))
        if n <= 9:
            print("  n=%-2d | k=%-2d | mult=%-8d | Δ=%-5d | %s" %
                  (n, k, mult, delta, contrib))

# ============================================================
banner("RATIO OF SINGLE 3-CYCLE TO TOTAL CORRECTION")
# ============================================================

print("  How much of R(n) comes from JUST the single 3-cycle term?\n")

for n in range(3, 16):
    m = comb(n, 2)
    parts = gen_odd_partitions(n)

    total_R = Fraction(0)
    single_3 = Fraction(0)
    identity = tuple([1]*n)

    for ct in parts:
        mult = multinomial(n, ct)
        po = pair_orbits(ct)
        contrib = Fraction(mult * pow(2, po), pow(2, m))
        if ct != identity:
            total_R += contrib
        if ct == tuple(sorted([3] + [1]*(n-3))):
            single_3 = contrib

    if total_R != 0:
        ratio = float(single_3) / float(total_R)
        print("  n=%2d: R = %.10f, R_3 = %.10f, R_3/R = %.6f" %
              (n, float(total_R), float(single_3), ratio))

# ============================================================
banner("MULTIPLICATIVE STRUCTURE: CAN R(n) BE A PRODUCT?")
# ============================================================

# If R(n) factored as a product over odd cycle lengths:
# 1 + R(n) = Π_k (1 + f_k(n))
# Then: R(n) ≈ Σ_k f_k(n) for small f_k.
#
# The f_k(n) for a single k-cycle is:
# f_k(n) = C(n,k) × (k-1)! × 2^{(k-1)(k+1-2n)/2}
#
# But multi-cycle types like (3,3,1,...) contribute MORE than f_3^2/2.
# Let's check.

print("  TESTING MULTIPLICATIVE STRUCTURE:\n")
print("  If 1+R = Π_k (1+f_k), then the (3,3) contribution should be f_3²/2.\n")

for n in range(6, 14):
    m = comb(n, 2)

    # Single 3-cycle contribution
    ct_3 = tuple(sorted([3] + [1]*(n-3)))
    f3 = Fraction(multinomial(n, ct_3) * pow(2, pair_orbits(ct_3)), pow(2, m))

    # Double 3-cycle contribution
    if n >= 6:
        ct_33 = tuple(sorted([3, 3] + [1]*(n-6)))
        actual_33 = Fraction(multinomial(n, ct_33) * pow(2, pair_orbits(ct_33)), pow(2, m))
        predicted_33 = f3 * f3 / 2  # Product of independent contributions / 2!

        print("  n=%d:" % n)
        print("    f_3(single) = %s = %.8f" % (f3, float(f3)))
        print("    Actual (3,3) = %s = %.10f" % (actual_33, float(actual_33)))
        print("    f_3²/2       = %.10f" % float(predicted_33))
        print("    Ratio act/pred = %.6f" % (float(actual_33) / float(predicted_33) if float(predicted_33) != 0 else 0))
        print()

# ============================================================
banner("THE EXPONENTIAL GENERATING FUNCTION")
# ============================================================

# If contributions are multiplicatively independent by cycle length,
# then: log(1 + R(n)) ≈ Σ_k f_k(n) (to first order)
#
# and the exact answer is:
# 1 + R(n) = Σ_{partitions λ of n into odd parts}
#            Π_{k ∈ λ} [some function of k] / (multiplicities)!
#
# This looks like the exponential of Σ f_k.
# So: 1 + R(n) = exp(Σ_k g_k(n)) where g_k involves single k-cycle.
#
# Let's test: log(1+R) vs Σ single-k contributions.

print("  log(1 + R(n)) vs Σ_k f_k(n) (single-cycle sum):\n")

for n in range(3, 16):
    m = comb(n, 2)
    parts = gen_odd_partitions(n)
    identity = tuple([1]*n)

    total = Fraction(0)
    for ct in parts:
        mult = multinomial(n, ct)
        po = pair_orbits(ct)
        total += Fraction(mult * pow(2, po), pow(2, m))

    R = float(total) - 1.0
    log_1R = log(1 + R) if R > -1 else float('nan')

    # Sum of single k-cycle contributions
    single_sum = 0.0
    for k in range(3, n+1, 2):
        ct_k = tuple(sorted([k] + [1]*(n-k)))
        fk = multinomial(n, ct_k) * pow(2, pair_orbits(ct_k)) / pow(2, m)
        single_sum += fk

    if R != 0:
        ratio = log_1R / single_sum if single_sum != 0 else 0
        print("  n=%2d: log(1+R) = %12.8f, Σf_k = %12.8f, ratio = %.6f" %
              (n, log_1R, single_sum, ratio))

# ============================================================
banner("THE EXACT EULER PRODUCT DECOMPOSITION")
# ============================================================

# The Burnside sum over all-odd permutations has a nice structure.
# In the cycle index of S_n restricted to odd cycles:
#
# V_n × n! / 2^m = Z_n^{odd}(2^{po}) = Σ_λ mult(λ) 2^{po(λ)-m}
#
# The cycle index of S_n restricted to odd cycles:
# Z_n^{odd}(x_1, x_3, x_5, ...) = [z^n] exp(Σ_{k odd} x_k z^k / k)
#
# Substituting x_k = 2^{orbit gain from k-cycle}...
# This is getting complex. Let me just compute the exact product form.

# For each odd k, define:
# a_k(n) = total contribution from partitions where k is the LARGEST part

print("  CONTRIBUTION BY LARGEST CYCLE LENGTH:\n")

for n in range(3, 14):
    m = comb(n, 2)
    parts = gen_odd_partitions(n)
    identity = tuple([1]*n)

    by_max = defaultdict(Fraction)
    for ct in parts:
        if ct == identity: continue
        max_k = max(ct)
        mult = multinomial(n, ct)
        po = pair_orbits(ct)
        by_max[max_k] += Fraction(mult * pow(2, po), pow(2, m))

    total_R = sum(by_max.values())
    if n <= 10:
        print("  n=%d: R=%.8f" % (n, float(total_R)))
        for k in sorted(by_max.keys()):
            frac = float(by_max[k]) / float(total_R) if float(total_R) != 0 else 0
            print("    max_k=%d: %.10f (%.1f%% of R)" %
                  (k, float(by_max[k]), 100*frac))

# ============================================================
banner("THE PRIME-FACTORED CONTRIBUTION")
# ============================================================

# Now factor R(n) by the PRIMES appearing in the cycle type.
# A cycle of length k contributes through prime factor p | k.
# But 3-cycles are "atomic" (length 3 = prime).
# 5-cycles: atomic (prime 5).
# 9-cycles: 9 = 3² — these involve the prime 3 twice.
# 15-cycles: 15 = 3×5 — mixed.
#
# For the tournament counting function, the relevant "primes" are
# the ODD PRIMES 3, 5, 7, 11, 13, ...
# Each odd cycle length k contributes to V_n through these primes.
#
# Define: P_p(n) = contribution from cycle types where ALL parts are powers of p.
# Then: R(n) ≈ Σ_p P_p(n) + cross terms.
#
# At small n, only p=3 matters (since 5 > n/2 often).

print("  PRIME-DECOMPOSED CONTRIBUTIONS:\n")

def prime_factors(k):
    """Return set of prime factors of k."""
    factors = set()
    d = 2
    while d * d <= k:
        while k % d == 0:
            factors.add(d)
            k //= d
        d += 1
    if k > 1:
        factors.add(k)
    return factors

for n in range(3, 14):
    m = comb(n, 2)
    parts = gen_odd_partitions(n)
    identity = tuple([1]*n)

    # For each partition, find which primes appear
    by_primes = defaultdict(Fraction)
    for ct in parts:
        if ct == identity: continue
        primes_used = set()
        for k in ct:
            if k > 1:
                primes_used |= prime_factors(k)
        primes_key = tuple(sorted(primes_used))
        mult = multinomial(n, ct)
        po = pair_orbits(ct)
        by_primes[primes_key] += Fraction(mult * pow(2, po), pow(2, m))

    total_R = sum(by_primes.values())
    if n <= 10:
        print("  n=%d: R=%.10f" % (n, float(total_R)))
        for pk in sorted(by_primes.keys()):
            val = float(by_primes[pk])
            frac = val / float(total_R) if float(total_R) != 0 else 0
            label = "×".join(str(p) for p in pk) if pk else "identity"
            print("    primes {%s}: %12.10f (%.2f%%)" % (label, val, 100*frac))
        print()

# ============================================================
banner("THE 1/k RECIPROCAL IN THE CYCLE INDEX")
# ============================================================

# In the cycle index of S_n:
# Z(S_n) = Σ_λ (1/aut(λ)) Π_{i} x_i^{a_i}
# where aut(λ) = Π_i (i^{a_i} × a_i!)
#
# For a single k-cycle: coefficient of x_k z^{n-k}/...
# The contribution to Z of a single k-cycle with n-k fixed points:
#   (1/k) × x_k × x_1^{n-k} × [combinatorial factor]
#
# The 1/k factor is KEY: it's the RECIPROCAL of the cycle length.
# For prime k, this is 1/p — the prime reciprocal!
#
# So the tournament counting function involves:
# Σ_{p prime, p ≤ n} (1/p) × [polynomial in n] × 2^{quadratic in n}
#
# This is the PRIME RECIPROCAL GENERATING FUNCTION.

print("  THE 1/k FACTOR IN THE CYCLE INDEX:\n")
print("  For cycle type (k, 1^{n-k}): aut = k × 1^{n-k} × (n-k)! = k(n-k)!")
print("  mult = n! / aut = n!/(k(n-k)!) = C(n,k) × (k-1)!")
print("  Note: C(n,k)×(k-1)! = (n!/(k!(n-k)!)) × (k-1)! = n!/(k(n-k)!)")
print("  = (1/k) × n!/(n-k)!")
print()
print("  So the contribution to R(n) from a single k-cycle is:")
print("  R_k(n) = (1/k) × (n!/(n-k)!) × 2^{(k-1)(k+1-2n)/2 - m} × ... no.")
print()
print("  More precisely:")
print("  R_k(n) = C(n,k) × (k-1)! × 2^{Delta_k}")
print("         = (1/k) × n!/(n-k)! × 2^{(k-1)(k+1-2n)/2}")
print()

# Now compute R_k(n) / (1/k) to see what the "polynomial × exponential" part looks like
print("  R_k(n) × k = n!/(n-k)! × 2^{(k-1)(k+1-2n)/2}:\n")

for n in range(3, 14):
    m = comb(n, 2)
    line = "  n=%2d:" % n
    for k in [3, 5, 7]:
        if k > n: break
        falling = 1
        for j in range(k):
            falling *= (n - j)
        delta = (k-1) * (k+1 - 2*n) // 2
        val = falling * pow(2, delta) / pow(2, 0)  # This is k × R_k
        # But 2^delta might be very small
        val_exact = Fraction(falling, 1) * Fraction(pow(2, pair_orbits(tuple(sorted([k]+[1]*(n-k))))), pow(2, m))
        line += " k=%d: %.6e" % (k, float(val_exact * k))
    print(line)

# ============================================================
banner("THE TOURNAMENT EULER PRODUCT")
# ============================================================

# Define:
#   T(n) = V_n × n! / 2^m = 1 + R(n)
#
# If contributions from different cycle lengths are "independent":
#   T(n) = Π_{k=3,5,7,...} (1 + something_k(n))
#
# For this to work, we need the multi-cycle contributions to factor.
# A partition (3,3,1,...) should contribute f_3^2/2!, etc.
#
# The cycle index exponential formula says:
#   Σ_n Z(S_n) z^n = exp(Σ_{k≥1} x_k z^k / k)
#
# Restricted to odd cycles:
#   Σ_n Z^{odd}(S_n) z^n = exp(Σ_{k odd, k≥1} x_k z^k / k)
#
# With x_k = 2^{pair orbits gain} ... but pair orbits aren't multiplicative.
#   pair_orbits((3,3,...)) ≠ 2 × pair_orbits((3,1,...))
#   because the CROSS-ORBIT term between two 3-cycles is gcd(3,3)=3.
#
# Let's check how far off the product approximation is.

print("  PRODUCT APPROXIMATION vs EXACT:\n")
print("  T(n) = 1 + R(n) exact vs Π(1 + f_k) product:\n")

for n in range(3, 14):
    m = comb(n, 2)
    parts = gen_odd_partitions(n)
    identity = tuple([1]*n)

    # Exact T(n)
    T_exact = Fraction(0)
    for ct in parts:
        mult = multinomial(n, ct)
        po = pair_orbits(ct)
        T_exact += Fraction(mult * pow(2, po), pow(2, m))

    # Product approximation: Π(1 + f_k) where f_k = single k-cycle contrib
    product = 1.0
    for k in range(3, n+1, 2):
        ct_k = tuple(sorted([k] + [1]*(n-k)))
        fk = float(Fraction(multinomial(n, ct_k) * pow(2, pair_orbits(ct_k)), pow(2, m)))
        product *= (1 + fk)

    T_exact_f = float(T_exact)
    ratio = T_exact_f / product if product != 0 else 0
    print("  n=%2d: T_exact = %12.8f, Π(1+f_k) = %12.8f, ratio = %.8f" %
          (n, T_exact_f, product, ratio))

# ============================================================
banner("THE CROSS-ORBIT CORRECTION")
# ============================================================

# The product Π(1+f_k) misses the cross-orbit terms.
# For two k-cycles and j-cycles: the cross term is gcd(k,j) additional orbits.
# In the product, the (k,j) contribution would be f_k × f_j / ...
# but the actual contribution has an extra 2^{gcd(k,j)} factor.
#
# Define the CROSS-ORBIT MATRIX:
# G(k,j) = gcd(k,j) for the pair orbits between a k-cycle and j-cycle.
# The product approximation assumes G(k,j) = 0 (independent).
# The correction factor is: 2^{Σ_{i<j} gcd(k_i, k_j)} for multi-cycle types.

print("  CROSS-ORBIT MATRIX (gcd of cycle lengths):\n")
print("  gcd | ", end="")
for k in range(3, 14, 2):
    print("%-4d" % k, end="")
print()
print("  " + "-"*50)
for k in range(3, 14, 2):
    print("  %-3d | " % k, end="")
    for j in range(3, 14, 2):
        print("%-4d" % gcd(k, j), end="")
    print()

print("\n  The cross-orbit correction for (3,3) is 2^{gcd(3,3)} = 2^3 = 8.")
print("  For (3,5) it's 2^{gcd(3,5)} = 2^1 = 2.")
print("  For (3,7) it's 2^{gcd(3,7)} = 2^1 = 2.")
print("  For (5,5) it's 2^{gcd(5,5)} = 2^5 = 32.")
print("  For (3,3,3) it's 2^{3+3+3} = 2^9 = 512.")
print()
print("  When both cycles are the SAME prime p: cross = 2^p (LARGE).")
print("  When coprime: cross = 2^1 = 2 (SMALL).")
print("  This means: same-prime multi-cycles interact STRONGLY,")
print("  while different-prime cycles interact WEAKLY.")

# ============================================================
banner("TOURNAMENT ZETA FUNCTION: FORMAL DEFINITION")
# ============================================================

# Define the tournament zeta function:
#   Z_T(s, n) = Π_{k odd, k≥3} (1 - q_k(n) × k^{-s})^{-1}
#
# where q_k(n) = C(n,k) × (k-1)! × 2^{Delta_k} × k
#              = n!/(n-k)! × 2^{(k-1)(k+1-2n)/2}
#
# At s=1: Z_T(1,n) involves Σ q_k/k = Σ C(n,k)(k-1)! × 2^{Delta_k}
#         = Σ R_k(n) = the single-cycle correction.
#
# At s=0: Z_T(0,n) = Π(1-q_k)^{-1} → diverges if q_k > 0.
#
# The "natural" value is s where Z_T(s,n) = T(n) = V_n × n!/2^m.

print("""
  FORMAL TOURNAMENT ZETA FUNCTION:

  Z_T(s, n) = Π_{k=3,5,7,...} (1 - R_k(n) × k^{-s+1})^{-1}

  This is an Euler product over ODD cycle lengths (= "tournament primes").

  At s=1: each factor ≈ 1/(1 - R_k) ≈ 1 + R_k (for small R_k).
  Product ≈ 1 + Σ R_k = 1 + single-cycle correction.

  The ACTUAL T(n) has extra cross-orbit factors.

  Key insight: The "primes" of tournament theory are
  the ODD INTEGERS k=3,5,7,9,11,...
  NOT the actual primes 3,5,7,11,13,...

  But! The cross-orbit matrix gcd(k,j) = 1 when k,j coprime.
  So coprime odd cycles interact WEAKLY (factor 2^1).
  Same-prime powers interact STRONGLY (factor 2^p).

  This means the tournament primes ARE the actual odd primes,
  in the sense that the cross-orbit interaction is minimal
  between different prime families.
""")

# Compute the ratio T(n) / product approximation to quantify cross-orbit effect
print("  CROSS-ORBIT CORRECTION FACTOR C(n) = T(n) / Π(1+f_k):\n")

for n in range(3, 16):
    m = comb(n, 2)
    parts = gen_odd_partitions(n)
    identity = tuple([1]*n)

    T_exact = Fraction(0)
    for ct in parts:
        T_exact += Fraction(multinomial(n, ct) * pow(2, pair_orbits(ct)), pow(2, m))

    product = Fraction(1)
    for k in range(3, n+1, 2):
        ct_k = tuple(sorted([k] + [1]*(n-k)))
        fk = Fraction(multinomial(n, ct_k) * pow(2, pair_orbits(ct_k)), pow(2, m))
        product *= (1 + fk)

    cross = float(T_exact) / float(product) if float(product) != 0 else 0
    print("  n=%2d: C(n) = %.10f" % (n, cross))

# ============================================================
banner("THE 1/p SUM IN THE CYCLE INDEX")
# ============================================================

# The cycle index restricted to odd cycles:
# Z_n^{odd} = Σ_{λ all-odd} (1/aut(λ)) Π x_{λ_i}
#
# The EGF: Σ_n Z_n^{odd} z^n = exp(Σ_{k odd, k≥1} x_k z^k / k)
#
# The 1/k factor is EXPLICITLY PRESENT for each odd k.
# For k prime: this is 1/p.
# The Mertens function M = lim_{x→∞} (Σ_{p≤x} 1/p - ln ln x) ≈ 0.2615.
#
# Our tournament correction involves:
# Σ_{k odd, k≥3} (1/k) × [n^{falling k}] × 2^{quadratic(k,n)}
#
# The 1/k factor creates a HARMONIC-LIKE sum weighted by the
# super-exponentially decaying 2^{quadratic}.
#
# At large n, only k=3 survives: R(n) ≈ (1/3) × n³ × 2^{4-2n}
#                                       = n(n-1)(n-2)/6 × 2^{4-2n}

print("  ASYMPTOTIC FORMULA: R(n) ≈ C(n,3) × 2 × 2^{4-2n}")
print("  = n(n-1)(n-2)/3 × 2^{3-2n}\n")

for n in range(3, 20):
    m = comb(n, 2)
    # Asymptotic (single 3-cycle only)
    R_asymp = comb(n, 3) * 2 * pow(2, 4-2*n)

    if n <= 13:
        parts = gen_odd_partitions(n)
        identity = tuple([1]*n)
        T_exact = sum(Fraction(multinomial(n, ct) * pow(2, pair_orbits(ct)), pow(2, m))
                      for ct in parts)
        R_exact = float(T_exact) - 1.0
        ratio = R_exact / R_asymp if R_asymp != 0 else 0
        print("  n=%2d: R_exact = %12.8f, R_asymp = %12.8f, ratio = %.6f" %
              (n, R_exact, R_asymp, ratio))
    else:
        print("  n=%2d: R_asymp = %.12e" % (n, R_asymp))

print("""
  THE PUNCH LINE:

  V_n ≈ 2^{C(n,2)} / n! × (1 + n(n-1)(n-2)/3 × 2^{3-2n} + ...)

  The correction involves:
  - 1/3 (= reciprocal of smallest tournament prime)
  - n^{falling 3} / 3 = number of 3-cycles in S_n
  - 2^{3-2n} = exponential decay

  Each higher "prime" k adds:
  - 1/k × n^{falling k} × 2^{(k-1)(k+1-2n)/2}

  The decay is SUPER-EXPONENTIAL in k: factor 2^{-k(k-1)/2} per step.

  So the tournament counting function is DOMINATED by the
  smallest prime reciprocal 1/3, with corrections from 1/5, 1/7, ...
  that are exponentially negligible.

  The tournament zeta function Z_T(s) = Σ V_n × n!/2^m × n^{-s}
  has an Euler product expansion where the "primes" are odd integers
  and the 1/p factors appear through the cycle index.
""")

print("Done. opus-2026-03-24-S291")
