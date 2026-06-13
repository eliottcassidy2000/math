#!/usr/bin/env python3
"""
opus-2026-04-05-S24c: MERGING 2-adic valuations of a(n) with QR-tournament theory.

Goals:
  1. Extend v₂(a(n)) table to n=30 using known A000568 values
  2. PROVE HYP-1712 (even-n v₂ conjecture) via Burnside partition analysis
  3. Connect v₂(a(n)) to orbit parity (HYP-1713)
  4. Analyze the odd part T(n)/2^{(n-1)/2}
  5. Prove HYP-1713 (orbit parity) via anti-automorphism

The key bridge: both HYP-1712 and HYP-1713 ultimately depend on
understanding which partition of n minimizes v₂ in a Burnside-type sum,
and how the anti-automorphism τ interacts with the 2-adic structure.
"""

import sys
from math import factorial, comb, gcd, log2
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
import time

sys.stdout.reconfigure(line_buffering=True)

# ═══════════════════════════════════════════════════════════════════
# Known A000568 values (from OEIS + our burnside_enum_v2.c extensions)
# ═══════════════════════════════════════════════════════════════════
A000568 = {
    1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456,
    8: 6880, 9: 191536, 10: 9733056, 11: 903753248,
    12: 154108311168, 13: 48542114686912,
    14: 28401423719122304, 15: 31021002160355166848,
    16: 63530415842308265100288,
    17: 244912778438520759443245824,
    18: 1783398846284777975419600287232,
    19: 24605641171260376770598003978281472,
    20: 645022068557873570931850526424042500096,
}

# ═══════════════════════════════════════════════════════════════════
# Utilities
# ═══════════════════════════════════════════════════════════════════

def v2(n):
    """2-adic valuation of n."""
    if n == 0:
        return float('inf')
    v = 0
    while n % 2 == 0:
        v += 1
        n //= 2
    return v

def odd_part(n):
    """n / 2^{v₂(n)}."""
    while n % 2 == 0:
        n //= 2
    return n

def euler_phi(n):
    """Euler's totient function."""
    result = n
    p = 2
    temp = n
    while p * p <= temp:
        if temp % p == 0:
            while temp % p == 0:
                temp //= p
            result -= result // p
        p += 1
    if temp > 1:
        result -= result // temp
    return result

def legendre(a, p):
    """Legendre symbol (a/p)."""
    if a % p == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1

def is_prime(n):
    if n < 2: return False
    if n < 4: return True
    if n % 2 == 0 or n % 3 == 0: return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def factor(n):
    """Return list of prime factors with multiplicity."""
    if n <= 1: return []
    factors = []
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors.append(d)
            n //= d
        d += 1
    if n > 1:
        factors.append(n)
    return factors

def all_odd_partitions(n, max_part=None):
    """Generate all partitions of n into odd parts, in non-increasing order."""
    if max_part is None:
        max_part = n
    if max_part % 2 == 0:
        max_part -= 1
    if n == 0:
        yield ()
        return
    if max_part < 1:
        return
    start = min(max_part, n)
    if start % 2 == 0:
        start -= 1  # ensure we start from an odd number
    for first in range(start, 0, -2):  # odd parts, decreasing
        for rest in all_odd_partitions(n - first, first):
            yield (first,) + rest

def edge_orbits(partition, n):
    """Number of orbits of the action of a permutation with given cycle type
    on the set of 2-element subsets of {1,...,n}.

    For cycle type (c_1, c_2, ..., c_k):
      orbits = sum_{i<j} gcd(c_i, c_j) + sum_i floor(c_i/2)
    """
    total = 0
    parts = list(partition)
    k = len(parts)
    # Pairs within same cycle: floor(c/2) orbits per cycle
    for c in parts:
        total += c // 2
    # Pairs across different cycles: gcd(c_i, c_j) orbits
    for i in range(k):
        for j in range(i + 1, k):
            total += gcd(parts[i], parts[j])
    return total

def centralizer_order(partition):
    """z_λ = prod(c_i) * prod(a_j!) where a_j is the multiplicity of part j."""
    prod_parts = 1
    for c in partition:
        prod_parts *= c
    mult = Counter(partition)
    prod_mult = 1
    for a in mult.values():
        prod_mult *= factorial(a)
    return prod_parts * prod_mult

print("=" * 78)
print("  MERGING 2-ADIC VALUATIONS WITH QR-TOURNAMENT THEORY")
print("  opus-2026-04-05-S24c")
print("=" * 78)

# ═══════════════════════════════════════════════════════════════════
# SECTION 1: Extended v₂(a(n)) table
# ═══════════════════════════════════════════════════════════════════
print("\n" + "═" * 78)
print("  SECTION 1: v₂(a(n)) FOR n = 1..20")
print("═" * 78)

print(f"\n  {'n':>3} {'a(n)':>30} {'v₂':>5} {'odd part':>25} {'n odd?':>6} "
      f"{'(n-1)/2':>7} {'pred':>5} {'match':>6}")
print("  " + "-" * 90)

for n in range(1, 21):
    a = A000568[n]
    v = v2(a)
    op = a >> v  # odd part

    if n % 2 == 1:
        pred = (n - 1) // 2
    else:
        # HYP-1712: v₂ = n/2 + v₂(c(n))
        # c(n) = #{odd k : 1 <= k < n/2, gcd(k,n) = 1}
        c_n = sum(1 for k in range(1, n // 2) if k % 2 == 1 and gcd(k, n) == 1)
        pred = n // 2 + v2(c_n) if c_n > 0 else n // 2

    match = "✓" if v == pred else "✗"

    # Truncate odd part for display
    op_str = str(op)
    if len(op_str) > 25:
        op_str = op_str[:10] + "..." + op_str[-10:]

    print(f"  {n:>3} {str(a):>30} {v:>5} {op_str:>25} "
          f"{'ODD' if n % 2 == 1 else 'even':>6} {(n-1)//2 if n%2==1 else '':>7} "
          f"{pred:>5} {match:>6}")

# ═══════════════════════════════════════════════════════════════════
# SECTION 2: PROVE THM-305 and HYP-1712 via Burnside partition analysis
# ═══════════════════════════════════════════════════════════════════
print("\n\n" + "═" * 78)
print("  SECTION 2: BURNSIDE PARTITION ANALYSIS — THE MINIMIZING PARTITION")
print("═" * 78)

print("""
  The Burnside formula: a(n) = (1/n!) × Σ_{λ all-odd ⊢ n} (n!/z_λ) × 2^{e(λ)}

  where e(λ) = edge_orbits(λ) and z_λ = centralizer_order(λ).

  Each term contributes: n!/z_λ × 2^{e(λ)} to n!·a(n).

  The v₂ of each term = v₂(n!/z_λ) + e(λ).

  Since v₂(n!) is fixed, the minimizer of v₂(term) maximizes
  v₂(z_λ) - e(λ), or equivalently minimizes e(λ) - v₂(z_λ).

  For ODD n:
  THM-305 shows λ = (n) minimizes with e = (n-1)/2, z = n, v₂(z) = 0.
  So min = e(n) - v₂(z_{(n)}) = (n-1)/2 - 0 = (n-1)/2.

  For EVEN n:
  The partition (n) is NOT all-odd. The candidate minimizers are:
  - λ = (n-1, 1): e = (n-2)/2 + 0 + gcd(n-1,1) = (n-2)/2 + 1 = n/2.
    z = (n-1)·1·1!·1! = n-1. v₂(z) = v₂(n-1).
    Contribution v₂ = v₂(n!) - v₂(n-1) + n/2 = v₂(n!) - v₂(n-1) + n/2.

  - λ = (1^n): e = C(n,2)/1 ... no. For (1^n), there are n fixed points.
    e = 0 (no pairs within a 1-cycle) + Σ_{i<j} gcd(1,1) = C(n,2).
    z = 1^n · n! = n!. v₂(z) = v₂(n!).
    Contribution: n!/n! × 2^{C(n,2)} = 2^{C(n,2)}.
    v₂ = C(n,2).

  So (n-1,1) has MUCH smaller e than (1^n). Let's check systematically.
""")

for n in [4, 6, 8, 10, 12, 14]:
    print(f"\n  n = {n}: all-odd partitions and their v₂ contributions:")
    print(f"    {'partition':>20} {'e(λ)':>6} {'z_λ':>12} {'v₂(z_λ)':>8} "
          f"{'e-v₂(z)':>8} {'v₂(n!/z)':>8} {'v₂(term)':>8}")
    print("    " + "-" * 76)

    v2_nfact = v2(factorial(n))
    terms = []

    for lam in all_odd_partitions(n):
        e = edge_orbits(lam, n)
        z = centralizer_order(lam)
        v2_z = v2(z)
        coeff = factorial(n) // z
        v2_coeff = v2(coeff)
        v2_term = v2_coeff + e  # = v₂(n!/z) + e = v₂(n!) - v₂(z) + e

        terms.append((v2_term, lam, e, z, v2_z))

    # Sort by v₂(term)
    terms.sort()

    # Show top 5 smallest
    for rank, (v2t, lam, e, z, v2z) in enumerate(terms[:8]):
        lam_str = str(lam)
        if len(lam_str) > 20:
            lam_str = lam_str[:17] + "..."
        v2_coeff = v2(factorial(n) // z)
        print(f"    {lam_str:>20} {e:>6} {z:>12} {v2z:>8} {e - v2z:>8} "
              f"{v2_coeff:>8} {v2t:>8}"
              f"{'  ← MIN' if rank == 0 else ''}")

    # The minimum
    min_v2 = terms[0][0]
    # a(n) = sum of terms / n!, so v₂(n!·a(n)) = min_v2 (if unique minimum)
    # v₂(a(n)) = min_v2 - v₂(n!)

    v2_an = min_v2 - v2_nfact
    actual_v2 = v2(A000568[n])

    # How many terms achieve the minimum?
    min_count = sum(1 for v, _, _, _, _ in terms if v == min_v2)
    second = terms[min_count][0] if min_count < len(terms) else "N/A"
    gap = second - min_v2 if isinstance(second, int) else "N/A"

    print(f"\n    v₂(n!) = {v2_nfact}, min v₂(term) = {min_v2}")
    print(f"    → v₂(a(n)) ≥ {v2_an}")
    print(f"    Actual v₂(a(n)) = {actual_v2}")
    print(f"    Terms at minimum: {min_count}, gap to next: {gap}")
    print(f"    Minimizing partition: {terms[0][1]}")

    # Check: is the minimizer always (n-1, 1)?
    if terms[0][1] == (n - 1, 1):
        print(f"    ✓ Minimizer is (n-1, 1) as predicted!")
    else:
        print(f"    ✗ Minimizer is {terms[0][1]}, NOT (n-1, 1)!")

# ═══════════════════════════════════════════════════════════════════
# SECTION 3: THE PROOF OF HYP-1712
# ═══════════════════════════════════════════════════════════════════
print("\n\n" + "═" * 78)
print("  SECTION 3: PROVING HYP-1712 — v₂(a(n)) FOR EVEN n")
print("═" * 78)

print("""
  THEOREM (HYP-1712 → THM): For even n ≥ 4:
    v₂(a(n)) = n/2 + v₂(c(n))

  where c(n) = #{odd k : 1 ≤ k < n/2, gcd(k,n) = 1}.

  PROOF STRUCTURE:

  Step 1: The Burnside sum for a(n) involves only all-odd partitions of n.
  For even n, these have an even number of parts (since sum of odd #parts = odd).

  Step 2: The partition λ = (n-1, 1) has:
    - e(λ) = (n-1-1)/2 + gcd(n-1, 1) = (n-2)/2 + 1 = n/2
    - z_λ = (n-1) · 1 · 1! · 1! = n-1
    - v₂(z_λ) = v₂(n-1)
    - v₂(term) = v₂(n!) - v₂(n-1) + n/2

  Step 3: For ALL other all-odd partitions λ' ≠ (n-1,1) of even n:
    v₂(term(λ')) > v₂(term((n-1,1))) + gap

  This needs to be verified. The critical competitor is (n-3, 1, 1, 1):
    - e = (n-4)/2 + 3·gcd(n-3,1) + 3·gcd(1,1) = (n-4)/2 + 3 + 3 = (n-4)/2 + 6
    Wait, that's wrong. Let me recompute.

  For λ = (n-3, 1, 1, 1):
    - e = floor((n-3)/2) + 0 + 0 + 0 [within-cycle for each part]
         + gcd(n-3,1) + gcd(n-3,1) + gcd(n-3,1) + gcd(1,1) + gcd(1,1) + gcd(1,1)
         [cross-cycle pairs]
    - = (n-4)/2 + 3·1 + 3·1 = (n-4)/2 + 6
    Wait, gcd(1,1) = 1, and there are C(3,2) = 3 pairs among the three 1-parts.
    And there are 3 pairs of (n-3)-part with each 1-part.
    - e = (n-4)/2 + 3 + 3 = (n-4)/2 + 6 = n/2 + 4
    - z = (n-3)·1·1·1 · 1!·3! = (n-3)·6
    - v₂(z) = v₂(n-3) + v₂(6) = v₂(n-3) + 1
    - v₂(term) = v₂(n!) - v₂(z) + e = v₂(n!) - v₂(n-3) - 1 + n/2 + 4

  Compare with (n-1,1):
    - v₂(term) = v₂(n!) - v₂(n-1) + n/2

  Difference: [v₂(n-1) - v₂(n-3)] + [4 - 1] = v₂(n-1) - v₂(n-3) + 3

  For n ≥ 4: v₂(n-1) ≥ 0 and v₂(n-3) ≥ 0.
  If n-1 and n-3 are both odd (n even, n ≡ 0 mod 2):
    n-1 odd ⟺ n even ✓
    n-3 odd ⟺ n even ✓
  So v₂(n-1) = 0 and v₂(n-3) = 0, giving difference = 0 + 3 = 3.

  So (n-3,1,1,1) has v₂ exactly 3 more than (n-1,1). The gap is ≥ 3.

  But wait — I said v₂(n-1) = 0 and v₂(n-3) = 0 for even n. That's because
  n-1 and n-3 are odd when n is even. So the (n-1,1) term has
  v₂(term) = v₂(n!) + n/2, and the (n-3,1,1,1) term has
  v₂(term) = v₂(n!) - 1 + n/2 + 4 = v₂(n!) + n/2 + 3.

  Now (n-1,1) itself: since n is even, n-1 is odd, v₂(n-1) = 0.
  So v₂(term) = v₂(n!) - 0 + n/2 = v₂(n!) + n/2.
  And v₂(a(n)) = v₂(term) - v₂(n!) = n/2.

  But HYP-1712 says v₂(a(n)) = n/2 + v₂(c(n)), NOT just n/2.
  The extra v₂(c(n)) comes from MULTIPLE terms achieving the same minimum!

  When multiple terms have v₂ = v₂(n!) + n/2, their sum might have
  higher v₂ if they cancel at the binary level.
""")

# Let's identify ALL terms at the minimum level for each even n
print("\n  Detailed minimum-level analysis for even n:\n")

for n in range(4, 22, 2):
    if n not in A000568:
        break

    v2_nfact = v2(factorial(n))
    min_terms = []
    all_terms_data = []

    for lam in all_odd_partitions(n):
        e = edge_orbits(lam, n)
        z = centralizer_order(lam)
        coeff = factorial(n) // z
        v2_term = v2(coeff) + e

        all_terms_data.append((v2_term, lam, coeff, e, z))

    all_terms_data.sort()
    min_v2_term = all_terms_data[0][0]

    # Collect all terms at the minimum level
    for v2t, lam, coeff, e, z in all_terms_data:
        if v2t == min_v2_term:
            # The actual coefficient in the Burnside sum (before dividing by n!)
            # is (n!/z_λ) × 2^{e(λ)}
            # At the minimum v₂ level, each such term contributes
            # odd_part(n!/z_λ) × 2^{e(λ) + v₂(n!/z_λ) - min_v2_term} ...
            # Actually, let's compute the odd coefficient at the minimum level.
            #
            # The term is (n!/z_λ) × 2^{e(λ)}.
            # v₂(this) = v₂(n!/z_λ) + e(λ) = min_v2_term.
            # So the odd part is (n!/z_λ)/2^{v₂(n!/z_λ)} = odd_part(n!/z_λ).
            # And 2^{e(λ) + v₂(n!/z_λ) - min_v2_term} = 2^0 = 1.
            # So the coefficient at the minimum level is odd_part(n!/z_λ).
            odd_coeff = odd_part(coeff)
            min_terms.append((lam, odd_coeff, coeff, e, z))

    # c(n) as defined in HYP-1712
    c_n = sum(1 for k in range(1, n // 2) if k % 2 == 1 and gcd(k, n) == 1)

    # Sum of odd coefficients at minimum level
    sum_odd = sum(oc for _, oc, _, _, _ in min_terms)
    v2_sum = v2(sum_odd)

    # Predicted v₂(a(n))
    predicted_v2 = n // 2 + (v2(c_n) if c_n > 0 else 0)
    actual_v2 = v2(A000568[n])

    # v₂(a(n)) = min_v2_term - v₂(n!) + v₂(sum_odd)
    computed_v2 = min_v2_term - v2_nfact + v2_sum

    print(f"  n={n:>2}: {len(min_terms)} terms at min level, "
          f"sum_odd_coeffs = {sum_odd}, v₂(sum) = {v2_sum}, "
          f"c(n) = {c_n}, v₂(c) = {v2(c_n) if c_n > 0 else 'N/A'}")
    print(f"        v₂(a(n)) = {min_v2_term} - {v2_nfact} + {v2_sum} "
          f"= {computed_v2}, actual = {actual_v2}, "
          f"predicted = {predicted_v2} {'✓' if computed_v2 == actual_v2 == predicted_v2 else '✗'}")

    if len(min_terms) <= 10:
        for lam, oc, coeff, e, z in min_terms:
            print(f"        λ={lam}, odd_coeff={oc}, e={e}, z={z}")
    else:
        # Show first 3 and last 2
        for lam, oc, coeff, e, z in min_terms[:3]:
            print(f"        λ={lam}, odd_coeff={oc}, e={e}, z={z}")
        print(f"        ... ({len(min_terms) - 5} more) ...")
        for lam, oc, coeff, e, z in min_terms[-2:]:
            print(f"        λ={lam}, odd_coeff={oc}, e={e}, z={z}")
    print()

# ═══════════════════════════════════════════════════════════════════
# SECTION 4: THE KEY INSIGHT — minimum level terms and c(n)
# ═══════════════════════════════════════════════════════════════════
print("\n" + "═" * 78)
print("  SECTION 4: WHY sum_odd = c(n) × (odd factor)")
print("═" * 78)

print("""
  For even n, the minimum-level partitions all have the form (n-1, 1) up to
  automorphism. But ARE there other partitions at the same minimum v₂ level?

  The partition (n-1, 1) has:
    v₂(term) = v₂(n!) + n/2  (since v₂(n-1) = 0 for even n)

  Other 2-part partitions (a, b) with a+b = n, a,b odd, a > b ≥ 1:
    e(a,b) = (a-1)/2 + (b-1)/2 + gcd(a,b) = (a+b-2)/2 + gcd(a,b) = (n-2)/2 + gcd(a,b)
    z(a,b) = a·b  (if a ≠ b) or a²·2! = 2a² (if a = b)

  For a ≠ b: v₂(term) = v₂(n!) - v₂(a·b) + (n-2)/2 + gcd(a,b)
  Since a,b odd: v₂(a·b) = 0.
  So v₂(term) = v₂(n!) + (n-2)/2 + gcd(a,b) = v₂(n!) + n/2 - 1 + gcd(a,b).

  This equals the minimum v₂(n!) + n/2 iff gcd(a,b) = 1.

  So the MINIMUM-LEVEL 2-PART PARTITIONS are exactly:
    (a, b) with a + b = n, a > b ≥ 1, a,b odd, gcd(a,b) = 1.

  The number of such partitions is:
    #{b : 1 ≤ b < n/2, b odd, gcd(b, n-b) = 1}
  = #{b : 1 ≤ b < n/2, b odd, gcd(b, n) = 1}  (since gcd(b,n-b) = gcd(b,n))

  THIS IS EXACTLY c(n)!

  So the minimum-level sum has c(n) terms, each contributing an odd coefficient.
  The v₂ of their sum is v₂(c(n)) plus v₂(average odd coefficient)...
  but the coefficients are not all 1. We need to check whether they are all
  ≡ 1 mod 2 (which would give v₂(sum) = v₂(c(n)) since c(n) odd × something).
""")

# Verify: are the minimum-level partitions exactly the coprime 2-part ones?
print("  Verification:\n")
for n in range(4, 22, 2):
    if n not in A000568: break

    v2_nfact = v2(factorial(n))

    # Find minimum v₂(term)
    min_v2_term = float('inf')
    for lam in all_odd_partitions(n):
        e = edge_orbits(lam, n)
        z = centralizer_order(lam)
        coeff = factorial(n) // z
        v2_term = v2(coeff) + e
        if v2_term < min_v2_term:
            min_v2_term = v2_term

    # Collect minimum-level partitions
    min_parts = []
    for lam in all_odd_partitions(n):
        e = edge_orbits(lam, n)
        z = centralizer_order(lam)
        coeff = factorial(n) // z
        v2_term = v2(coeff) + e
        if v2_term == min_v2_term:
            min_parts.append(lam)

    # Check: are they all 2-part?
    all_2part = all(len(lam) == 2 for lam in min_parts)

    # Count c(n) directly
    c_n = sum(1 for k in range(1, n // 2) if k % 2 == 1 and gcd(k, n) == 1)

    # The 2-part coprime partitions
    coprime_parts = [(n - k, k) for k in range(1, n // 2)
                     if k % 2 == 1 and gcd(k, n) == 1]

    match = set(min_parts) == set(coprime_parts)

    print(f"  n={n:>2}: min-level = {len(min_parts)} partitions, "
          f"c(n) = {c_n}, all 2-part? {all_2part}, "
          f"match coprime? {match}")

    if not match and len(min_parts) <= 20:
        extra = set(min_parts) - set(coprime_parts)
        missing = set(coprime_parts) - set(min_parts)
        if extra:
            print(f"        Extra min-level: {extra}")
        if missing:
            print(f"        Missing from min: {missing}")

# ═══════════════════════════════════════════════════════════════════
# SECTION 5: THE ODD COEFFICIENTS
# ═══════════════════════════════════════════════════════════════════
print("\n\n" + "═" * 78)
print("  SECTION 5: ODD COEFFICIENTS AT THE MINIMUM LEVEL")
print("═" * 78)

print("""
  For a 2-part partition (a, b) with a + b = n, a > b, gcd(a,b) = 1:
    n!/z_λ = n!/(a·b) = n!/(a·b)

  Since a, b are odd and coprime:
    odd_part(n!/z_λ) = odd_part(n!) / (a·b)
    (because a·b is odd, so v₂(n!/z) = v₂(n!))

  Wait: v₂(n!/(a·b)) = v₂(n!) - v₂(a·b) = v₂(n!) - 0 = v₂(n!).
  So the odd part of n!/(a·b) is n!/(2^{v₂(n!)} · a · b).

  The minimum-level coefficient for partition (a,b) with gcd(a,b)=1 is:
    odd_part(n!/(a·b)) = odd_part(n!) / (a·b)

  The SUM of all minimum-level coefficients is:
    Σ_{(a,b)} odd_part(n!) / (a·b)  where the sum is over coprime pairs

  Since odd_part(n!) is a common factor, v₂ of the sum is:
    v₂(sum) = v₂(Σ_{(a,b)} 1/(a·b) × odd_part(n!))
            = v₂(odd_part(n!)) + v₂(Σ 1/(a·b))
  But odd_part(n!) is ODD, so v₂(odd_part(n!)) = 0.
  And 1/(a·b) is not an integer...

  Let me reconsider. The sum of the n!/(a·b) terms at the minimum level is:
    S = Σ_{(a,b): a+b=n, gcd(a,b)=1, a,b odd, a>b} n!/(a·b)

  Factor out n!: S = n! × Σ 1/(a·b).
  This is an integer because n!/(a·b) is always an integer.

  v₂(S) = v₂(n!) + v₂(Σ 1/(a·b))

  The full term is S × 2^{n/2} (since e(λ) = n/2 for coprime 2-part).
  Wait, actually e(λ) = (n-2)/2 + gcd(a,b) = (n-2)/2 + 1 = n/2.
  The term in n!·a(n) is: (n!/z_λ) × 2^{e(λ)} = (n!/(a·b)) × 2^{n/2}.

  Sum over all minimum-level: S × 2^{n/2} where S = n! × Σ 1/(a·b).

  v₂(S × 2^{n/2}) = v₂(S) + n/2.

  Then v₂(a(n)) = v₂(S × 2^{n/2} + higher) - v₂(n!)
               = v₂(S) + n/2 - v₂(n!)  (if higher terms don't affect)
               = v₂(Σ n!/(a·b)) + n/2 - v₂(n!)
               = [v₂(n!) + v₂(Σ 1/(a·b))] + n/2 - v₂(n!)
               = v₂(Σ 1/(a·b)) + n/2

  So HYP-1712 reduces to: v₂(Σ_{coprime (a,b)} 1/(a·b)) = v₂(c(n)).
  i.e., v₂(Σ 1/(a·(n-a))) = v₂(c(n)) where sum is over odd a < n/2 with gcd(a,n) = 1.
""")

# Verify this reduction
print("  Verification of the reduction v₂(Σ 1/(a(n-a))) = v₂(c(n)):\n")

for n in range(4, 22, 2):
    # Compute Σ 1/(a·(n-a)) as a fraction
    total = Fraction(0)
    count = 0
    for a in range(1, n // 2):
        if a % 2 == 1 and gcd(a, n) == 1:
            b = n - a
            total += Fraction(1, a * b)
            count += 1

    c_n = count
    v2_cn = v2(c_n) if c_n > 0 else 0

    # v₂ of the fraction total = v₂(numerator) - v₂(denominator)
    if total.numerator != 0:
        v2_total = v2(total.numerator) - v2(total.denominator)
    else:
        v2_total = float('inf')

    print(f"  n={n:>2}: c(n) = {c_n}, v₂(c) = {v2_cn}, "
          f"Σ 1/(a(n-a)) = {total}, v₂ = {v2_total}, "
          f"match: {'✓' if v2_total == v2_cn else '✗'}")

# ═══════════════════════════════════════════════════════════════════
# SECTION 6: v₂(Σ 1/(a(n-a))) = v₂(c(n)) — THE ALGEBRAIC PROOF
# ═══════════════════════════════════════════════════════════════════
print("\n\n" + "═" * 78)
print("  SECTION 6: PROVING v₂(Σ 1/(a(n-a))) = v₂(c(n))")
print("═" * 78)

print("""
  We need: v₂(Σ_{a ∈ C} 1/(a(n-a))) = v₂(|C|)
  where C = {a : 1 ≤ a < n/2, a odd, gcd(a,n) = 1}.

  Write Σ 1/(a(n-a)) = (1/n) Σ (1/a + 1/(n-a)) = (1/n) Σ n/(a(n-a))
  Actually: 1/(a(n-a)) = (1/n)(1/a + 1/(n-a)).

  So Σ 1/(a(n-a)) = (1/n) Σ_{a ∈ C} (1/a + 1/(n-a)).

  Since n-a runs over {n-a : a ∈ C}, and n-a is also odd (since n even, a odd)
  and gcd(n-a, n) = gcd(a, n) = 1, the set {n-a : a ∈ C} is a DIFFERENT
  subset of the odd residues coprime to n.

  Actually, if a < n/2, then n-a > n/2, so {n-a : a ∈ C} are the elements
  in (n/2, n) that are odd and coprime to n. Together, C and {n-a : a ∈ C}
  form ALL odd numbers in [1, n-1] coprime to n.

  So Σ (1/a + 1/(n-a)) = Σ_{a ∈ C} 1/a + Σ_{a ∈ C} 1/(n-a)
                        = Σ_{d : 1 ≤ d ≤ n-1, d odd, gcd(d,n)=1} 1/d

  And Σ 1/(a(n-a)) = (1/n) × Σ_{d ∈ D} 1/d
  where D = {d : 1 ≤ d ≤ n-1, d odd, gcd(d,n)=1}, |D| = 2|C| = 2c(n).

  So the question reduces to: v₂(Σ_{d ∈ D} 1/d) = v₂(n) + v₂(c(n)).
  Wait, let me redo: v₂(Σ 1/(a(n-a))) = v₂(1/n × Σ 1/d) = -v₂(n) + v₂(Σ 1/d).
  We want this = v₂(c(n)).
  So we need: v₂(Σ_{d ∈ D} 1/d) = v₂(n) + v₂(c(n)).

  But v₂(n) ≥ 1 (n is even). And Σ 1/d is a fraction.
  Hmm, this is getting complicated. Let me just verify numerically.
""")

# Verify the partial fraction decomposition
print("  Detailed verification:\n")
for n in range(4, 16, 2):
    C = [a for a in range(1, n // 2) if a % 2 == 1 and gcd(a, n) == 1]
    D = [d for d in range(1, n) if d % 2 == 1 and gcd(d, n) == 1]

    sum_reciprocal_D = sum(Fraction(1, d) for d in D)
    sum_product = sum(Fraction(1, a * (n - a)) for a in C)
    pf_check = sum_product * n  # should equal sum_reciprocal_D

    print(f"  n={n}: C = {C}, |C| = {len(C)}")
    print(f"        D = {D}, |D| = {len(D)}")
    print(f"        Σ_D 1/d = {sum_reciprocal_D}")
    print(f"        Σ_C 1/(a(n-a)) = {sum_product}")
    print(f"        n × Σ_C = {pf_check}, equals Σ_D? {pf_check == sum_reciprocal_D}")
    print(f"        v₂(Σ_C) = {v2(sum_product.numerator) - v2(sum_product.denominator)}, "
          f"v₂(c(n)) = {v2(len(C))}")
    print()

# ═══════════════════════════════════════════════════════════════════
# SECTION 7: ODD PARTS T(n)/2^{(n-1)/2} FOR ODD n
# ═══════════════════════════════════════════════════════════════════
print("\n" + "═" * 78)
print("  SECTION 7: ODD PARTS a(n)/2^{(n-1)/2} FOR ODD n")
print("═" * 78)

print(f"\n  {'n':>3} {'a(n)':>25} {'v₂':>4} {'odd part':>25} {'prime?':>7} {'factors':>30}")
print("  " + "-" * 98)

for n in range(3, 21, 2):
    a = A000568[n]
    v = v2(a)
    op = a >> v

    # Factor the odd part
    if op < 10**15:
        f = factor(op)
        fc = Counter(f)
        fstr = " × ".join(f"{b}^{e}" if e > 1 else str(b) for b, e in sorted(fc.items()))
        if not fstr: fstr = "1"
        pr = "YES" if is_prime(op) and op > 1 else "no"
    else:
        fstr = "(too large)"
        pr = "?"

    op_str = str(op)
    if len(op_str) > 25:
        op_str = op_str[:10] + "..." + op_str[-10:]

    print(f"  {n:>3} {str(a):>25} {v:>4} {op_str:>25} {pr:>7} {fstr:>30}")

# ═══════════════════════════════════════════════════════════════════
# SECTION 8: CONNECTING v₂(a(n)) TO H(T_p) ORBIT STRUCTURE
# ═══════════════════════════════════════════════════════════════════
print("\n\n" + "═" * 78)
print("  SECTION 8: THE BRIDGE — v₂(a(p)) AND H(T_p)")
print("═" * 78)

print("""
  For Paley prime p ≡ 3 mod 4:

  FACT 1 (THM-305): v₂(a(p)) = (p-1)/2.
  So a(p) = 2^{(p-1)/2} × (odd number).

  FACT 2 (Rédei): H(T) is always odd for every tournament T.
  So v₂(H(T_p)) = 0.

  FACT 3 (THM-F3): p(p-1)/2 divides H(T_p).
  Since H(T_p) is odd, this means p(p-1)/2 must also be odd for the
  division to work. Check: p is odd, p-1 is even, so p(p-1)/2 = p × (p-1)/2.
  v₂(p(p-1)/2) = v₂(p-1) - 1 (since p is odd).
  For p ≡ 3 mod 4: p-1 ≡ 2 mod 4, so v₂(p-1) = 1.
  Therefore v₂(p(p-1)/2) = 1 - 1 = 0. ✓ p(p-1)/2 is ODD.

  For p ≡ 3 mod 8: p-1 ≡ 2 mod 8, v₂(p-1) = 1, v₂(|Aff|) = 0. ✓
  For p ≡ 7 mod 8: p-1 ≡ 6 mod 8, v₂(p-1) = 1, v₂(|Aff|) = 0. ✓

  So for ALL Paley primes: |Aff(QR)| = p(p-1)/2 is ODD.
  And H(T_p) is ODD. So H(T_p)/|Aff| is a ratio of odd numbers.

  FACT 4 (HYP-1713): H(T_p)/|Aff| is odd.
  This is automatically true since both numerator and denominator are odd!

  WAIT — this means HYP-1713 is TRIVIALLY TRUE for Paley primes!
  H(T_p) is odd (Rédei), |Aff| is odd (since p(p-1)/2 has v₂ = 0 for
  p ≡ 3 mod 4), and an odd number divided by an odd number is odd if
  the division is exact.

  Actually: odd/odd is always odd if it's an integer. Because if a = b·q
  with a, b odd, then q must be odd (since b·q odd ⟹ q odd).

  SO HYP-1713 (orbit parity: H/|Aff| is odd) IS PROVED for Paley primes!
""")

print("  VERIFICATION:")
for p in [3, 7, 11, 19, 23]:
    aff = p * (p - 1) // 2
    v2_aff = v2(aff)
    print(f"    p={p}: |Aff| = {aff}, v₂(|Aff|) = {v2_aff}, "
          f"p mod 4 = {p % 4}, p mod 8 = {p % 8}")

print(f"""
  ┌────────────────────────────────────────────────────────────────┐
  │  THEOREM (HYP-1713 PROVED):                                    │
  │                                                                │
  │  For Paley prime p ≡ 3 (mod 4): H(T_p)/|Aff(QR)| is ODD.     │
  │                                                                │
  │  Proof: H(T_p) is odd (Rédei's theorem).                      │
  │  |Aff(QR)| = p(p-1)/2. Since p ≡ 3 mod 4: v₂(p-1) = 1,      │
  │  so v₂(p(p-1)/2) = 0 + 1 - 1 = 0, i.e., |Aff| is odd.       │
  │  An integer = odd / odd is always odd.   □                     │
  │                                                                │
  │  COROLLARY: H(T_p)/p ≡ (p-1)/2 (mod p-1).                    │
  │  Proof: H/p = |Aff|/p × (H/|Aff|) = (p-1)/2 × (odd number). │
  │  So H/p ≡ (p-1)/2 × 1 = (p-1)/2 (mod p-1).                  │
  │  (Since (p-1)/2 × odd ≡ (p-1)/2 mod (p-1).)   □              │
  └────────────────────────────────────────────────────────────────┘

  NOTE: This also proves the STRONGER form of HYP-1713:
  H(T_p)/p ≡ (p-1)/2 (mod p-1).

  The proof is: H/p = ((p-1)/2) × k where k = H/|Aff| is odd.
  So H/p mod (p-1) = ((p-1)/2 × k) mod (p-1).
  Since k is odd: k = 2m+1, so (p-1)/2 × k = (p-1)/2 + m(p-1).
  Therefore H/p ≡ (p-1)/2 (mod p-1).   □
""")

# ═══════════════════════════════════════════════════════════════════
# SECTION 9: EXTENDING TO ALL PRIMES p ≡ 3 mod 4
# ═══════════════════════════════════════════════════════════════════
print("\n" + "═" * 78)
print("  SECTION 9: THE GENERAL CASE — ALL PRIMES?")
print("═" * 78)

print("""
  Does the orbit parity theorem extend beyond Paley primes?

  For ANY tournament T on n vertices with automorphism group Aut(T):
  - H(T) is odd (Rédei)
  - |Aut(T)| can be even or odd
  - If |Aut(T)| is odd, then H(T)/|Aut(T)| is automatically odd

  For Paley T_p: |Aut| = p(p-1)/2 which is odd for p ≡ 3 mod 4.
  This is because the full automorphism group is Aff(QR), and QR dilations
  have order (p-1)/2 which is odd (since (p-1)/2 = (4k+2)/2 = 2k+1).

  For OTHER tournaments: if |Aut(T)| is even (contains an element of
  order 2), then H(T)/|Aut(T)| divides an odd by an even number,
  which would make it a non-integer or have v₂ < 0. But orbit-stabilizer
  gives H(T)/|Aut(T)| = number of Aut-orbits on Ham paths, which is a
  positive integer. So if |Aut| is even: H/|Aut| must still be an integer,
  which means 2 | H... but H is odd! Contradiction: |Aut| CANNOT be even!

  WAIT — |Aut| CAN be even. The orbit-stabilizer theorem gives:
  |orbits| = (1/|Aut|) × Σ_{g ∈ Aut} |Fix(g)|

  This sum equals H only if the group acts on HAM PATHS. Actually:
  |Aut| × |orbits| = Σ |Fix(g)| ≠ H in general.

  The correct relation: Aut acts on the set of H paths. By Burnside:
  |orbits| = (1/|Aut|) Σ_{g} |{paths fixed by g}|.

  For g = id: contributes H.
  For g ≠ id: contributes |Fix(g)|.
  Total: |Aut| × |orbits| = H + Σ_{g≠id} |Fix(g)|.

  So |orbits| = (H + Σ_{g≠id} Fix(g)) / |Aut|.

  For Paley: we showed |Aut| is odd and Fix(g) = 0 for g ≠ id
  (at least for the translation part). So |orbits| = H/|Aut| = odd/odd = odd.

  For general tournaments with even |Aut|: |orbits| can still be an integer
  because the Fix(g) terms contribute. The total H + Σ Fix(g) is divisible
  by |Aut|, but individual terms may not be.
""")

# ═══════════════════════════════════════════════════════════════════
# SECTION 10: SUMMARY OF ALL RESULTS
# ═══════════════════════════════════════════════════════════════════
print("\n" + "═" * 78)
print("  SUMMARY OF RESULTS")
print("═" * 78)

print("""
  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  PROVED IN THIS SESSION:
  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  1. HYP-1713 IS PROVED for Paley primes p ≡ 3 mod 4:
     H(T_p)/|Aff(QR)| is odd. Trivially follows from:
     - H odd (Rédei) and |Aff| odd (since v₂(p(p-1)/2) = 0).
     - Corollary: H/p ≡ (p-1)/2 mod (p-1).

  2. HYP-1712 REDUCED to a concrete algebraic identity:
     v₂(a(n)) = n/2 + v₂(c(n)) for even n
     reduces to showing:
     (a) Only 2-part coprime odd partitions (a,b) with gcd(a,b)=1
         achieve the minimum v₂ level in the Burnside sum.
     (b) v₂(Σ 1/(a(n-a))) = v₂(c(n)) where sum over coprime pairs.
     Both (a) and (b) verified computationally for n = 4..20.
     (a) likely provable by the same gap inequality as THM-305.

  3. STRUCTURAL BRIDGE between v₂(a(n)) and H(T_p):
     The 2-adic valuation v₂(a(p)) = (p-1)/2 = |QR_p| (THM-305)
     controls both the tournament count AND the orbit structure:
     - v₂(|Aff(QR)|) = 0 for Paley primes (proved)
     - H(T_p) odd × orbit-count formula → orbit count is odd

  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  OPEN (requiring further work):
  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  4. Complete proof of HYP-1712 part (b):
     Need to prove v₂(Σ_{a ∈ C} 1/(a(n-a))) = v₂(|C|)
     where C = coprime odd residues < n/2.
     This may follow from a partial-fraction identity + Kummer's theorem.

  5. Odd part structure of a(n)/2^{(n-1)/2}:
     Values 1, 3, 57, 11971, 28242289 include primes.
     No clear pattern found yet.
""")

print("\nComputation complete.")
