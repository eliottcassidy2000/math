#!/usr/bin/env python3
"""
wn_deep_investigation_s116.py — Deep investigation of W(n) as a sequence
opus-2026-03-21-S116

W(n) = Σ_T H(T) where sum is over all labeled tournaments on n vertices
     = Σ_{σ ∈ NUD(n)} 2^{adj1(σ)}
     = 2^{C(n,2)} × E[H]   ... wait, that's wrong.

Actually:
  W(n) = Σ_T H(T) = sum of ALL Hamiltonian path counts across ALL tournaments
  E[H] = W(n) / 2^{C(n,2)} = n! / 2^{n-1}

So W(n) = n! × 2^{C(n,2)} / 2^{n-1} = n! × 2^{C(n,2)-(n-1)} = n! × 2^{(n-1)(n-2)/2}

Wait, let me check: C(n,2) - (n-1) = n(n-1)/2 - (n-1) = (n-1)(n-2)/2.

So W(n) = n! × 2^{(n-1)(n-2)/2}?

Let me verify:
  n=3: W(3)=8. n!×2^{1} = 6×2 = 12. NO!
  n=3: From the data, W(3) = 8. And 2^{C(3,2)} = 8 tournaments.
  E[H] = 3/2. So Σ_T H(T) = 8 × 3/2 = 12. But W(3) = 8?

Hmm, maybe W(n) is NOT the sum of H(T) over all tournaments.
Let me re-read the definition.

From hook-C: W(n) = Σ_{σ ∈ NUD(n)} 2^{adj1(σ)}
where NUD(n) = permutations with no unit descent (= A000255(n-1))
and adj1(σ) = #{unit ascents, i.e., positions where σ(i+1)=σ(i)+1}

So W(n) is a WEIGHTED SUM over NUD permutations, NOT the sum of H(T).

The connection: W(n)/n! = 1 + CV²(H) where CV² = Var(H)/E[H]²

So: W(n) = n! × (1 + Var(H)/E[H]²) = n! + n! × Var(H)/E[H]²
         = n! + (2^{n-1})² × Var(H) / (n!)

Wait: E[H]² = (n!/2^{n-1})², and Var(H) × n! = n! × Var(H).
W(n) = n!(1 + Var(H)/E[H]²) = n! + n! × Var(H) × 4^{n-1} / (n!)²
     = n! + Var(H) × 4^{n-1} / n!

Hmm. Let me just verify numerically.
"""

from math import comb, factorial, log, sqrt, exp
from fractions import Fraction
from collections import Counter
import time

print("=" * 72)
print("  DEEP INVESTIGATION OF W(n)")
print("  opus-2026-03-21-S116")
print("=" * 72)
print()

# The known values
W = {
    1: 1, 2: 2, 3: 8, 4: 32, 5: 158, 6: 928, 7: 6350, 8: 49752,
    9: 439670, 10: 4327904, 11: 46963358, 12: 556953448,
    13: 7166360054, 14: 99428495088, 15: 1479600188798,
}

# =============================================================================
# PART 1: BASIC PROPERTIES
# =============================================================================

print("=" * 72)
print("  PART 1: BASIC PROPERTIES OF W(n)")
print("=" * 72)
print()

print("  W(n) = Σ_{σ ∈ NUD(n)} 2^{adj1(σ)}")
print("  NUD(n) = permutations of [n] with no unit descent (σ(i+1) ≠ σ(i)-1)")
print("  adj1(σ) = #{unit ascents: positions where σ(i+1) = σ(i)+1}")
print()

def A000255(n):
    if n <= 1: return 1
    a, b = 1, 1
    for k in range(2, n+1):
        a, b = b, k*b + (k-1)*a
    return b

print(f"  {'n':>3} | {'W(n)':>18} | {'n!':>12} | {'NUD(n)':>10} | {'W/n!':>12} | {'W/NUD':>10} | {'W/W(n-1)':>10}")
print("  " + "-" * 85)

for n in range(1, 16):
    nf = factorial(n)
    nud = A000255(n)
    ratio_nf = W[n] / nf
    ratio_nud = W[n] / nud
    ratio_prev = W[n] / W[n-1] if n > 1 else 0

    print(f"  {n:3d} | {W[n]:18d} | {nf:12d} | {nud:10d} | {ratio_nf:12.6f} | {ratio_nud:10.6f} | {ratio_prev:10.4f}")

print()

# =============================================================================
# PART 2: W(n)/n! = 1 + CV²(H) CONNECTION
# =============================================================================

print("=" * 72)
print("  PART 2: W(n)/n! = 1 + CV²(H)")
print("=" * 72)
print()

print("  CV²(H) = Var(H) / E[H]²")
print("  E[H] = n! / 2^{n-1}")
print()

# Known exact Var(H) values
Var_H = {
    3: Fraction(3, 4),
    4: Fraction(3),
    5: Fraction(285, 16),
    6: Fraction(585, 4),
    7: Fraction(206325, 128),
}

print(f"  {'n':>3} | {'W(n)/n!':>14} | {'1+CV²':>14} | {'CV²':>12} | {'Var(H)/E[H]²':>14}")
print("  " + "-" * 65)

for n in range(3, 8):
    W_over_nf = Fraction(W[n], factorial(n))
    E_H = Fraction(factorial(n), 2**(n-1))
    if n in Var_H:
        CV2 = Var_H[n] / (E_H * E_H)
        one_plus_CV2 = 1 + CV2
        print(f"  {n:3d} | {float(W_over_nf):14.8f} | {float(one_plus_CV2):14.8f} | {float(CV2):12.8f} | {float(Var_H[n]):14.6f}/{float(E_H**2):.4f}")
        match = W_over_nf == one_plus_CV2
        if not match:
            print(f"       MISMATCH: W/n! = {W_over_nf}, 1+CV² = {one_plus_CV2}")

print()

# =============================================================================
# PART 3: RATIOS AND GROWTH RATE
# =============================================================================

print("=" * 72)
print("  PART 3: GROWTH RATE OF W(n)")
print("=" * 72)
print()

print("  W(n)/W(n-1) should approach n (since W(n) ~ n! × constant)")
print()
print(f"  {'n':>3} | {'W(n)/W(n-1)':>12} | {'n':>4} | {'ratio/n':>10} | {'ratio-n':>10}")
print("  " + "-" * 50)

for n in range(2, 16):
    r = W[n] / W[n-1]
    print(f"  {n:3d} | {r:12.6f} | {n:4d} | {r/n:10.6f} | {r-n:10.4f}")

print()
print("  W(n)/W(n-1) ≈ n - 1 + 2/n + ...")
print()

# More precise: W(n)/W(n-1) - n
print(f"  {'n':>3} | {'W(n)/W(n-1)-n':>14} | {'1/n':>8} | {'(W/W-n)*n':>12}")
print("  " + "-" * 45)
for n in range(3, 16):
    r = W[n] / W[n-1]
    diff = r - n
    print(f"  {n:3d} | {diff:14.6f} | {1/n:8.6f} | {diff*n:12.6f}")

print()

# =============================================================================
# PART 4: FACTORIZATION OF W(n)
# =============================================================================

print("=" * 72)
print("  PART 4: PRIME FACTORIZATION OF W(n)")
print("=" * 72)
print()

def factorize(n):
    if n <= 1: return {}
    factors = {}
    d = 2
    m = abs(n)
    while d * d <= m:
        while m % d == 0:
            factors[d] = factors.get(d, 0) + 1
            m //= d
        d += 1
    if m > 1:
        factors[m] = factors.get(m, 0) + 1
    return factors

for n in range(3, 13):
    f = factorize(W[n])
    primes = sorted(f.keys())
    print(f"  W({n:2d}) = {W[n]:>15d} = ", end="")
    parts = []
    for p in primes:
        if f[p] == 1:
            parts.append(str(p))
        else:
            parts.append(f"{p}^{f[p]}")
    print(" × ".join(parts))

print()

# =============================================================================
# PART 5: W(n) mod small primes
# =============================================================================

print("=" * 72)
print("  PART 5: W(n) mod SMALL PRIMES")
print("=" * 72)
print()

for p in [2, 3, 5, 7, 11, 13]:
    mods = [W[n] % p for n in range(1, 16)]
    print(f"  W(n) mod {p:2d}: {mods}")

print()

# Check W(n) mod n
print("  W(n) mod n:")
for n in range(2, 16):
    print(f"    W({n:2d}) mod {n:2d} = {W[n] % n}")

print()

# =============================================================================
# PART 6: THE N(n,j) DISTRIBUTION — UNIT ASCENT COUNTS
# =============================================================================

print("=" * 72)
print("  PART 6: N(n,j) — THE UNIT ASCENT DISTRIBUTION AMONG NUD PERMS")
print("=" * 72)
print()

Nnj = {
    1: [1],
    2: [0, 1],
    3: [0, 2, 1],
    4: [2, 5, 3, 1],
    5: [14, 20, 14, 4, 1],
    6: [90, 115, 72, 26, 5, 1],
    7: [646, 790, 467, 168, 41, 6, 1],
    8: [5242, 6217, 3557, 1285, 319, 59, 7, 1],
}

print("  N(n,j) = #{NUD permutations of [n] with exactly j unit ascents}")
print("  W(n) = Σ_j N(n,j) × 2^j")
print()

for n in range(1, 9):
    coeffs = Nnj[n]
    nud_total = sum(coeffs)
    w_check = sum(c * 2**j for j, c in enumerate(coeffs))
    print(f"  n={n}: N(n,j) = {coeffs}")
    print(f"    NUD({n}) = {nud_total} = A000255({n}) = {A000255(n)}: {nud_total == A000255(n)}")
    print(f"    W({n}) = Σ N(n,j)×2^j = {w_check}: {w_check == W[n]}")
    print()

# =============================================================================
# PART 7: THE j=0 COLUMN — HERTZSPRUNG NUMBERS
# =============================================================================

print("=" * 72)
print("  PART 7: N(n,0) — HERTZSPRUNG NUMBERS")
print("=" * 72)
print()

print("  N(n,0) = #{NUD perms with NO unit ascents AND no unit descents}")
print("  = #{perms with no succession of either type}")
print("  = Hertzsprung numbers (OEIS A000179)")
print()

hertz = [Nnj[n][0] for n in range(1, 9)]
print(f"  Hertzsprung sequence: {hertz}")
print()

# These are sequence A000179: 1, 0, 0, 2, 14, 90, 646, 5242, ...
# Which counts permutations with no 3-letter arithmetic progression
# (no i,i+1 with σ(i+1)-σ(i) = ±1)

print("  A000179: 1, 0, 0, 2, 14, 90, 646, 5242, 47622, ...")
print(f"  Our N(n,0): {hertz}")
print(f"  Match with A000179: {hertz == [1, 0, 0, 2, 14, 90, 646, 5242]}")
print()

# =============================================================================
# PART 8: THE GENERATING FUNCTION APPROACH
# =============================================================================

print("=" * 72)
print("  PART 8: GENERATING FUNCTION STRUCTURE")
print("=" * 72)
print()

print("""
  From hook-C: the bivariate GF of the x-tribonacci is:
    Σ F(N,x) y^N = (1 + 2xy + xy²) / (1 - y - xy² - xy³)

  At x=1: denominator 1 - y - y² - y³ has tribonacci root τ ≈ 1.8393.
  At x=2: denominator 1 - y - 2y² - 2y³.

  The W(n) sequence is related to F(n, 2) via the NUD transfer matrix.

  Let me compute the characteristic polynomial at x=2:
    p(y) = 1 - y - 2y² - 2y³ = 0
    2y³ + 2y² + y - 1 = 0
""")

import numpy as np

# Roots of 2y³ + 2y² + y - 1 = 0
coeffs_poly = [2, 2, 1, -1]  # 2y³ + 2y² + y - 1
roots = np.roots(coeffs_poly)
print(f"  Roots of 2y³ + 2y² + y - 1 = 0:")
for i, r in enumerate(roots):
    print(f"    r{i+1} = {r:.10f}")

# The dominant root controls the growth
dominant = max(roots, key=lambda r: abs(r))
print(f"\n  Dominant root: {dominant:.10f}")
print(f"  |dominant| = {abs(dominant):.10f}")
print()

# Check: does W(n) / (dominant^n × n!) approach a constant?
print("  W(n) / (dominant^n × n!):")
for n in range(3, 16):
    ratio = W[n] / (abs(dominant)**n * factorial(n))
    print(f"    n={n:2d}: {ratio:.8f}")
print()

# =============================================================================
# PART 9: W(n) AND THE SECOND MOMENT E[H²]
# =============================================================================

print("=" * 72)
print("  PART 9: W(n) AND E[H²]")
print("=" * 72)
print()

print("""
  W(n)/n! = 1 + CV²(H) = 1 + Var(H)/E[H]²
  = E[H²]/E[H]²

  So: W(n)/n! = E[H²] / E[H]²
  And: E[H²] = W(n) × E[H]² / n! = W(n) × (n!)² / (4^{n-1} × n!)
             = W(n) × n! / 4^{n-1}

  Wait: E[H] = n!/2^{n-1}, so E[H]² = (n!)² / 4^{n-1}.
  W(n)/n! = E[H²]/E[H]², so E[H²] = W(n) × E[H]² / n!
  = W(n) × (n!)² / (4^{n-1} × n!) = W(n) × n! / 4^{n-1}.

  Let me verify:
""")

for n in range(3, 8):
    E_H = Fraction(factorial(n), 2**(n-1))
    E_H2_from_W = Fraction(W[n] * factorial(n), 4**(n-1))

    # Known E[H²]
    E_H2_known = {
        3: Fraction(3),
        4: Fraction(12),
        5: Fraction(1185, 16),
        6: Fraction(1305, 2),
        7: Fraction(206325, 128) + Fraction(factorial(7), 2**(7-1))**2,
    }
    # Actually E[H²] = Var(H) + E[H]²
    if n in Var_H:
        E_H2 = Var_H[n] + E_H**2
        print(f"  n={n}: E[H²] from Var+E² = {E_H2} = {float(E_H2):.4f}")
        print(f"         E[H²] from W(n)   = {E_H2_from_W} = {float(E_H2_from_W):.4f}")
        print(f"         Match: {E_H2 == E_H2_from_W}")
        print()

# =============================================================================
# PART 10: OEIS SEARCHES AND CONNECTIONS
# =============================================================================

print("=" * 72)
print("  PART 10: CONNECTIONS AND PATTERNS")
print("=" * 72)
print()

# W(n) mod 2
print("  W(n) mod 2: ", [W[n] % 2 for n in range(1, 16)])
print("  Parity: ", ['E' if W[n] % 2 == 0 else 'O' for n in range(1, 16)])
print()

# W(n)/2 for even W(n)
print("  W(n)/2 when W(n) is even:")
for n in range(1, 13):
    if W[n] % 2 == 0:
        print(f"    W({n})/2 = {W[n]//2}")
    else:
        print(f"    W({n}) is odd")
print()

# The sequence W(n)/n!
print("  W(n)/n! (= 1 + CV²):")
for n in range(1, 13):
    r = Fraction(W[n], factorial(n))
    print(f"    n={n}: {r} = {float(r):.10f}")
print()

# Second differences of W(n)/n!
print("  Differences of W(n)/n!:")
ratios = [Fraction(W[n], factorial(n)) for n in range(1, 13)]
first_diff = [ratios[i+1] - ratios[i] for i in range(len(ratios)-1)]
print(f"    First differences: {[float(d) for d in first_diff[:8]]}")
second_diff = [first_diff[i+1] - first_diff[i] for i in range(len(first_diff)-1)]
print(f"    Second differences: {[float(d) for d in second_diff[:7]]}")
print()

# =============================================================================
# PART 11: W(n) AND THE OVERLAP POLYNOMIAL F(2)
# =============================================================================

print("=" * 72)
print("  PART 11: W(n) AND THE OVERLAP POLYNOMIAL")
print("=" * 72)
print()

print("""
  From S106: E[H²] = F_n(2) / 4^{n-1} where F_n(x) = conflict-free overlap GF.
  And: E[H²] = W(n) × n! / 4^{n-1}.
  So: F_n(2) = W(n) × n!.

  THE OVERLAP POLYNOMIAL EVALUATED AT 2 IS W(n) × n!.

  This is the KEY CONNECTION: W(n) encodes the overlap structure at x=2.

  F_n(2) = n! × W(n):
""")

# Verify with known F_n(2) values from S106
# F_3(2) = 66? No, from the corrected computation: F_3(2) with conflict exclusion
# gave E[H²] = 3 at n=3. So F_3(2) = 3 × 4^2 = 48? Wait, 4^{n-1} = 4^2 = 16 at n=3.
# E[H²] = F_3(2)/16 = 3. So F_3(2) = 48.
# n!×W(n) = 6×8 = 48. YES!

for n in range(3, 8):
    fn2 = factorial(n) * W[n]
    E_H = Fraction(factorial(n), 2**(n-1))
    if n in Var_H:
        E_H2 = Var_H[n] + E_H**2
        fn2_from_EH2 = E_H2 * 4**(n-1)
        print(f"  n={n}: n!×W(n) = {fn2}")
        print(f"         F_n(2) from E[H²] = {fn2_from_EH2}")
        print(f"         Match: {Fraction(fn2) == fn2_from_EH2}")
        print()

# =============================================================================
# PART 12: DEEP CONNECTIONS
# =============================================================================

print("=" * 72)
print("  PART 12: DEEP CONNECTIONS")
print("=" * 72)
print()

print("""
  W(n) sits at the INTERSECTION of three structures:

  1. PERMUTATION COMBINATORICS:
     W(n) = Σ_{σ ∈ NUD(n)} 2^{adj1(σ)}
     It weights NUD permutations by 2^{unit ascents}.
     The NUD condition = no unit descent = succession-free (up to direction).

  2. TOURNAMENT STATISTICS:
     W(n)/n! = 1 + CV²(H) = E[H²]/E[H]²
     It encodes the second moment of Hamiltonian path counts.

  3. OVERLAP THEORY:
     n! × W(n) = F_n(2) = the conflict-free overlap GF at x=2
     It counts compatible permutation pairs weighted by 2^{shared arcs}.

  THESE THREE ARE THE SAME QUANTITY seen from different angles:
    Permutation → Tournament → Overlap

  The GENERATING FUNCTION has tribonacci structure at x=1
  and a modified tribonacci at x=2.

  THE OCR CONNECTION:
    OCR = 1 - Var(H|scores)/Var(H)
    Var(H) = E[H²] - E[H]² = (W(n)/n! - 1) × E[H]²
    So Var(H) = (W(n) - n!) × E[H]² / n! = (W(n) - n!) × (n!)² / (4^{n-1} × n!)
              = (W(n) - n!) × n! / 4^{n-1}

  And: (W(n) - n!) / 4^{n-1} = Var(H) / n!

  The RESIDUAL (W(n) - n!) encodes the variance!
  W(n) - n!: {W[n] - factorial(n) for n in range(3, 8)}
""")

for n in range(3, 10):
    residual = W[n] - factorial(n)
    print(f"  n={n}: W(n) - n! = {residual}")
    f = factorize(residual) if residual > 0 else {}
    print(f"    Factors: {f}")
print()

# The ratio (W(n)-n!)/n!
print("  (W(n)-n!)/n! = CV²(H):")
for n in range(3, 13):
    cv2 = Fraction(W[n] - factorial(n), factorial(n))
    print(f"    n={n}: CV² = {cv2} = {float(cv2):.8f}")
print()

# CV² should approach 2/n as n → ∞
print("  CV² × n (should approach 2):")
for n in range(3, 16):
    cv2 = (W[n] - factorial(n)) / factorial(n)
    print(f"    n={n}: CV²×n = {cv2*n:.6f}")
print()

print("=" * 72)
print("  SUMMARY")
print("=" * 72)
print()

print("""
  W(n) IS THREE THINGS AT ONCE:

  1. A weighted sum over NUD permutations: W = Σ 2^{adj1}
  2. The second moment normalizer: W/n! = E[H²]/E[H]²
  3. The overlap polynomial at 2: F(2) = n! × W

  The KEY FORMULA linking everything:
    Var(H) = (W(n) - n!) × n! / 4^{n-1}

  And the OCR:
    1 - OCR = Var(H|scores) / Var(H)
    = Var(H|scores) × 4^{n-1} / ((W(n) - n!) × n!)

  So: 1 - OCR = 4^{n-1} × Var(H|scores) / ((W(n) - n!) × n!)

  The OCR is COMPLETELY DETERMINED by W(n) and Var(H|scores).
  W(n) gives the total variance. Var(H|scores) gives the residual.
  Their ratio is the OCR.
""")

print("opus-2026-03-21-S116: DEEP INVESTIGATION OF W(n)")
