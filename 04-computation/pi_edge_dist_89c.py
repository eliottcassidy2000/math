#!/usr/bin/env python3
"""
pi_edge_dist_89c.py — Common edge distribution analysis
opus-2026-03-14-S89c

The common edge distribution c(n,k) counts compatible permutation pairs
(π,σ) with exactly k common undirected edges, divided by n!.

Known data:
n=3: [2, 1]                          sum = 3 = A000255(1)
n=4: [2, 5, 3, 1]                    sum = 11 = A000255(2)  [wait, n=4 has 3 edges max]
Actually for n vertices, max common edges = n-1, and k ranges from 0 to n-1.

Let me recompute carefully. The c(n,k)/n! table is:
n=3, k=0..2: [0, 2, 1]  (k=0: 0 common edges, k=2: all edges common)
n=4, k=0..3: [2, 5, 3, 1]
n=5, k=0..4: [14, 20, 14, 4, 1]
n=6, k=0..5: [90, 115, 72, 26, 5, 1]
n=7, k=0..6: [646, 790, 467, 168, 41, 6, 1]

Goal: find recurrences or closed forms for these.
"""

from itertools import permutations
from math import factorial, comb
from fractions import Fraction

# Compute all distributions
all_dist = {}

for n in range(3, 8):
    if n <= 6:
        all_perms = list(permutations(range(n)))
        dist = {}
        for pi in all_perms:
            epi_d = set((pi[i], pi[i+1]) for i in range(n-1))
            epi_u = set((min(pi[i],pi[i+1]), max(pi[i],pi[i+1])) for i in range(n-1))
            for sig in all_perms:
                esig_d = set((sig[i], sig[i+1]) for i in range(n-1))
                if any((v,u) in esig_d for u,v in epi_d):
                    continue
                esig_u = set((min(sig[i],sig[i+1]), max(sig[i],sig[i+1])) for i in range(n-1))
                k = len(epi_u & esig_u)
                dist[k] = dist.get(k, 0) + 1
        # Normalize by n!
        cnts = [dist.get(k, 0) // factorial(n) for k in range(n)]
        all_dist[n] = cnts
    else:
        # n=7 from computation
        raw = {0: 3255840, 1: 3981600, 2: 2353680, 3: 846720, 4: 206640, 5: 30240, 6: 5040}
        cnts = [raw[k] // factorial(7) for k in range(7)]
        all_dist[n] = cnts

# Display
print("Common edge distribution c(n,k)/n!:")
print(f"{'n':>3}", end="")
for k in range(7):
    print(f"  {'k='+str(k):>8}", end="")
print(f"  {'Σ':>8}")

for n in sorted(all_dist.keys()):
    print(f"{n:>3}", end="")
    for k in range(7):
        if k < n:
            print(f"  {all_dist[n][k]:>8}", end="")
        else:
            print(f"  {'·':>8}", end="")
    print(f"  {sum(all_dist[n]):>8}")

# k = n-1 column (max overlap): always 1
print("\n=== Column analysis ===")

# k = n-1: 1, 1, 1, 1, 1 (always 1)
print("k = n-1: always 1 (only the identity pair)")

# k = n-2: 2, 3, 4, 5, 6 = n-1
print(f"k = n-2: {[all_dist[n][n-2] for n in range(3,8)]} = n-1")
print("  Why: one edge differs. Choose which of n-1 edges to drop (n-1 choices),")
print("  and the replacement must be compatible. Exactly 1 valid replacement each time.")

# k = n-3: -, 5, 14, 26, 41 (starting at n=4)
col_nm3 = [all_dist[n][n-3] for n in range(4, 8)]
print(f"k = n-3: {col_nm3}")
# Differences: 9, 12, 15 → second diff 3, 3
print(f"  Differences: {[col_nm3[i+1]-col_nm3[i] for i in range(len(col_nm3)-1)]}")
print(f"  Second differences: {[col_nm3[i+2]-2*col_nm3[i+1]+col_nm3[i] for i in range(len(col_nm3)-2)]}")
# So quadratic: fit
# 5, 14, 26, 41 at n = 4, 5, 6, 7
# a(n) = An² + Bn + C
# 16A + 4B + C = 5
# 25A + 5B + C = 14
# → 9A + B = 9
# 36A + 6B + C = 26
# → 11A + B = 12
# → 2A = 3, A = 3/2
# B = 9 - 27/2 = -9/2
# C = 5 - 24 + 18 = -1
print(f"  Fit: c(n, n-3)/n! = (3n² - 9n - 2)/2")
for n in range(4, 8):
    predicted = (3*n**2 - 9*n - 2) // 2
    print(f"    n={n}: predicted={predicted}, actual={all_dist[n][n-3]}, match={'✓' if predicted == all_dist[n][n-3] else '✗'}")

# k = n-4: -, -, 14, 72, 168 (starting at n=5)
if 5 in all_dist:
    col_nm4 = [all_dist[n][n-4] for n in range(5, 8)]
    print(f"\nk = n-4: {col_nm4}")
    diffs = [col_nm4[i+1]-col_nm4[i] for i in range(len(col_nm4)-1)]
    print(f"  Differences: {diffs}")

# k = 0 column: 0, 2, 14, 90, 646
col_0 = [all_dist[n][0] for n in range(3, 8)]
print(f"\nk = 0 (no common edges): {col_0}")
# Ratios
for i in range(1, len(col_0)):
    if col_0[i-1] > 0:
        print(f"  c({i+3},0)/c({i+2},0) = {col_0[i]/col_0[i-1]:.4f}")

print()
print("=" * 60)
print("E[H²] FROM COMMON EDGE DISTRIBUTION")
print("=" * 60)

# E[H²] = Σ_k c(n,k) × n! × 2^{k - 2(n-1)}
# = n! × 2^{-2(n-1)} × Σ_k c(n,k) × 2^k

for n in sorted(all_dist.keys()):
    dist = all_dist[n]
    # E[H²] = Σ_k (count_k) × 2^{-(2(n-1)-k)}
    # count_k = c(n,k) × n!
    eh2 = Fraction(0)
    for k in range(n):
        constrained = 2*(n-1) - k
        count_k = dist[k] * factorial(n)
        eh2 += Fraction(count_k, 2**constrained)

    eh = Fraction(factorial(n), 2**(n-1))
    cv2 = eh2 / eh**2 - 1

    # Alternative: E[H²] = 2^{-2(n-1)} × n! × Σ c(n,k) 2^k
    weighted = sum(dist[k] * 2**k for k in range(n))
    eh2_alt = Fraction(factorial(n) * weighted, 4**(n-1))

    print(f"\n  n={n}: E[H²] = {eh2} = {float(eh2):.6f}")
    print(f"    CV² = {cv2} = {float(cv2):.10f}")
    print(f"    Σ c(n,k) 2^k = {weighted}")
    print(f"    E[H²] alt = {eh2_alt} = {float(eh2_alt):.6f}")
    print(f"    Match: {'✓' if eh2 == eh2_alt else '✗'}")

# The weighted sum Σ c(n,k) 2^k
print("\n" + "=" * 60)
print("WEIGHTED SUM W(n) = Σ c(n,k) 2^k")
print("=" * 60)

for n in sorted(all_dist.keys()):
    dist = all_dist[n]
    W = sum(dist[k] * 2**k for k in range(n))
    # E[H²] = n! × W / 4^{n-1}
    # CV² = E[H²]/E[H]² - 1 = n!×W/4^{n-1} / (n!/2^{n-1})² - 1
    # = W × 4^{n-1} / (n! × 4^{n-1}) - 1 ... wait
    # E[H]² = (n!/2^{n-1})² = n!²/4^{n-1}
    # E[H²]/E[H]² = W×n!/(4^{n-1}) / (n!²/4^{n-1}) = W/n!
    # So CV² = W/n! - 1 !!!

    cv2 = Fraction(W, factorial(n)) - 1
    print(f"  n={n}: W = {W}, W/n! = {Fraction(W, factorial(n))}, CV² = {cv2} = {float(cv2):.10f}")

print()
print("KEY INSIGHT: CV² = W(n)/n! - 1 where W(n) = Σ c(n,k) 2^k")
print("And E[H²]/E[H]² = W(n)/n!")
print()

# So we need: W(n) = Σ_{k=0}^{n-1} c(n,k) × 2^k
# W values: let me compute
W_vals = {}
for n in sorted(all_dist.keys()):
    dist = all_dist[n]
    W = sum(dist[k] * 2**k for k in range(n))
    W_vals[n] = W
    print(f"  W({n}) = {W}")

# W/n!:
print()
for n in sorted(W_vals.keys()):
    r = Fraction(W_vals[n], factorial(n))
    print(f"  W({n})/n! = {r}")

# W sequence: 4, 16, 79, 522, 5080
# Check OEIS?
print(f"\nW sequence: {[W_vals[n] for n in sorted(W_vals.keys())]}")

# Also: Σ c(n,k) = A000255(n-1) and Σ c(n,k) 2^k = W(n)
# So W(n) is a "2-weighted" version of A000255

# Check: does W(n) satisfy a recurrence?
# W(3) = 4
# W(4) = 16
# W(5) = 79
# W(6) = 522
# W(7) = 5080

# Try: W(n) = a*n*W(n-1) + b*(n-1)*W(n-2)?
# 79 = 5a*16 + 4b*4 = 80a + 16b
# 522 = 6a*79 + 5b*16 = 474a + 80b
# From first: 5a + b = 79/16... not integer, so this recurrence doesn't work directly.

# Try: W(n) = a*W(n-1) + b*W(n-2)?
# 79 = 16a + 4b
# 522 = 79a + 16b
# From first: b = (79 - 16a)/4
# Sub: 522 = 79a + 4(79 - 16a) = 79a + 316 - 64a = 15a + 316
# 15a = 206, a = 206/15... not integer.

# Try: W(n) = (n+1)*W(n-1) - (n-2)*W(n-2)?
# n=5: 6*16 - 3*4 = 96 - 12 = 84 ≠ 79. No.

# Try: W(n) = n*W(n-1) + something?
# 16 = 4*4 + 0
# 79 = 5*16 - 1 = 79 ✓
# 522 = 6*79 + 48 = 474 + 48 = 522 ✓
# 5080 = 7*522 + 1426 = 3654 + 1426 = 5080 ✓
# Remainders: 0, -1, 48, 1426
# Hmm, those aren't clean. Let me try W(n) - n*W(n-1):
print("\nW(n) - n*W(n-1):")
for n in range(4, 8):
    diff = W_vals[n] - n * W_vals[n-1]
    print(f"  n={n}: W({n}) - {n}*W({n-1}) = {W_vals[n]} - {n*W_vals[n-1]} = {diff}")

# W(n) - n*W(n-1) - (n-1)*W(n-2)?
print("\nW(n) - n*W(n-1) - (n-1)*W(n-2):")
for n in range(5, 8):
    diff = W_vals[n] - n * W_vals[n-1] - (n-1) * W_vals[n-2]
    print(f"  n={n}: {diff}")

# Compare with A000255 recurrence: a(n) = n*a(n-1) + (n-1)*a(n-2)
# Our W seems to NOT follow the same recurrence.

# Let me try: W(n) = c*n*W(n-1) + d*(n-1)*W(n-2) for constants c, d
# n=5: 79 = 5c*16 + 4d*4 = 80c + 16d
# n=6: 522 = 6c*79 + 5d*16 = 474c + 80d
# From first: 5c + d = 79/16 — not integer. So no simple (c,d) recurrence.

# Let me try the EGF approach.
# A000255 has EGF exp(-x)/(1-x)^2
# W(n)/n! = Σ c(n,k)/n! × 2^k where c(n,k)/n! is the number of compatible
# σ with identity that have exactly k common edges with identity.
# So W(n)/n! = E_σ[2^{common(id, σ)} | compatible]... no, it's a sum not expectation.
# Actually W(n) = Σ_σ 2^{common(id, σ)} where σ ranges over compatible perms.

print("\n" + "=" * 60)
print("INTERPRETATION: W(n) = Σ_{σ compat with id} 2^{common(id,σ)}")
print("=" * 60)
print("W(n)/n! = E[H²]/E[H]² = 1 + CV²")
print()

# Let's try the EGF of W(n):
# W(n)/n! values: 4/6, 16/24, 79/120, 522/720, 5080/5040
# = 2/3, 2/3, 79/120, 29/40, 635/504 (? let me check)
for n in sorted(W_vals.keys()):
    r = Fraction(W_vals[n], factorial(n))
    print(f"  W({n})/n! = {r} = {float(r):.10f}")

# These are exactly E[H²]/E[H]² = 4/3, 4/3, 79/60, 58/45, 635/504
# Wait: W(3)/3! = 4/6 = 2/3, but E[H²]/E[H]² = 4/3 at n=3.
# Discrepancy! Let me recheck.

# E[H²] = n! × W / 4^{n-1}
# E[H]² = (n!/2^{n-1})² = n!² / 4^{n-1}
# E[H²]/E[H]² = W / n!
# At n=3: W = 4, n! = 6, so E[H²]/E[H]² = 4/6 = 2/3.
# But we KNOW E[H²]/E[H]² = 4/3 at n=3!
# Something is wrong.

# Let me recheck: E[H²] = Σ_{compat (π,σ)} 2^{-constrained}
# where constrained = |E(π) ∪ E(σ)| = 2(n-1) - common
# So E[H²] = Σ_{compat} 2^{-(2(n-1) - common)} = 2^{-2(n-1)} Σ_{compat} 2^{common}
# And Σ_{compat} 2^{common} = Σ_π Σ_{σ compat with π} 2^{common(π,σ)}
# = n! × W(n)  [by symmetry, each π contributes W(n)]
# So E[H²] = n! × W(n) / 4^{n-1}

# At n=3: E[H²] = 6 × 4 / 4² = 24/16 = 3/2
# But E[H²] should be 3 (from pi_harmonic)!
# Check: E[H] = 3/2, E[H²]/E[H]² = 4/3, so E[H²] = 4/3 × (3/2)² = 4/3 × 9/4 = 3.
# But our formula gives 3/2. Off by factor 2!

# The issue: at n=3, we got c(3, k=0) = 0, c(3, k=1) = 2, c(3, k=2) = 1.
# But the TOTAL compat pairs = n! × Σ c(n,k) = 6 × (0+2+1) = 18.
# E[H²] = (1/4^2) × Σ_{all compat (π,σ)} 2^{common(π,σ)}
# = (1/16) × n! × W(n) ... wait, let me redo.

# Actually c(n,k) is defined as count(k) / n! where count(k) is the number
# of compatible PAIRS. So the total sum is Σ_k count(k) 2^k.
# n=3: count = {0: 0, 1: 12, 2: 6}
# Σ count(k) 2^k = 0 + 12×2 + 6×4 = 24 + 24 = 48??? Let me recheck.
# Wait, c(n,k)/n! was {0, 2, 1} which means count(k) = {0, 12, 6}.
# Total pairs = 0+12+6 = 18 ✓
# Σ count(k) 2^k = 0 + 12*2 + 6*4 = 0 + 24 + 24 = 48.
# But n! × W(n) = 6 × 4 = 24. That's half!

# Hmm, my W(n) = Σ c(n,k) 2^k = 0 + 2×2 + 1×4 = 0+4+4 = 8???
# Wait no, W(3) = sum(dist[k]*2^k for k in range(3)) where dist = [0, 2, 1]
# = 0*1 + 2*2 + 1*4 = 0 + 4 + 4 = 8. But I got W(3) = 4 above!

# Let me re-run and check.
print("\n=== RECHECK W VALUES ===")
for n in sorted(all_dist.keys()):
    dist = all_dist[n]
    W = sum(dist[k] * 2**k for k in range(n))
    print(f"  n={n}: dist={dist}, W = Σ c(n,k)×2^k = {' + '.join(f'{dist[k]}×{2**k}' for k in range(n))} = {W}")

print("\nDone!")
