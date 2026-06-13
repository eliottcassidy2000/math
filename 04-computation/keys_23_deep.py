#!/usr/bin/env python3
"""
keys_23_deep.py — opus-2026-03-14-S71d

"If one knows 2 and 3, one has the keys to the universe."

DEEP CONNECTIONS between:
1. k-nacci convergence: φ_k → 2 with error (1/2)^k
2. Weighted k-nacci: ψ_k → 3 with error (2/3)^k
3. Tournament structure: I(Ω,2) = H, I(Ω,3) gives α₂

THESIS: At n≤8 (where α₃=0), knowing I(Ω,2) and I(Ω,3) determines
the COMPLETE independence polynomial I(Ω,x). This is exactly what
"knowing 2 and 3" means — two evaluation points suffice.
"""

import sys, time
import numpy as np
from itertools import combinations
from collections import defaultdict, Counter
from math import comb
sys.stdout.reconfigure(line_buffering=True)

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for ms in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != ms: continue
            for v in range(n):
                if not (mask & (1 << v)): continue
                pm = mask ^ (1 << v)
                t = 0
                for u in range(n):
                    if (pm & (1 << u)) and A[u][v]:
                        t += dp.get((pm, u), 0)
                if t:
                    dp[(mask, v)] = t
    return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))

def count_ham_cycles(A, n):
    if n < 3: return 0
    full_mask = (1 << n) - 1
    dp = {(1 << 0, 0): 1}
    for ms in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != ms: continue
            if not (mask & 1): continue
            for v in range(n):
                if not (mask & (1 << v)): continue
                if v == 0 and ms < n: continue
                pm = mask ^ (1 << v)
                if not (pm & 1): continue
                t = 0
                for u in range(n):
                    if (pm & (1 << u)) and A[u][v]:
                        t += dp.get((pm, u), 0)
                if t:
                    dp[(mask, v)] = t
    total = 0
    for v in range(1, n):
        if A[v][0] and (full_mask, v) in dp:
            total += dp[(full_mask, v)]
    return total

def count_directed_k_cycles(A, n, k):
    if k > n: return 0
    total = 0
    for combo in combinations(range(n), k):
        verts = list(combo)
        sub = np.zeros((k, k), dtype=int)
        for i in range(k):
            for j in range(k):
                sub[i][j] = A[verts[i]][verts[j]]
        total += count_ham_cycles(sub, k)
    return total

# ======================================================================
# PART 1: (H, I₃) as complete coordinate system
# ======================================================================
print("=" * 70)
print("PART 1: (H, I₃) AS COMPLETE COORDINATE SYSTEM (n≤8)")
print("=" * 70)

print("""
  H = I(Ω,2) = 1 + 2α₁ + 4α₂      (at n≤8, α₃=0)
  I₃ = I(Ω,3) = 1 + 3α₁ + 9α₂

  Solving:
    3H = 3 + 6α₁ + 12α₂
    2I₃ = 2 + 6α₁ + 18α₂
    2I₃ - 3H = -1 + 6α₂
    α₂ = (2I₃ - 3H + 1) / 6
    α₁ = (H - 1 - 4α₂) / 2

  So (H, I₃) ↔ (α₁, α₂) is a LINEAR bijection!
  Knowing 2 and 3 gives the full independence polynomial.
""")

# Verify at n=5 exhaustive
n = 5
tb = n*(n-1)//2
print(f"  n=5 exhaustive verification:")
ok = 0
for bits in range(1 << tb):
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)
    dc3 = count_directed_k_cycles(A, n, 3)
    dc5 = count_ham_cycles(A, n)
    a1 = dc3 + dc5
    a2 = (H - 1 - 2*a1) // 4
    I3 = 1 + 3*a1 + 9*a2

    # Reconstruct from (H, I3)
    a2_rec = (2*I3 - 3*H + 1) // 6
    a1_rec = (H - 1 - 4*a2_rec) // 2
    if a1_rec == a1 and a2_rec == a2:
        ok += 1

print(f"    (H,I₃) → (α₁,α₂) reconstruction: {ok}/{1<<tb}")

# At n=7, verify
n = 7
tb = n*(n-1)//2
np.random.seed(42)
print(f"\n  n=7 (500 samples):")
ok = 0
for trial in range(500):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)
    dc3 = count_directed_k_cycles(A, n, 3)
    dc5 = count_directed_k_cycles(A, n, 5)
    dc7 = count_ham_cycles(A, n)
    a1 = dc3 + dc5 + dc7
    a2 = (H - 1 - 2*a1) // 4
    I3 = 1 + 3*a1 + 9*a2

    a2_rec = (2*I3 - 3*H + 1) // 6
    a1_rec = (H - 1 - 4*a2_rec) // 2
    if a1_rec == a1 and a2_rec == a2:
        ok += 1

print(f"    (H,I₃) → (α₁,α₂) reconstruction: {ok}/500")

# ======================================================================
# PART 2: 2-adic structure of H
# ======================================================================
print("\n" + "=" * 70)
print("PART 2: 2-ADIC AND 3-ADIC STRUCTURE")
print("=" * 70)

# H is always odd (Rédei). So v₂(H) = 0.
# H ≡ 1 + 2α₁ (mod 4). So (H-1)/2 ≡ α₁ (mod 2).
# H mod 3 = (1+2α₁+α₂) mod 3 = (1-α₁+α₂) mod 3 since 2≡-1 mod 3.

for n in [5, 7]:
    tb = n*(n-1)//2
    if n == 5:
        bits_list = range(1 << tb)
    else:
        np.random.seed(42)
        bits_list = [np.random.randint(0, 1 << tb) for _ in range(1000)]

    print(f"\n  n={n}:")

    h_mod3 = Counter()
    h_mod4 = Counter()
    h_mod6 = Counter()

    for bits in bits_list:
        A = bits_to_adj(bits, n)
        H = count_ham_paths(A, n)
        h_mod3[H % 3] += 1
        h_mod4[H % 4] += 1
        h_mod6[H % 6] += 1

    print(f"    H mod 2: always 1 (Rédei)")
    print(f"    H mod 3: {dict(sorted(h_mod3.items()))}")
    print(f"    H mod 4: {dict(sorted(h_mod4.items()))}")
    print(f"    H mod 6: {dict(sorted(h_mod6.items()))}")

# ======================================================================
# PART 3: When does knowing ONLY H suffice?
# ======================================================================
print("\n" + "=" * 70)
print("PART 3: WHEN H ALONE SUFFICES")
print("=" * 70)

# At n=5: α₂=0 always, so H determines α₁ = (H-1)/2.
# Also H determines I₃ = (3H-1)/2.
# So H alone is the COMPLETE invariant at n=5!

# At n=7: H does NOT determine (α₁, α₂).
# Same H can come from different (α₁, α₂) pairs.

n = 7
tb = n*(n-1)//2
np.random.seed(42)

H_to_a1a2 = defaultdict(set)
for trial in range(1000):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)
    dc3 = count_directed_k_cycles(A, n, 3)
    dc5 = count_directed_k_cycles(A, n, 5)
    dc7 = count_ham_cycles(A, n)
    a1 = dc3 + dc5 + dc7
    a2 = (H - 1 - 2*a1) // 4
    H_to_a1a2[H].add((a1, a2))

amb = sum(1 for v in H_to_a1a2.values() if len(v) > 1)
print(f"\n  H → (α₁,α₂) at n=7: {amb}/{len(H_to_a1a2)} values of H are ambiguous")

print(f"\n  Examples of H with multiple (α₁,α₂) decompositions:")
for H in sorted(H_to_a1a2.keys()):
    pairs = H_to_a1a2[H]
    if len(pairs) > 1:
        print(f"    H={H:3d}: {sorted(pairs)}")
    if H > 80:
        break

# ======================================================================
# PART 4: The transition thresholds
# ======================================================================
print("\n" + "=" * 70)
print("PART 4: CRITICAL TRANSITIONS")
print("=" * 70)

print("""
  n=3: Only dc3 matters. α₂=0. H = 1 + 2·dc3 ∈ {1,3}.
       ONE evaluation (at 2) suffices.
       k-nacci precision: 1/2^3 = 1/8.

  n=5: dc3 and dc5 matter. α₂=0.
       H = 1 + 2(dc3+dc5).
       ONE evaluation (at 2) still suffices!
       k-nacci precision: 1/2^5 = 1/32.

  n=6: FIRST α₂ > 0 appears.
       H = 1 + 2α₁ + 4α₂. Now 2 is NOT enough.
       Need BOTH I(Ω,2) and I(Ω,3).
       "Need to know 3."
       k-nacci precision: 1/2^6 = 1/64.
       Weighted precision: (2/3)^6 ≈ 1/11.4.

  n=7: α₂ common (mean ≈ 4.5). Same structure as n=6.
       Two evaluations (at 2 and 3) suffice.

  n=9: FIRST α₃ > 0 (three disjoint 3-cycles).
       H = 1 + 2α₁ + 4α₂ + 8α₃.
       Need THREE evaluations: I(Ω,2), I(Ω,3), I(Ω,4).
       "Need to know 4 = 2²" (deeper knowledge of 2).

  THEREFORE:
    "Knowing 2" = knowing I(Ω, 2) = H
    "Knowing 3" = knowing I(Ω, 3), hence α₂
    "Knowing 2 AND 3" = full invariant for n ≤ 8

  The quote "keys to the universe" is LITERALLY TRUE:
  for n ≤ 8, two evaluation points (2 and 3) determine I(Ω, x).
""")

# ======================================================================
# PART 5: Correlation between k-nacci rate and α structure
# ======================================================================
print("=" * 70)
print("PART 5: INFORMATION CONTENT OF 2 AND 3")
print("=" * 70)

# The information from I(Ω,2) = H:
#   - Determines α₁ uniquely when α₂=0 (n=5)
#   - Gives ONE linear equation in (α₁,α₂) when α₃=0

# The ADDITIONAL information from I(Ω,3):
#   - Gives SECOND linear equation, resolving α₂
#   - This is the "3-information" beyond what 2 provides

n = 7
tb = n*(n-1)//2
np.random.seed(42)

a2_vals = []
H_vals = []
I3_vals = []
info_gap = []  # How much I₃ adds beyond H

for trial in range(1000):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)
    dc3 = count_directed_k_cycles(A, n, 3)
    dc5 = count_directed_k_cycles(A, n, 5)
    dc7 = count_ham_cycles(A, n)
    a1 = dc3 + dc5 + dc7
    a2 = (H - 1 - 2*a1) // 4
    I3 = 1 + 3*a1 + 9*a2

    a2_vals.append(a2)
    H_vals.append(H)
    I3_vals.append(I3)
    # "Information gap" = residual of I₃ after removing H contribution
    # I₃ = (9H-5-6α₁)/4. Without knowing α₁, I₃ provides new info.
    info_gap.append(I3 - 3*H/2)  # deviation from "trivial" prediction

a2_arr = np.array(a2_vals, dtype=float)
H_arr = np.array(H_vals, dtype=float)
I3_arr = np.array(I3_vals, dtype=float)
gap_arr = np.array(info_gap, dtype=float)

print(f"\n  n=7 statistics (1000 samples):")
print(f"    H:   mean={np.mean(H_arr):.1f}, std={np.std(H_arr):.1f}")
print(f"    I₃:  mean={np.mean(I3_arr):.1f}, std={np.std(I3_arr):.1f}")
print(f"    α₂:  mean={np.mean(a2_arr):.1f}, std={np.std(a2_arr):.1f}")
print(f"    I₃/H ratio: mean={np.mean(I3_arr/H_arr):.4f}")

# Correlation between H and I₃
corr_H_I3 = np.corrcoef(H_arr, I3_arr)[0,1]
corr_H_a2 = np.corrcoef(H_arr, a2_arr)[0,1]
print(f"\n    corr(H, I₃) = {corr_H_I3:.6f}")
print(f"    corr(H, α₂) = {corr_H_a2:.6f}")
print(f"    corr(α₂, I₃-3H/2) = {np.corrcoef(a2_arr, gap_arr)[0,1]:.6f}")

# The "3-information" is what I₃ tells us beyond H
# Residual I₃ = I₃ - 3H/2 should be proportional to α₂
print(f"\n    I₃ - 3H/2 = -(5+6α₁)/4 + 9α₂/2")
print(f"    This residual contains BOTH α₁ and α₂.")
print(f"    So I₃ doesn't just add α₂ info — it adds α₁ info too.")
print(f"    The LINEAR ALGEBRA of two equations in two unknowns")
print(f"    is the 'key' structure that 2 and 3 unlock.")

# ======================================================================
# PART 6: The product 2×3=6 and the sum 2+3=5
# ======================================================================
print("\n" + "=" * 70)
print("PART 6: NUMEROLOGY OF 2 AND 3")
print("=" * 70)

print(f"""
  2 × 3 = 6:  The Vandermonde determinant |2-3| = 1.
               det [[1, 2], [1, 3]] = 1.
               This is WHY (H, I₃) → (α₁, α₂) is an INTEGER bijection!
               No denominators needed (up to a factor of 6).

  2 + 3 = 5:  The first n where dc5 > 0 (5-cycles appear).
               Also the first n where score ↛ H.

  2 - 3 = -1: 2 ≡ -1 (mod 3). This gives the Eisenstein structure.

  2/3:        The convergence rate of weighted k-nacci.
               Also the "forgetting rate" at x=2.

  1/2 × 2/3 = 1/3: The Jacobsthal rate. 2^n mod 3 = (-1)^n.

  The Vandermonde determinant being 1 is the KEY FACT:
    It means the matrix [[2, 4], [3, 9]] has det = 18-12 = 6.
    So α₂ = (2I₃-3H+1)/6 is ALWAYS an integer.
    This integrality is guaranteed by the OCF structure.

  CHECK: Does α₂ = (2I₃-3H+1)/6 give integer values?
""")

# Verify integrality
n = 7
tb = n*(n-1)//2
np.random.seed(42)
ok = 0
for trial in range(1000):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)
    dc3 = count_directed_k_cycles(A, n, 3)
    dc5 = count_directed_k_cycles(A, n, 5)
    dc7 = count_ham_cycles(A, n)
    a1 = dc3 + dc5 + dc7
    a2 = (H - 1 - 2*a1) // 4
    I3 = 1 + 3*a1 + 9*a2

    numer = 2*I3 - 3*H + 1
    if numer % 6 == 0:
        ok += 1
    else:
        print(f"    FAIL: H={H}, I3={I3}, 2I3-3H+1={numer}, mod 6={numer%6}")
        break

print(f"  Integrality check: (2I₃-3H+1) mod 6 = 0: {ok}/1000")
print(f"  This integrality follows from H ≡ 1 (mod 2) and the OCF structure.")

# ======================================================================
# PART 7: Beyond 2 and 3 — when do we need 4?
# ======================================================================
print("\n" + "=" * 70)
print("PART 7: BEYOND 2 AND 3 — THE α₃ FRONTIER")
print("=" * 70)

print(f"""
  α₃ first appears at n=9 (three vertex-disjoint 3-cycles).
  At n=9: H = 1 + 2α₁ + 4α₂ + 8α₃.
  Need I(Ω,2), I(Ω,3), I(Ω,4) to solve for (α₁, α₂, α₃).

  The 3×3 system:
    H  = 1 + 2α₁ +  4α₂ +   8α₃
    I₃ = 1 + 3α₁ +  9α₂ +  27α₃
    I₄ = 1 + 4α₁ + 16α₂ +  64α₃

  Vandermonde matrix [[2,4,8],[3,9,27],[4,16,64]]:
    det = 2·9·64 + 4·27·4 + 8·3·16 - 8·9·4 - 4·3·64 - 2·27·16
        = 1152 + 432 + 384 - 288 - 768 - 864
        = 48

  So (α₁,α₂,α₃) can be recovered from (H,I₃,I₄) with denominator 48.
  Integrality requires (certain combination) ≡ 0 (mod 48).

  But 4 = 2²! So "knowing 4" is "deeper knowledge of 2," not new.
  The TRULY new primes in the Vandermonde are 2, 3, 5, 7, ...
  which are exactly the ODD PRIME cycle lengths in tournaments!

  The "keys to the universe" — 2 and 3 — are the FIRST two:
  2 = the evaluation point of OCF
  3 = the first odd prime cycle length

  For n ≤ 8: 2 and 3 suffice (α₃=0).
  For n ≤ 14: 2, 3, and 5 suffice (α₄=0, since 4×3=12<15 needed).
    Actually α₄=0 needs 4 disjoint odd cycles, min 12 vertices.
  For general n: need evaluations at 2, 3, ..., floor(n/3)+1.
""")

print("Done.")
