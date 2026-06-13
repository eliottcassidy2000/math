#!/usr/bin/env python3
"""
prove_linear_H_S2_s97.py — opus-2026-03-21-S97

PROVING THE LINEAR FORMULA: H = H_max - α·S₂

Discovery from spectral_kurtosis_deep_s97.py:
  n=5: H = 15.0 - 1.5·S₂  (R²=0.947, EXACT regression)
  n=6: H = 45.0 - 3.0·S₂  (R²=0.923, EXACT regression)

The coefficient pattern:
  n=5: α = 1.5 = 3/2, H_max = 15
  n=6: α = 3.0, H_max = 45

Is this H_max/S₂_max?
  n=5: S₂_max = 10 (transitive), H_max/S₂_max = 15/10 = 3/2 ✓
  n=6: S₂_max = 15 (transitive), H_max/S₂_max = 45/15 = 3 ✓

If H = H_max·(1 - S₂/S₂_max), then α = H_max/S₂_max.

Check: at S₂ = S₂_max (transitive), H = 0? No, H(transitive) = 1.
So the formula is approximate: H ≈ H_max·(1 - S₂/S₂_max) ≈ 0 at transitive.
But H(transitive) = 1, not 0. And 15·(1 - 10/10) = 0 ≠ 1.

So the linear formula is NOT exact. It's the BEST LINEAR FIT.

Let's investigate what the EXACT formula is.

KEY INSIGHT FROM OCF:
  H = I(Ω, 2) = 1 + 2α₁ + 4α₂ + ...
  α₁ = #{odd cycles}
  For 3-cycles only: α₁ ≈ c₃

  And c₃ is directly related to S₂ via:
  c₃ = C(n,3)/4 - S₂/4 · something?

  Actually, the classical formula:
  c₃ = C(n,3)/4 - (1/4)Σᵢsᵢ(sᵢ-1)/2... no.

  The standard formula for c₃:
  c₃ = C(n,3) - Σᵢ C(sᵢ, 2)

  where the sum counts the number of 3-element subsets that are NOT 3-cycles.
  (Each 3-subset is either a 3-cycle or a transitive triple, and transitive
  triples contribute C(sᵢ, 2) = sᵢ(sᵢ-1)/2 per vertex via dominance.)

  Wait, more carefully:
  Number of transitive triples = Σᵢ C(sᵢ, 2) = (1/2)Σᵢ sᵢ(sᵢ-1)
  And c₃ = C(n,3) - Σᵢ C(sᵢ, 2)

  Now Σᵢ C(sᵢ, 2) = (1/2)(Σsᵢ² - Σsᵢ) = (1/2)(Σsᵢ² - C(n,2))
                    = (1/2)(S₂ + n((n-1)/2)² - C(n,2))

  Since S₂ = Σ(sᵢ - (n-1)/2)² = Σsᵢ² - n(n-1)²/4.
  So Σsᵢ² = S₂ + n(n-1)²/4.
  And Σᵢ C(sᵢ,2) = (1/2)(S₂ + n(n-1)²/4 - n(n-1)/2)
                   = (1/2)(S₂ + n(n-1)(n-3)/4)

  c₃ = C(n,3) - (1/2)(S₂ + n(n-1)(n-3)/4)
     = n(n-1)(n-2)/6 - S₂/2 - n(n-1)(n-3)/8

  So c₃ is LINEAR in S₂ with coefficient -1/2!

  c₃ = [n(n-1)(n-2)/6 - n(n-1)(n-3)/8] - S₂/2
     = n(n-1)·[4(n-2) - 3(n-3)] / 24 - S₂/2
     = n(n-1)·(4n-8-3n+9) / 24 - S₂/2
     = n(n-1)·(n+1) / 24 - S₂/2
     = n(n-1)(n+1)/24 - S₂/2

  BEAUTIFUL! So c₃ = C(n+1, 3)/4 - S₂/2.

  Wait: n(n-1)(n+1)/24 = (n+1)n(n-1)/24 = (n+1)!/((n-2)!·24) = C(n+1,3)/4.
  Yes! C(n+1,3) = (n+1)n(n-1)/6, so C(n+1,3)/4 = (n+1)n(n-1)/24.

  THEOREM: c₃ = C(n+1, 3)/4 - S₂/2.

  At n=5: c₃ = C(6,3)/4 - S₂/2 = 20/4 - S₂/2 = 5 - S₂/2.
  Check: S₂=0 (regular) → c₃=5 ✓. S₂=10 (transitive) → c₃=0 ✓.

  At n=6: c₃ = C(7,3)/4 - S₂/2 = 35/4 - S₂/2 = 8.75 - S₂/2.
  Check: S₂=0 → c₃=8.75? But c₃ must be integer! And at n=6 (even),
  there is no regular tournament (S₂ > 0 always).
  Min S₂ at n=6: score (2,2,2,3,3,3) has S₂ = 6·(1/4) = 3/2? Let me compute.
  Actually at n=6, (n-1)/2 = 5/2. Scores sum to C(6,2) = 15.
  Most balanced: (2,2,2,3,3,3), S₂ = 3·(2-5/2)² + 3·(3-5/2)² = 3·1/4 + 3·1/4 = 3/2.
  c₃ = 8.75 - 3/4 = 8. Check: c₃ for a specific (2,2,2,3,3,3) tournament?

  Let me verify numerically.
"""

import sys
import numpy as np
from itertools import combinations
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

def all_tournaments(n):
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(pairs)
    for bits in range(2**m):
        A = np.zeros((n, n), dtype=int)
        for k, (i, j) in enumerate(pairs):
            if (bits >> k) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        yield A

def count_hp(A, n):
    dp = defaultdict(int)
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[(mask, v)] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    dp[(mask | (1 << w), w)] += dp[(mask, v)]
    full = (1 << n) - 1
    return sum(dp[(full, v)] for v in range(n))

def count_c3(A, n):
    c3 = 0
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]:
            c3 += 1
        if A[i][k] and A[k][j] and A[j][i]:
            c3 += 1
    return c3

print("=" * 72)
print("  PROVING THE LINEAR FORMULA: H = H_max - α·S₂")
print("  opus-2026-03-21-S97")
print("=" * 72)

# ========================================================================
# PART 1: Verify c₃ = C(n+1,3)/4 - S₂/2
# ========================================================================
print("\n" + "=" * 72)
print("  PART 1: THEOREM — c₃ = C(n+1,3)/4 - S₂/2")
print("=" * 72)

from math import comb

for n in [3, 4, 5, 6, 7]:
    print(f"\n--- n = {n} ---")
    predicted_formula = comb(n + 1, 3) / 4
    print(f"  c₃ = {predicted_formula} - S₂/2")

    if n <= 6:
        all_ok = True
        for A in all_tournaments(n):
            c3 = count_c3(A, n)
            scores = A.sum(axis=1)
            S2 = np.sum((scores - (n - 1) / 2) ** 2)
            predicted = predicted_formula - S2 / 2
            if abs(c3 - predicted) > 0.01:
                all_ok = False
                print(f"  FAIL: c₃={c3}, S₂={S2}, predicted={predicted}")
                break

        if all_ok:
            print(f"  VERIFIED for all {2**comb(n,2)} tournaments ✓")
    else:
        # Sample at n=7
        import random
        random.seed(42)
        all_ok = True
        for _ in range(1000):
            A = np.zeros((n, n), dtype=int)
            for i in range(n):
                for j in range(i + 1, n):
                    if random.random() < 0.5:
                        A[i][j] = 1
                    else:
                        A[j][i] = 1
            c3 = count_c3(A, n)
            scores = A.sum(axis=1)
            S2 = np.sum((scores - (n - 1) / 2) ** 2)
            predicted = predicted_formula - S2 / 2
            if abs(c3 - predicted) > 0.01:
                all_ok = False
                break

        print(f"  VERIFIED for 1000 random tournaments ✓" if all_ok
              else f"  FAILED!")

# ========================================================================
# PART 2: From c₃ to H via OCF
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 2: FROM c₃ TO H — THE LINEAR APPROXIMATION")
print("=" * 72)
print("""
  OCF: H = 1 + 2α₁ + 4α₂ + ...
  where α₁ = #{odd cycles}, α₂ = #{independent pairs}, etc.

  For SMALL n (n=3,4,5), the only odd cycles are 3-cycles and 5-cycles.
  α₁ = c₃ + c₅ + c₇ + ...

  At n=5: α₁ = c₃ + c₅. Max c₅ per tournament = 0 or more.
  Actually the n=5 regular tournament has c₅ = 24/5 = ... let me check.

  The key question: how does H relate to c₃ linearly?

  We know: c₃ = C(n+1,3)/4 - S₂/2  (proved above)
  And regression found: H ≈ H_max - α·S₂ = H_max - 2α·c₃ + ... ?

  If H is approximately linear in c₃:
    H ≈ a + b·c₃
  Then H ≈ a + b·(C(n+1,3)/4 - S₂/2) = (a + b·C(n+1,3)/4) - b/2·S₂

  So α = b/2 where b is the c₃ coefficient of H.

  Let's fit H vs c₃ directly.
""")

for n in [3, 4, 5, 6]:
    print(f"\n--- n = {n} ---")
    data = []
    for A in all_tournaments(n):
        H = count_hp(A, n)
        c3 = count_c3(A, n)
        scores = A.sum(axis=1)
        S2 = np.sum((scores - (n - 1) / 2) ** 2)
        data.append((H, c3, S2))

    Hs = np.array([d[0] for d in data])
    c3s = np.array([d[1] for d in data])
    S2s = np.array([d[2] for d in data])

    # Linear fit H vs c3
    X = np.column_stack([c3s, np.ones(len(data))])
    beta, _, _, _ = np.linalg.lstsq(X, Hs, rcond=None)
    H_pred = X @ beta
    R2 = 1 - np.sum((Hs - H_pred) ** 2) / np.sum((Hs - np.mean(Hs)) ** 2)

    print(f"  H ≈ {beta[0]:.4f}·c₃ + {beta[1]:.4f}")
    print(f"  R² = {R2:.6f}")
    print(f"  Coefficient of c₃: {beta[0]:.4f}")
    print(f"  Expected α = coeff/2 = {beta[0]/2:.4f}")

    # Substituting c₃ = C(n+1,3)/4 - S₂/2:
    # H ≈ b·c₃ + a = b·(C(n+1,3)/4 - S₂/2) + a
    # = (a + b·C(n+1,3)/4) - (b/2)·S₂
    predicted_intercept = beta[1] + beta[0] * comb(n + 1, 3) / 4
    predicted_slope = -beta[0] / 2
    print(f"  Predicted H = {predicted_intercept:.4f} + ({predicted_slope:.4f})·S₂")

    # Compare with direct S₂ fit
    X2 = np.column_stack([S2s, np.ones(len(data))])
    beta2, _, _, _ = np.linalg.lstsq(X2, Hs, rcond=None)
    print(f"  Direct S₂ fit: H = {beta2[0]:.4f}·S₂ + {beta2[1]:.4f}")
    print(f"  Match: slope {abs(predicted_slope - beta2[0]) < 0.01}, "
          f"intercept {abs(predicted_intercept - beta2[1]) < 0.01}")

# ========================================================================
# PART 3: The coefficient b = slope of H vs c₃
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 3: THE H-vs-c₃ COEFFICIENT")
print("=" * 72)
print("""
  From the regression, the coefficient of c₃ in H is:
""")

for n in [3, 4, 5, 6, 7]:
    print(f"\n--- n = {n} ---")
    data = []

    if n <= 6:
        for A in all_tournaments(n):
            data.append((count_hp(A, n), count_c3(A, n)))
    else:
        import random
        random.seed(42)
        for _ in range(5000):
            A = np.zeros((n, n), dtype=int)
            for i in range(n):
                for j in range(i + 1, n):
                    if random.random() < 0.5:
                        A[i][j] = 1
                    else:
                        A[j][i] = 1
            data.append((count_hp(A, n), count_c3(A, n)))

    Hs = np.array([d[0] for d in data], dtype=float)
    c3s = np.array([d[1] for d in data], dtype=float)

    X = np.column_stack([c3s, np.ones(len(data))])
    beta, _, _, _ = np.linalg.lstsq(X, Hs, rcond=None)
    R2 = 1 - np.sum((Hs - X @ beta) ** 2) / np.sum((Hs - np.mean(Hs)) ** 2)

    # b should be close to 2 (from OCF: H ≈ 1 + 2·c₃ + higher order)
    # But H = 1 + 2α₁ + 4α₂ + ..., and α₁ ≥ c₃ (includes 5-cycles etc.)
    # So the effective coefficient should be > 2 if longer cycles correlate with c₃.

    c3_max = max(c3s)
    H_at_c3_max = np.mean([h for h, c in data if c == c3_max])

    print(f"  coeff(c₃) = {beta[0]:.4f}, intercept = {beta[1]:.4f}, R² = {R2:.6f}")
    print(f"  c₃_max = {c3_max:.0f}, H at c₃_max ≈ {H_at_c3_max:.0f}")

    # OCF lower bound: H ≥ 1 + 2c₃
    # Check at c₃ = 0: H should be 1 (transitive)
    H_at_c3_0 = np.mean([h for h, c in data if c == 0])
    print(f"  H at c₃=0: {H_at_c3_0:.1f} (should be 1 for transitive)")
    print(f"  OCF lower bound 1+2·c₃_max = {1 + 2*c3_max:.0f}, "
          f"actual H_max ≈ {max(Hs):.0f}")

# ========================================================================
# PART 4: The exact formula chain
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 4: THE EXACT FORMULA CHAIN")
print("=" * 72)
print("""
  SUMMARY OF EXACT RESULTS:

  1. THEOREM (c₃-S₂): c₃ = C(n+1,3)/4 - S₂/2
     Proof: c₃ = C(n,3) - Σ_i C(s_i, 2), expand using S₂ = Σ(s_i-(n-1)/2)².

  2. EMPIRICAL: H ≈ b·c₃ + a (linear with high R²)
     Coefficient b increases with n.

  3. COMBINING: H ≈ b·(C(n+1,3)/4 - S₂/2) + a
                   = (a + b·C(n+1,3)/4) - (b/2)·S₂
                   ≡ H_max_pred - α·S₂

  4. OCF EXPLANATION:
     H = 1 + 2α₁ + 4α₂ + ... where α₁ ≈ c₃ + (correction from longer cycles)
     The linear approximation H ≈ 1 + 2c₃ + (nonlinear correction)
     works because α₂, α₃, etc. are CORRELATED with c₃.

     WHY: More 3-cycles → more odd cycles of all lengths → more independent
     sets in Ω → larger I(Ω, 2) → larger H.

  5. THE KURTOSIS REFINEMENT:
     Within a fixed score class (fixed S₂), the residual H depends on
     tr(B⁴). This captures the DISTRIBUTION of cycles, not just the count.
     Conditional corr(H, -tr(B⁴) | S₂ fixed) ≈ 0.5-0.9.

  THE TWO-PARAMETER HIERARCHY:
    Level 0: H ≈ mean(H)                           (no info)
    Level 1: H ≈ f(S₂) = H_max - α·S₂             (92-95% R²)
    Level 2: H ≈ f(S₂, tr(B⁴))                     (refines within score class)
    Level 3: H ≈ f(S₂, tr(B⁴), tr(B⁶))            (further refinement)
    ...
    Level ∞: H = I(Ω(T), 2)                        (exact, OCF)

  The LINEAR formula H ≈ H_max - α·S₂ captures >90% of H's variance
  because c₃ (which dominates H via OCF) is EXACTLY linearly related to S₂.
""")

# Exact values
print("\n  Summary table:")
print(f"  {'n':>3s} {'C(n+1,3)/4':>12s} {'α (S₂ coeff)':>14s} {'b (c₃ coeff)':>14s} {'R²':>8s}")

for n in [3, 4, 5, 6]:
    data = []
    for A in all_tournaments(n):
        H = count_hp(A, n)
        c3 = count_c3(A, n)
        scores = A.sum(axis=1)
        S2 = np.sum((scores - (n - 1) / 2) ** 2)
        data.append((H, c3, S2))

    Hs = np.array([d[0] for d in data], dtype=float)
    c3s = np.array([d[1] for d in data], dtype=float)
    S2s = np.array([d[2] for d in data], dtype=float)

    Xc = np.column_stack([c3s, np.ones(len(data))])
    bc, _, _, _ = np.linalg.lstsq(Xc, Hs, rcond=None)

    Xs = np.column_stack([S2s, np.ones(len(data))])
    bs, _, _, _ = np.linalg.lstsq(Xs, Hs, rcond=None)
    R2 = 1 - np.sum((Hs - Xs @ bs) ** 2) / np.sum((Hs - np.mean(Hs)) ** 2)

    c_val = comb(n + 1, 3) / 4
    print(f"  {n:3d} {c_val:12.2f} {bs[0]:14.4f} {bc[0]:14.4f} {R2:8.4f}")

print("\n\nDone. opus-2026-03-21-S97")
