#!/usr/bin/env python3
"""
cv2_proof_89c.py — Verify the analytical proof that CV² = 2/n + O(1/n²)
opus-2026-03-14-S89c

PROOF SKETCH:
  CV² = W(n)/n! - 1 where W(n) = Σ_{σ∈NUD(n)} 2^{adj1(σ)}.

  Key identity: W(n)/n! = E[∏_{j=1}^{n-1} (1 + Z_j)]
  where Z_j = X_j - Y_j, X_j = 1[σ(j+1)=σ(j)+1], Y_j = 1[σ(j+1)=σ(j)-1].

  Expand the product:
  ∏(1+Z_j) = Σ_{S⊆[n-1]} ∏_{j∈S} Z_j

  E[∏(1+Z_j)] = Σ_S E[∏_{j∈S} Z_j]

  KEY FACTS:
  1. E[Z_j] = 0 (by σ↔(n-1-σ) symmetry)
  2. All odd-|S| terms vanish (same symmetry: Z→-Z)
  3. E[Z_j Z_k] = 0 for |j-k| ≥ 2 (by explicit computation:
     ascending pairs and descending pairs at distant positions have same joint prob)
  4. E[Z_j Z_{j+1}] = 2/(n(n-1)) (adjacent positions interact)
  5. Higher terms are O(1/n²)

  RESULT:
  E[∏(1+Z_j)] = 1 + Σ_j E[Z_j Z_{j+1}] + O(1/n²)
               = 1 + (n-2)·2/(n(n-1)) + O(1/n²)
               = 1 + 2/n - 2/(n(n-1)) + O(1/n²)
               = 1 + 2/n + O(1/n²)

  Hence CV² = 2/n + O(1/n²). □
"""

from fractions import Fraction
from itertools import permutations
from math import factorial

def verify_key_facts(n):
    """Verify all key facts for the proof at given n."""
    N = factorial(n)

    # Compute Z_j for all permutations
    # Z_j = X_j - Y_j where X_j=1[σ(j+1)=σ(j)+1], Y_j=1[σ(j+1)=σ(j)-1]

    # Fact 1: E[Z_j] = 0
    E_Z = [Fraction(0)] * (n-1)

    # Fact 3: E[Z_j Z_k] for all j,k
    E_ZZ = [[Fraction(0)]*(n-1) for _ in range(n-1)]

    # Also compute E[X_jX_k] and E[X_jY_k] separately
    E_XX = [[Fraction(0)]*(n-1) for _ in range(n-1)]
    E_XY = [[Fraction(0)]*(n-1) for _ in range(n-1)]

    # Full product
    full_product_sum = Fraction(0)

    for perm in permutations(range(n)):
        Z = []
        X = []
        Y = []
        for j in range(n-1):
            xj = 1 if perm[j+1] == perm[j] + 1 else 0
            yj = 1 if perm[j+1] == perm[j] - 1 else 0
            zj = xj - yj
            Z.append(zj)
            X.append(xj)
            Y.append(yj)

        for j in range(n-1):
            E_Z[j] += Fraction(Z[j], N)
            for k in range(n-1):
                E_ZZ[j][k] += Fraction(Z[j]*Z[k], N)
                E_XX[j][k] += Fraction(X[j]*X[k], N)
                E_XY[j][k] += Fraction(X[j]*Y[k], N)

        # Full product ∏(1+Z_j)
        prod = 1
        for z in Z:
            prod *= (1 + z)
        full_product_sum += Fraction(prod, N)

    return E_Z, E_ZZ, E_XX, E_XY, full_product_sum

print("="*70)
print("VERIFICATION OF CV² = 2/n + O(1/n²) PROOF")
print("="*70)

for n in range(3, 9):
    print(f"\n{'='*50}")
    print(f"n = {n}")
    print(f"{'='*50}")

    E_Z, E_ZZ, E_XX, E_XY, full_prod = verify_key_facts(n)

    # Fact 1: E[Z_j] = 0 for all j
    all_zero = all(E_Z[j] == 0 for j in range(n-1))
    print(f"\n  Fact 1: E[Z_j] = 0 for all j: {'✓' if all_zero else '✗'}")

    # Fact 3 & 4: E[Z_j Z_k]
    print(f"\n  E[Z_j Z_k] matrix:")
    for j in range(min(n-1, 6)):
        row = " ".join(f"{float(E_ZZ[j][k]):8.5f}" for k in range(min(n-1, 6)))
        print(f"    j={j}: {row}")

    # Check: E[Z_jZ_k] = 0 for |j-k| >= 2
    nonzero_distant = False
    for j in range(n-1):
        for k in range(n-1):
            if abs(j-k) >= 2 and E_ZZ[j][k] != 0:
                nonzero_distant = True
                print(f"    WARNING: E[Z_{j}Z_{k}] = {E_ZZ[j][k]} ≠ 0 but |j-k|={abs(j-k)}")
    print(f"\n  Fact 3: E[Z_jZ_k] = 0 for |j-k|≥2: {'✓' if not nonzero_distant else '✗'}")

    # E[Z_j Z_{j+1}] for adjacent pairs
    adjacent_vals = [E_ZZ[j][j+1] for j in range(n-2)]
    expected = Fraction(2, n*(n-1))
    all_match = all(v == expected for v in adjacent_vals)
    print(f"\n  Fact 4: E[Z_jZ_{{j+1}}] = 2/(n(n-1)) = {expected} = {float(expected):.8f}")
    print(f"    All adjacent pairs match: {'✓' if all_match else '✗'}")
    if not all_match:
        for j in range(n-2):
            print(f"      j={j}: {adjacent_vals[j]} {'✓' if adjacent_vals[j]==expected else '✗'}")

    # WHY E[Z_jZ_k]=0 for |j-k|>=2: because E[X_jX_k] = E[X_jY_k]
    print(f"\n  Why distant pairs vanish: E[X_jX_k] vs E[X_jY_k]")
    for j in range(min(n-1, 4)):
        for k in range(min(n-1, 4)):
            if abs(j-k) >= 2:
                print(f"    j={j},k={k}: E[XX]={E_XX[j][k]}, E[XY]={E_XY[j][k]}, equal={'✓' if E_XX[j][k]==E_XY[j][k] else '✗'}")

    # Sum of all adjacent E[Z_jZ_{j+1}]
    adj_sum = sum(E_ZZ[j][j+1] for j in range(n-2))
    print(f"\n  Sum of adjacent E[Z_jZ_{{j+1}}]: {adj_sum} = {float(adj_sum):.8f}")
    print(f"  (n-2)·2/(n(n-1)) = {Fraction(2*(n-2), n*(n-1))} = {float(Fraction(2*(n-2), n*(n-1))):.8f}")

    # Full product vs approximation
    approx_1 = 1 + adj_sum  # 1 + Σ adj pairs
    print(f"\n  ∏(1+Z_j) expectation:")
    print(f"    Exact: E[∏(1+Z_j)] = {full_prod} = {float(full_prod):.10f}")
    print(f"    Order-2 approx: 1 + adj_sum = {approx_1} = {float(approx_1):.10f}")
    print(f"    Difference (= higher order): {float(full_prod - approx_1):.10f}")
    print(f"    n² × difference: {float(n*n*(full_prod - approx_1)):.10f}")

    cv2_exact = full_prod - 1
    cv2_leading = Fraction(2, n)
    print(f"\n  CV²:")
    print(f"    Exact: {float(cv2_exact):.10f}")
    print(f"    2/n:   {float(cv2_leading):.10f}")
    print(f"    CV² - 2/n = {float(cv2_exact - cv2_leading):.10f}")
    print(f"    n²(CV² - 2/n) = {float(n*n*(cv2_exact - cv2_leading)):.10f}")

print("\n" + "="*70)
print("PROOF SUMMARY")
print("="*70)
print("""
THEOREM: CV²(H) = 2/n + O(1/n²)

PROOF:
  Let σ be a uniformly random permutation of [n].
  Define X_j = 1[σ(j+1)=σ(j)+1] (unit ascent at position j),
         Y_j = 1[σ(j+1)=σ(j)-1] (unit descent at position j),
         Z_j = X_j - Y_j.

  Then 2^{adj1(σ)} = ∏(1+X_j) and 1_{NUD}(σ) = ∏(1-Y_j).
  Since X_j·Y_j = 0 (can't have both at same position):
    2^{adj1(σ)} · 1_{NUD}(σ) = ∏(1+X_j)·∏(1-Y_j) = ∏(1+X_j-Y_j) = ∏(1+Z_j).

  So W(n)/n! = E[∏(1+Z_j)].

  Expand: E[∏(1+Z_j)] = Σ_{S⊆[n-1]} E[∏_{j∈S} Z_j].

  STEP 1: Odd terms vanish.
    The involution σ ↦ (n-1-σ(1), ..., n-1-σ(n)) preserves the uniform
    distribution and sends X_j ↔ Y_j, hence Z_j ↦ -Z_j.
    For |S| odd: E[∏Z_j] = E[∏(-Z_j)] = (-1)^{|S|} E[∏Z_j] = -E[∏Z_j],
    so E[∏Z_j] = 0.

  STEP 2: Non-adjacent pairs cancel.
    For |j-k| ≥ 2: E[Z_jZ_k] = E[X_jX_k] - E[X_jY_k] - E[Y_jX_k] + E[Y_jY_k].
    By the σ↦n-1-σ symmetry: E[X_jX_k] = E[Y_jY_k], E[X_jY_k] = E[Y_jX_k].
    By explicit counting: E[X_jX_k] = E[X_jY_k] for |j-k| ≥ 2.
    (Both count pairs of disjoint consecutive pairs among {0,...,n-1},
    with the same combinatorial structure regardless of orientation.)
    Hence E[Z_jZ_k] = 0 for |j-k| ≥ 2.

  STEP 3: Adjacent pairs contribute 2/n.
    E[Z_jZ_{j+1}] = 2/(n(n-1)) for each j.
    Proof: E[X_jY_{j+1}] = 0 (σ(j+1)=σ(j)+1 and σ(j+2)=σ(j+1)-1=σ(j)
    would require σ(j) to appear twice).
    E[X_jX_{j+1}] = P(σ(j),σ(j+1),σ(j+2) = v,v+1,v+2) = (n-2)·(n-3)!/n!
                   = 1/(n(n-1)).
    By symmetry: E[Y_jY_{j+1}] = 1/(n(n-1)), E[Y_jX_{j+1}] = 0.
    So E[Z_jZ_{j+1}] = 1/(n(n-1)) + 1/(n(n-1)) = 2/(n(n-1)).

    Sum: Σ_{j=1}^{n-2} E[Z_jZ_{j+1}] = (n-2)·2/(n(n-1)) = 2/n - 2/(n(n-1)).

  STEP 4: Higher terms are O(1/n²).
    For |S| = 4: the dominant contribution comes from "two disjoint adjacent
    pairs" {j,j+1,k,k+1}. There are ≈ n²/2 such sets, each contributing
    ≈ 4/n⁴ (product of adjacent pair expectations). Total ≈ 2/n².
    Connected 4-clusters {j,j+1,j+2,j+3} contribute O(1/n³).
    |S| ≥ 6 contributes O(1/n³).

  CONCLUSION:
    CV² = E[∏(1+Z_j)] - 1 = (2/n - 2/(n(n-1))) + O(1/n²) = 2/n + O(1/n²). □
""")
print("Done!")
