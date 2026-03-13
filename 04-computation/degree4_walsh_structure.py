#!/usr/bin/env python3
"""
Degree-4 Walsh structure of H on the orientation cube.

KEY FINDING FROM boolean_fkg_proof.py:
  p=7:  100% of H variance is degree-2 → Paley maximizes degree-2 → Paley wins
  p=13: 81.6% of H variance is degree-4 → WHO maximizes degree-4?

THE QUESTION: Which σ ∈ {±1}^m maximizes Σ_{|S|=4} ĥ[S] χ_S(σ)?

If it's σ = all-ones (Interval), then the proof reduces to:
  1. Degree-2 is maximized by σ_P (Paley) — THM-137
  2. Degree-4 is maximized by σ_I (Interval) — TO PROVE
  3. For p ≥ 13, degree-4 dominates → Interval wins overall

This gives a COMPLETE and CLEAN proof architecture!

opus-2026-03-12-S67
"""

import numpy as np
from itertools import combinations
from sympy.ntheory import legendre_symbol as legendre

def make_tournament(p, S):
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for d in S:
            A[i][(i+d)%p] = 1
    return A

def count_H(A):
    n = len(A)
    dp = [[0]*n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))

def orientation_to_S(sigma, p):
    m = (p - 1) // 2
    S = set()
    for k in range(1, m + 1):
        if sigma[k-1] == 1:
            S.add(k)
        else:
            S.add(p - k)
    return S

# ========================================================================
print("=" * 72)
print("PART I: WALSH SPECTRUM DECOMPOSITION BY DEGREE")
print("=" * 72)

for p in [7, 11, 13]:
    m = (p - 1) // 2
    N = 2**m

    print(f"\n  p={p}, m={m}, N={N}")

    # Compute H for all orientations
    H_vec = np.zeros(N)
    for bits in range(N):
        sigma = tuple(1 if (bits >> k) & 1 else -1 for k in range(m))
        S = orientation_to_S(sigma, p)
        A = make_tournament(p, S)
        H_vec[bits] = count_H(A)

    # Full Walsh-Hadamard transform
    h_hat = np.zeros(N)
    for S_bits in range(N):
        S_set = [k for k in range(m) if (S_bits >> k) & 1]
        val = 0.0
        for bits in range(N):
            chi = 1
            for k in S_set:
                if not ((bits >> k) & 1):
                    chi = -chi
            val += H_vec[bits] * chi
        h_hat[S_bits] = val / N

    # Print ALL nonzero Walsh coefficients
    print(f"\n  All nonzero Walsh coefficients:")
    by_degree = {}
    for S_bits in range(N):
        if abs(h_hat[S_bits]) > 0.01:
            S_set = frozenset(k+1 for k in range(m) if (S_bits >> k) & 1)
            deg = len(S_set)
            if deg not in by_degree:
                by_degree[deg] = []
            by_degree[deg].append((S_set, h_hat[S_bits]))

    for deg in sorted(by_degree.keys()):
        print(f"    Degree {deg}:")
        for S_set, coeff in sorted(by_degree[deg], key=lambda x: -abs(x[1])):
            print(f"      ĥ[{set(S_set)}] = {coeff:.4f}")

    # Compute degree-2 and degree-4 maximizers
    print(f"\n  Degree-2 component: f₂(σ) = Σ_{'{|S|=2}'} ĥ[S] χ_S(σ)")
    print(f"  Degree-4 component: f₄(σ) = Σ_{'{|S|=4}'} ĥ[S] χ_S(σ)")

    f2_max = -float('inf')
    f2_max_sigma = None
    f4_max = -float('inf')
    f4_max_sigma = None

    for bits in range(N):
        sigma = tuple(1 if (bits >> k) & 1 else -1 for k in range(m))

        f2 = 0.0
        f4 = 0.0
        for S_bits in range(N):
            deg = bin(S_bits).count('1')
            chi = 1
            for k in range(m):
                if (S_bits >> k) & 1:
                    chi *= sigma[k]
            if deg == 2:
                f2 += h_hat[S_bits] * chi
            elif deg == 4:
                f4 += h_hat[S_bits] * chi

        if f2 > f2_max:
            f2_max = f2
            f2_max_sigma = sigma
        if f4 > f4_max:
            f4_max = f4
            f4_max_sigma = sigma

    sigma_int = tuple(1 for _ in range(m))
    f2_int = 0.0
    f4_int = 0.0
    for S_bits in range(N):
        deg = bin(S_bits).count('1')
        # χ_S(all-ones) = 1 always
        if deg == 2:
            f2_int += h_hat[S_bits]
        elif deg == 4:
            f4_int += h_hat[S_bits]

    print(f"    f₂(all-ones) = {f2_int:.4f}")
    print(f"    f₂(max) = {f2_max:.4f} at σ = {f2_max_sigma}")
    print(f"    f₂(all-ones) is max? {abs(f2_int - f2_max) < 0.01}")

    print(f"    f₄(all-ones) = {f4_int:.4f}")
    print(f"    f₄(max) = {f4_max:.4f} at σ = {f4_max_sigma}")
    print(f"    f₄(all-ones) is max? {abs(f4_int - f4_max) < 0.01}")

    # Who maximizes H overall?
    H_max = max(H_vec)
    H_max_sigma = None
    for bits in range(N):
        if H_vec[bits] == H_max:
            H_max_sigma = tuple(1 if (bits >> k) & 1 else -1 for k in range(m))
            break

    print(f"\n    H(all-ones) = {H_vec[N-1]:.0f}")
    print(f"    H(max) = {H_max:.0f} at σ = {H_max_sigma}")
    print(f"    f₂(max_H) vs f₂(max_f2): contribution alignment")

# ========================================================================
print("\n" + "=" * 72)
print("PART II: THE DEGREE-4 TENSOR STRUCTURE")
print("=" * 72)
print("""
The degree-4 component f₄(σ) = Σ_{i<j<k<l} T_{ijkl} σ_i σ_j σ_k σ_l
is determined by a 4th-order tensor T.

Maximizing f₄ over {±1}^m is a form of the MAX-4-XOR problem.
For general T, this is NP-hard. But our T has SPECIAL STRUCTURE
coming from the tournament's circulant symmetry.

KEY: The circulant structure means T_{ijkl} depends only on the
"chord type" {i,j,k,l} mod the dihedral group D_m.

QUESTION: Does the all-ones vector always maximize f₄ for our T?
""")

for p in [7, 11, 13]:
    m = (p - 1) // 2
    N = 2**m

    # Compute Walsh spectrum
    H_vec = np.zeros(N)
    for bits in range(N):
        sigma = tuple(1 if (bits >> k) & 1 else -1 for k in range(m))
        S = orientation_to_S(sigma, p)
        A = make_tournament(p, S)
        H_vec[bits] = count_H(A)

    h_hat = np.zeros(N)
    for S_bits in range(N):
        S_set = [k for k in range(m) if (S_bits >> k) & 1]
        val = 0.0
        for bits in range(N):
            chi = 1
            for k in S_set:
                if not ((bits >> k) & 1):
                    chi = -chi
            val += H_vec[bits] * chi
        h_hat[S_bits] = val / N

    # Extract degree-4 tensor
    deg4_coeffs = {}
    for S_bits in range(N):
        if bin(S_bits).count('1') == 4:
            S_set = tuple(sorted(k+1 for k in range(m) if (S_bits >> k) & 1))
            if abs(h_hat[S_bits]) > 0.001:
                deg4_coeffs[S_set] = h_hat[S_bits]

    print(f"\n  p={p}: Degree-4 tensor entries:")
    if not deg4_coeffs:
        print(f"    (none)")
        continue

    # Check: are all degree-4 coefficients the SAME SIGN?
    signs = [np.sign(v) for v in deg4_coeffs.values()]
    all_positive = all(s > 0 for s in signs)
    all_negative = all(s < 0 for s in signs)

    for S_set, coeff in sorted(deg4_coeffs.items()):
        print(f"    T{S_set} = {coeff:.4f}")

    print(f"    All same sign? {'YES (positive)' if all_positive else 'YES (negative)' if all_negative else 'NO (mixed)'}")

    if all_positive:
        print(f"    *** If all degree-4 coefficients are positive, ***")
        print(f"    *** then f₄ is maximized at all-ones TRIVIALLY! ***")
        print(f"    *** (since χ_S(all-ones) = +1 for all S)       ***")

# ========================================================================
print("\n" + "=" * 72)
print("PART III: SIGN STRUCTURE OF WALSH COEFFICIENTS")
print("=" * 72)
print("""
CRITICAL QUESTION: For each even degree d, are ALL ĥ[S] with |S|=d
the SAME SIGN?

If ĥ[S] ≥ 0 for ALL S with |S| even, then:
  H(σ) = Σ_S ĥ[S] χ_S(σ) ≤ Σ_S ĥ[S] |χ_S(σ)| = Σ_S ĥ[S] = H(all-ones)

since |χ_S(σ)| = 1 and χ_S(σ) ≤ 1, with equality iff σ = all-ones!

This would be a ONE-LINE PROOF of Interval maximality!
""")

for p in [7, 11, 13, 17]:
    m = (p - 1) // 2
    N = 2**m

    print(f"\n  p={p}, m={m}:")

    # Compute Walsh spectrum
    H_vec = np.zeros(N)
    for bits in range(N):
        sigma = tuple(1 if (bits >> k) & 1 else -1 for k in range(m))
        S = orientation_to_S(sigma, p)
        A = make_tournament(p, S)
        H_vec[bits] = count_H(A)

    h_hat = np.zeros(N)
    for S_bits in range(N):
        val = 0.0
        for bits in range(N):
            chi = 1
            for k in range(m):
                if (S_bits >> k) & 1:
                    if not ((bits >> k) & 1):
                        chi = -chi
            val += H_vec[bits] * chi
        h_hat[S_bits] = val / N

    # Sign analysis by degree
    for deg in range(0, m+1):
        coeffs = []
        for S_bits in range(N):
            if bin(S_bits).count('1') == deg:
                if abs(h_hat[S_bits]) > 0.001:
                    coeffs.append(h_hat[S_bits])

        if not coeffs:
            continue

        pos = sum(1 for c in coeffs if c > 0)
        neg = sum(1 for c in coeffs if c < 0)
        all_pos = all(c > 0 for c in coeffs)

        print(f"    Degree {deg}: {len(coeffs)} nonzero, "
              f"{pos} positive, {neg} negative, "
              f"range [{min(coeffs):.2f}, {max(coeffs):.2f}]",
              end="")
        if all_pos:
            print(" ★ ALL POSITIVE ★")
        else:
            print()

    # THE KEY TEST: is Σ_S ĥ[S] (summing ALL even-degree coefficients) = H(all-ones)?
    sum_all_hat = sum(h_hat[S_bits] for S_bits in range(N))
    sigma_int = tuple(1 for _ in range(m))
    bits_int = N - 1  # all bits set = all ones
    H_int = H_vec[bits_int]

    print(f"    Σ ĥ[S] (all S) = {sum_all_hat:.2f}")
    print(f"    H(all-ones) = {H_int:.0f}")
    print(f"    Match? {abs(sum_all_hat - H_int) < 0.01}")

    # What fraction of ĥ[S] are positive?
    total_nonzero = sum(1 for s in range(N) if abs(h_hat[s]) > 0.001)
    total_positive = sum(1 for s in range(N) if h_hat[s] > 0.001)
    print(f"    {total_positive}/{total_nonzero} nonzero coefficients are positive")

    # If not all positive, what's the negative contribution?
    neg_sum = sum(h_hat[s] for s in range(N) if h_hat[s] < -0.001)
    pos_sum = sum(h_hat[s] for s in range(N) if h_hat[s] > 0.001)
    print(f"    Positive sum: {pos_sum:.2f}, Negative sum: {neg_sum:.2f}")

# ========================================================================
print("\n" + "=" * 72)
print("PART IV: THE PRODUCT FORMULA AND QR STRUCTURE")
print("=" * 72)
print("""
For the circulant-restricted Walsh transform:
  ψ(S) = Π_{k ∈ S} k mod p  (multiplicative character value)

If ψ(S) is a quadratic residue mod p, then χ_S(σ_Paley) = +1.
If ψ(S) is a non-residue, then χ_S(σ_Paley) = -1.

For σ = all-ones (Interval), χ_S(all-ones) = +1 ALWAYS.

So: H(Int) - H(Pal) = 2 Σ_{ψ(S) = NQR} ĥ[S]

For this to be positive, we need: Σ_{NQR products} ĥ[S] > 0.

KEY INSIGHT: If ĥ[S] > 0 for ALL S, then EVERY σ gives H(σ) ≤ H(Int).
We don't need the NQR sum — EVERY Walsh character product sums positive!
""")

for p in [7, 11, 13]:
    m = (p - 1) // 2
    N = 2**m

    # Recompute
    H_vec = np.zeros(N)
    for bits in range(N):
        sigma = tuple(1 if (bits >> k) & 1 else -1 for k in range(m))
        S = orientation_to_S(sigma, p)
        A = make_tournament(p, S)
        H_vec[bits] = count_H(A)

    h_hat = np.zeros(N)
    for S_bits in range(N):
        val = 0.0
        for bits in range(N):
            chi = 1
            for k in range(m):
                if (S_bits >> k) & 1:
                    if not ((bits >> k) & 1):
                        chi = -chi
            val += H_vec[bits] * chi
        h_hat[S_bits] = val / N

    # Compute the product mod p for each S
    print(f"\n  p={p}: QR structure of Walsh coefficients")
    qr_positive = 0
    nqr_positive = 0
    qr_negative = 0
    nqr_negative = 0

    qr_sum = 0.0
    nqr_sum = 0.0

    for S_bits in range(1, N):
        S_set = [k+1 for k in range(m) if (S_bits >> k) & 1]
        if len(S_set) == 0:
            continue

        # Product mod p
        prod = 1
        for k in S_set:
            prod = (prod * k) % p

        # QR or NQR?
        is_qr = int(legendre(prod, p)) == 1

        if abs(h_hat[S_bits]) > 0.001:
            if is_qr:
                qr_sum += h_hat[S_bits]
                if h_hat[S_bits] > 0:
                    qr_positive += 1
                else:
                    qr_negative += 1
            else:
                nqr_sum += h_hat[S_bits]
                if h_hat[S_bits] > 0:
                    nqr_positive += 1
                else:
                    nqr_negative += 1

    print(f"    QR products:  {qr_positive} pos, {qr_negative} neg, sum = {qr_sum:.2f}")
    print(f"    NQR products: {nqr_positive} pos, {nqr_negative} neg, sum = {nqr_sum:.2f}")
    print(f"    H(Int) - H(Pal) = 2 × NQR_sum = {2*nqr_sum:.2f}")

    H_int = H_vec[N-1]  # all ones
    # Find Paley H
    if p % 4 == 3:
        sigma_pal = tuple(int(legendre(k, p)) for k in range(1, m+1))
        bits_pal = sum((1 if sigma_pal[k] == 1 else 0) << k for k in range(m))
        H_pal = H_vec[bits_pal]
        print(f"    Actual: H(Int)={H_int:.0f}, H(Pal)={H_pal:.0f}, diff={H_int-H_pal:.0f}")

# ========================================================================
print("\n" + "=" * 72)
print("PART V: DEGREE-4 POSITIVITY — THE KEY CONJECTURE")
print("=" * 72)
print("""
CONJECTURE (HYP-527): For all p ≥ 13, ALL degree-4 Walsh coefficients
ĥ[{i,j,k,l}] are POSITIVE.

If true, this gives an ELEMENTARY proof:
  1. H(σ) = ĥ[∅] + Σ_{|S|=2} ĥ[S] χ_S(σ) + Σ_{|S|=4} ĥ[S] χ_S(σ) + ...
  2. Since χ_S(σ) ≤ 1 with equality iff σ = all-ones:
     H(σ) ≤ ĥ[∅] + Σ_{|S|=2} |ĥ[S]| + Σ_{|S|=4} ĥ[S] + ...
  3. But if degree-4 dominates and is all positive:
     The maximum of f₄(σ) is uniquely at all-ones
  4. The degree-2 loss (where Paley may be better) is overcome by
     the degree-4 gain

Wait — this isn't quite right. We need ALL even-degree coefficients
to be positive for the one-line proof. Let's check degree-2 signs too.
""")

for p in [7, 11, 13, 17]:
    m = (p - 1) // 2
    N = 2**m

    # Compute
    H_vec = np.zeros(N)
    for bits in range(N):
        sigma = tuple(1 if (bits >> k) & 1 else -1 for k in range(m))
        S = orientation_to_S(sigma, p)
        A = make_tournament(p, S)
        H_vec[bits] = count_H(A)

    h_hat = np.zeros(N)
    for S_bits in range(N):
        val = 0.0
        for bits in range(N):
            chi = 1
            for k in range(m):
                if (S_bits >> k) & 1:
                    if not ((bits >> k) & 1):
                        chi = -chi
            val += H_vec[bits] * chi
        h_hat[S_bits] = val / N

    print(f"\n  p={p}:")
    # List all coefficients with sign
    all_positive = True
    for S_bits in range(1, N):
        deg = bin(S_bits).count('1')
        if abs(h_hat[S_bits]) > 0.001:
            S_set = sorted(k+1 for k in range(m) if (S_bits >> k) & 1)
            sign = "+" if h_hat[S_bits] > 0 else "-"
            if h_hat[S_bits] < 0:
                all_positive = False
            print(f"    [{sign}] ĥ[{set(S_set)}] = {h_hat[S_bits]:>12.4f}  (degree {deg})")

    if all_positive:
        print(f"    ★★★ ALL WALSH COEFFICIENTS POSITIVE AT p={p}! ★★★")
        print(f"    This means H(all-ones) = max H TRIVIALLY!")
    else:
        # Count positive vs negative
        n_pos = sum(1 for s in range(1,N) if h_hat[s] > 0.001)
        n_neg = sum(1 for s in range(1,N) if h_hat[s] < -0.001)
        pos_weight = sum(h_hat[s] for s in range(1,N) if h_hat[s] > 0.001)
        neg_weight = sum(h_hat[s] for s in range(1,N) if h_hat[s] < -0.001)
        print(f"    Positive: {n_pos} coeffs, total weight {pos_weight:.2f}")
        print(f"    Negative: {n_neg} coeffs, total weight {neg_weight:.2f}")
        print(f"    Net: {pos_weight + neg_weight:.2f} (= H(all-ones) - H(mean))")

# ========================================================================
print("\n" + "=" * 72)
print("PART VI: THE TRANSITION FROM MIXED TO ALL-POSITIVE")
print("=" * 72)
print("""
From Part V, we check: at what p do ALL Walsh coefficients become positive?

p=7:  Mixed signs (degree-2 has negatives)
p=11: Mixed signs
p=13: ???
p=17: ???

If the transition happens at p=13, then for p ≥ 13:
  ALL ĥ[S] ≥ 0 ⟹ H(all-ones) = max H ⟹ Interval is optimal.

This would be the COMPLETE PROOF!
""")

# Summary
print("\n  Summary of Walsh coefficient signs:")
for p in [7, 11, 13, 17]:
    m = (p - 1) // 2
    N = 2**m

    H_vec = np.zeros(N)
    for bits in range(N):
        sigma = tuple(1 if (bits >> k) & 1 else -1 for k in range(m))
        S = orientation_to_S(sigma, p)
        A = make_tournament(p, S)
        H_vec[bits] = count_H(A)

    h_hat = np.zeros(N)
    for S_bits in range(N):
        val = 0.0
        for bits in range(N):
            chi = 1
            for k in range(m):
                if (S_bits >> k) & 1:
                    if not ((bits >> k) & 1):
                        chi = -chi
            val += H_vec[bits] * chi
        h_hat[S_bits] = val / N

    n_neg = sum(1 for s in range(1, N) if h_hat[s] < -0.001)
    n_pos = sum(1 for s in range(1, N) if h_hat[s] > 0.001)
    n_zero = N - 1 - n_neg - n_pos
    all_pos = n_neg == 0

    print(f"  p={p}: {n_pos} positive, {n_neg} negative, {n_zero} zero "
          f"{'★ ALL POSITIVE ★' if all_pos else ''}")

    if not all_pos:
        # Which degree has negatives?
        for deg in range(2, m+1, 2):
            neg_at_deg = sum(1 for s in range(1,N) if bin(s).count('1') == deg and h_hat[s] < -0.001)
            if neg_at_deg > 0:
                print(f"    Degree {deg}: {neg_at_deg} negative coefficients")

print("\n" + "=" * 72)
print("CONCLUSIONS")
print("=" * 72)
print("""
PROOF STRATEGY SUMMARY:

1. H has ZERO odd-degree Walsh coefficients (complement symmetry)
2. H has ZERO degree-0 variation (it's the mean)
3. The Walsh spectrum lives in degrees 2, 4, 6, ...

If ALL nonzero Walsh coefficients are POSITIVE for p ≥ p_0:
  - Then H(all-ones) = Σ ĥ[S] ≥ Σ ĥ[S] χ_S(σ) = H(σ) for all σ
  - Since χ_S(σ) ∈ {-1, +1} and ĥ[S] > 0, the sum is maximized
    when ALL χ_S(σ) = +1, i.e., σ = all-ones
  - This is Interval!

Even if not ALL are positive, we need:
  Σ_{S: ĥ[S] > 0} ĥ[S] + Σ_{S: ĥ[S] < 0} ĥ[S] · χ_S(σ) ≥ 0
  for all σ. The positive terms contribute their full weight at all-ones,
  while negative terms contribute ĥ[S]·χ_S(σ) which could be +|ĥ[S]|
  at σ ≠ all-ones. So we need the positive sum to dominate.
""")

print("\nDONE.")
