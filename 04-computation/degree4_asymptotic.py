#!/usr/bin/env python3
"""
Asymptotic analysis of degree-4 Walsh coefficients.

The key quantity is ĥ[{i,j,k,l}] = (1/2^m) Σ_σ H(σ) σ_i σ_j σ_k σ_l.

Since H(σ) depends on σ only through the imaginary parts b_t = Σ_j σ_j sin(2πjt/p),
and b_t is LINEAR in σ, the degree-4 Walsh coefficient measures the 4th-order
NONLINEAR coupling between chords i,j,k,l through the Hamiltonian path count.

APPROACH: Express ĥ[S] via the MULTILINEAR EXTENSION of H.

H(σ) = Σ_T c_T Π_{j∈T} σ_j (Walsh expansion)
ĥ[S] = c_S (the Walsh coefficient)

We want to understand c_S for |S|=4 in terms of the eigenvalue structure.

KEY IDEA: Use the SECOND DERIVATIVE test.

∂⁴H / (∂σ_i ∂σ_j ∂σ_k ∂σ_l) evaluated at σ = 1_m (all-ones)
gives exactly 4! · ĥ[{i,j,k,l}] (from the multilinear extension).

Wait — H is already multilinear (function of ±1 variables), so:
  ĥ[S] = (1/2^m) Σ_σ H(σ) χ_S(σ)

This is just the standard Walsh-Hadamard coefficient. There's no derivative
to take since σ is discrete.

But we CAN relate it to the CONTINUOUS eigenvalue landscape.

Define H̃(b) = H-value for the tournament whose eigenvalue imaginary parts are b.
Then: H̃(Bσ) = H(σ) where B is the sine matrix.

The degree-4 Walsh coefficient then involves the 4th-order Taylor coefficient
of H̃ around the all-ones eigenvalue point.

opus-2026-03-12-S67
"""

import numpy as np
from itertools import combinations
from sympy.ntheory import legendre_symbol as legendre
import time

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
print("PART I: THE B-MATRIX ROTATION VIEW")
print("=" * 72)
print("""
The sine matrix B[t,j] = sin(2πjt/p) maps σ ∈ {±1}^m to b ∈ R^{p-1}.
Since B has all singular values = √(p/4), the mapping σ → b = Bσ
is a SCALED ROTATION.

The degree-4 Walsh coefficient ĥ[{i,j,k,l}] can be written as:
  ĥ[S] = (1/2^m) Σ_σ H̃(Bσ) · σ_i σ_j σ_k σ_l

where H̃ is the function that maps eigenvalue imaginary parts to H.

In the B-coordinates:
  σ_i = Σ_t B†[i,t] b_t / (p/4)   (pseudo-inverse)

So the product σ_i σ_j σ_k σ_l involves a sum over 4 eigenvalue indices.

The CONTINUOUS approximation (valid for large m):
  ĥ[{i,j,k,l}] ≈ (4/p)^4 Σ_{t₁,t₂,t₃,t₄} B[t₁,i]B[t₂,j]B[t₃,k]B[t₄,l]
                    × (∂⁴H̃/∂b_{t₁}∂b_{t₂}∂b_{t₃}∂b_{t₄})

This connects the Walsh coefficient to the 4th-order sensitivity of H
to eigenvalue perturbations!
""")

for p in [13]:
    m = (p - 1) // 2

    # Sine matrix
    B = np.zeros((p-1, m))
    for t in range(1, p):
        for j in range(1, m+1):
            B[t-1, j-1] = np.sin(2 * np.pi * j * t / p)

    # SVD of B
    U, S_vals, Vt = np.linalg.svd(B, full_matrices=False)
    print(f"  p={p}: B is {p-1} x {m}")
    print(f"  Singular values: {[f'{s:.4f}' for s in S_vals]}")
    print(f"  Expected: √(p/4) = {np.sqrt(p/4):.4f}")
    print(f"  All equal? {np.allclose(S_vals, np.sqrt(p/4))}")

    # B^T B = (p/4) I_m (Parseval)
    BtB = B.T @ B
    print(f"  B^T B / (p/4) ≈ I? {np.allclose(BtB, (p/4)*np.eye(m))}")

# ========================================================================
print("\n" + "=" * 72)
print("PART II: THE 4-POINT EIGENVALUE CORRELATOR")
print("=" * 72)
print("""
The 4-point correlator:
  C₄(i,j,k,l) = Σ_t sin(2πit/p) sin(2πjt/p) sin(2πkt/p) sin(2πlt/p)

This is the "spectral 4-point function" — it measures how much
the 4 sine waves (corresponding to chords i,j,k,l) overlap in
the eigenvalue spectrum.

Using product-to-sum formulas:
  4 sin A sin B sin C sin D =
  cos(A-B-C+D) + cos(A-B+C-D) + cos(A+B-C-D)
  - cos(A+B+C+D) - cos(A-B-C-D) - cos(A+B-C+D) - cos(A-B+C+D)
  + cos(A+B+C-D)

Actually, let's just use:
  sin a · sin b = [cos(a-b) - cos(a+b)]/2

Applied twice:
  sin a sin b sin c sin d = (1/4)[cos(a-b) - cos(a+b)][cos(c-d) - cos(c+d)]

So C₄ involves sums of cos(2π(±i±j±k±l)t/p).

By orthogonality: Σ_t cos(2πnt/p) = p if n≡0 mod p, else 0.

So C₄ = (p/4) · #{±i±j±k±l ≡ 0 mod p} / 4.

Wait, let me be more careful.
""")

for p in [13]:
    m = (p - 1) // 2

    print(f"\n  p={p}: 4-point correlator analysis")

    # Compute C₄ for all 4-subsets
    B = np.zeros((p-1, m))
    for t in range(1, p):
        for j in range(1, m+1):
            B[t-1, j-1] = np.sin(2 * np.pi * j * t / p)

    for S_tuple in [(1,2,3,4), (1,2,5,6), (1,3,4,5), (2,4,5,6)]:
        i, j, k, l = [s-1 for s in S_tuple]
        C4 = sum(B[t,i]*B[t,j]*B[t,k]*B[t,l] for t in range(p-1))

        # Analytic: count how many of ±i±j±k±l ≡ 0 mod p
        # where the indices are the chord numbers (1-indexed)
        a, b, c, d = S_tuple
        zero_count = 0
        for s1 in [1, -1]:
            for s2 in [1, -1]:
                for s3 in [1, -1]:
                    for s4 in [1, -1]:
                        if (s1*a + s2*b + s3*c + s4*d) % p == 0:
                            zero_count += 1

        C4_analytic = (p-1) * zero_count / 16  # approximate
        print(f"    C₄{S_tuple}: numeric={C4:.4f}, "
              f"#(±i±j±k±l≡0 mod {p})={zero_count}, "
              f"(p-1)·count/16={C4_analytic:.4f}")

    # The exact formula:
    # C₄(i,j,k,l) = (1/4) Σ_{t=1}^{p-1} [cos(2π(i-j)t/p) - cos(2π(i+j)t/p)]
    #                                       [cos(2π(k-l)t/p) - cos(2π(k+l)t/p)]
    # = (1/4) Σ_{e∈{i-j,-(i+j)}} Σ_{f∈{k-l,-(k+l)}} Σ_t cos(2π(e+f)t/p)·sign

    print(f"\n    Exact computation via product-to-sum:")
    for S_tuple in combinations(range(1, m+1), 4):
        a, b, c, d = S_tuple

        # sin(a)sin(b) = [cos(a-b) - cos(a+b)]/2
        # sin(c)sin(d) = [cos(c-d) - cos(c+d)]/2
        # Product = (1/4)[cos(a-b)cos(c-d) - cos(a-b)cos(c+d)
        #                 - cos(a+b)cos(c-d) + cos(a+b)cos(c+d)]
        # Each cos·cos = [cos(sum)+cos(diff)]/2
        # So we get 8 cosine terms

        C4_exact = 0.0
        pairs_ab = [(a-b, 1), (-(a+b), -1)]
        pairs_cd = [(c-d, 1), (-(c+d), -1)]

        for (e, se) in pairs_ab:
            for (f, sf) in pairs_cd:
                sign = se * sf
                # Σ_{t=1}^{p-1} cos(2π(e+f)t/p)
                n = (e + f) % p
                cossum = (p - 1) if n == 0 else -1  # orthogonality
                C4_exact += sign * cossum / 4.0

        # Also need to account for: there are 3 pairings of 4 elements into 2 pairs
        # (ab)(cd), (ac)(bd), (ad)(bc)

        C4_exact_full = 0.0
        for (pair1, pair2) in [((a,b),(c,d)), ((a,c),(b,d)), ((a,d),(b,c))]:
            x1, x2 = pair1
            y1, y2 = pair2
            for (e, se) in [(x1-x2, 1), (-(x1+x2), -1)]:
                for (f, sf) in [(y1-y2, 1), (-(y1+y2), -1)]:
                    sign = se * sf
                    n = (e + f) % p
                    cossum = (p - 1) if n == 0 else -1
                    C4_exact_full += sign * cossum / 4.0

        # Hmm, the product sin(a)sin(b)sin(c)sin(d) uses just ONE pairing.
        # We need to pick the specific pairing (ab)(cd) and then permute.
        # Actually, sin·sin·sin·sin is NOT decomposed by pairings — it's
        # a single product of 4 terms.

        # Let me just compute directly:
        i0, j0, k0, l0 = a-1, b-1, c-1, d-1
        C4_direct = sum(B[t,i0]*B[t,j0]*B[t,k0]*B[t,l0] for t in range(p-1))

        # Zero-sum count
        zc = 0
        for s1 in [1, -1]:
            for s2 in [1, -1]:
                for s3 in [1, -1]:
                    for s4 in [1, -1]:
                        if (s1*a + s2*b + s3*c + s4*d) % p == 0:
                            zc += 1

        if abs(C4_direct) > 0.01 or zc > 0:
            print(f"      {S_tuple}: C₄={C4_direct:>8.4f}, "
                  f"zero-sum={zc}")

# ========================================================================
print("\n" + "=" * 72)
print("PART III: THE ZERO-SUM CHARACTERIZATION")
print("=" * 72)
print("""
From the product-to-sum analysis:
  C₄(i,j,k,l) = (p/16) · W(i,j,k,l) - corrections

where W(i,j,k,l) = #{sign patterns (s₁,s₂,s₃,s₄) ∈ {±1}⁴ : s₁i+s₂j+s₃k+s₄l ≡ 0 mod p}

The zero-sum patterns tell us about the ADDITIVE STRUCTURE of {i,j,k,l} mod p.

W = 0: "generic" quadruple (no additive relation)
W = 2: one additive relation (e.g., i+j ≡ k+l mod p)
W = 4: two relations (arithmetic progression)
W = 8: highly structured (e.g., i=j=k=l)

The SIGN of ĥ[S] should correlate with W: more structure → ???
""")

for p in [13]:
    m = (p - 1) // 2
    N = 2**m

    # Compute H and Walsh spectrum
    H_vec = np.zeros(N, dtype=np.int64)
    for bits in range(N):
        sigma = tuple(1 if (bits >> k) & 1 else -1 for k in range(m))
        S = orientation_to_S(sigma, p)
        A = make_tournament(p, S)
        H_vec[bits] = count_H(A)

    h_hat = np.zeros(N)
    for S_bits in range(N):
        val = 0.0
        for bits in range(N):
            common = bin(S_bits & bits).count('1')
            deg = bin(S_bits).count('1')
            chi = (-1)**(deg - common)
            val += H_vec[bits] * chi
        h_hat[S_bits] = val / N

    print(f"\n  p={p}: Zero-sum classification of degree-4 coefficients:")
    print(f"  {'Subset':>15} {'ĥ':>12} {'W':>4} {'|ĥ| class':>12} {'sign':>5}")

    # Classify ALL degree-4 coefficients by W
    from collections import defaultdict
    w_groups = defaultdict(list)

    for S_tuple in combinations(range(1, m+1), 4):
        S_bits = sum(1 << (k-1) for k in S_tuple)
        coeff = h_hat[S_bits]
        if abs(coeff) < 0.001:
            continue

        a, b, c, d = S_tuple
        W = 0
        for s1 in [1, -1]:
            for s2 in [1, -1]:
                for s3 in [1, -1]:
                    for s4 in [1, -1]:
                        if (s1*a + s2*b + s3*c + s4*d) % p == 0:
                            W += 1

        # Magnitude class
        abs_coeff = abs(coeff)
        if abs_coeff > 5000:
            mag_class = "large"
        elif abs_coeff > 3000:
            mag_class = "medium"
        else:
            mag_class = "small"

        sign = "+" if coeff > 0 else "-"
        print(f"  {str(S_tuple):>15} {coeff:>12.2f} {W:>4} {mag_class:>12} {sign:>5}")

        w_groups[W].append((coeff, S_tuple))

    print(f"\n  Summary by W:")
    for W in sorted(w_groups.keys()):
        entries = w_groups[W]
        pos = sum(1 for c, _ in entries if c > 0)
        neg = len(entries) - pos
        pos_sum = sum(c for c, _ in entries if c > 0)
        neg_sum = sum(c for c, _ in entries if c < 0)
        print(f"    W={W}: {len(entries)} coefficients, "
              f"{pos} pos (sum={pos_sum:.2f}), {neg} neg (sum={neg_sum:.2f}), "
              f"net={pos_sum+neg_sum:.2f}")

# ========================================================================
print("\n" + "=" * 72)
print("PART IV: COMPARISON ACROSS PRIMES")
print("=" * 72)

# For p=7 and p=11, same analysis
for p in [7, 11]:
    m = (p - 1) // 2
    N = 2**m

    H_vec = np.zeros(N, dtype=np.int64)
    for bits in range(N):
        sigma = tuple(1 if (bits >> k) & 1 else -1 for k in range(m))
        S = orientation_to_S(sigma, p)
        A = make_tournament(p, S)
        H_vec[bits] = count_H(A)

    h_hat = np.zeros(N)
    for S_bits in range(N):
        val = 0.0
        for bits in range(N):
            common = bin(S_bits & bits).count('1')
            deg = bin(S_bits).count('1')
            chi = (-1)**(deg - common)
            val += H_vec[bits] * chi
        h_hat[S_bits] = val / N

    # Degree-4 coefficients
    deg4_exists = False
    for S_tuple in combinations(range(1, m+1), 4):
        S_bits = sum(1 << (k-1) for k in S_tuple)
        if abs(h_hat[S_bits]) > 0.001:
            deg4_exists = True

    if not deg4_exists:
        print(f"\n  p={p}: No degree-4 coefficients (m={m} < 4)")
        continue

    print(f"\n  p={p}: Degree-4 coefficients with W:")
    for S_tuple in combinations(range(1, m+1), 4):
        S_bits = sum(1 << (k-1) for k in S_tuple)
        coeff = h_hat[S_bits]
        if abs(coeff) < 0.001:
            continue

        a, b, c, d = S_tuple
        W = 0
        for s1 in [1, -1]:
            for s2 in [1, -1]:
                for s3 in [1, -1]:
                    for s4 in [1, -1]:
                        if (s1*a + s2*b + s3*c + s4*d) % p == 0:
                            W += 1

        sign = "+" if coeff > 0 else "-"
        print(f"    [{sign}] {S_tuple}: ĥ={coeff:>10.4f}, W={W}")

# ========================================================================
print("\n" + "=" * 72)
print("PART V: THE CONNECTION TO ADDITIVE COMBINATORICS")
print("=" * 72)
print("""
SUMMARY:

The degree-4 Walsh coefficient ĥ[{i,j,k,l}] depends on:
  1. The ADDITIVE structure of {i,j,k,l} mod p (via zero-sum count W)
  2. The POSITION of the 4-subset within {1,...,m} (left/right bias)
  3. The FULL nonlinear H function (not just cycle counts)

The zero-sum count W captures how many additive relations exist
among ±i, ±j, ±k, ±l modulo p. This connects to:

  - SUMSET theory (Freiman, Ruzsa): structured sets have more relations
  - ADDITIVE ENERGY: E(X) = #{(a,b,c,d) : a+b=c+d} is related to W
  - EXPONENTIAL SUMS: W involves character sums over F_p

The PROOF STRATEGY:
  1. Express ĥ[{i,j,k,l}] as a sum over eigenvalue 4-tuples
  2. Show the dominant contribution comes from W > 0 quadruples
  3. Prove the hyperplane condition holds using character sum bounds

This connects the tournament H-maximization problem to DEEP questions
in additive combinatorics and analytic number theory.
""")

print("\nDONE.")
