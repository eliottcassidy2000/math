#!/usr/bin/env python3
"""
FOURIER DECOMPOSITION for path homology of circulant digraphs

Implements the Tang-Yau technique (arXiv:2602.04140):
The shift automorphism τ: v → v+1 mod n decomposes the chain complex
into eigenspaces indexed by n-th roots of unity λ.

For each λ with λ^n = 1:
  The boundary map ∂_m restricted to eigenspace λ has a "symbol matrix"
  whose entries are determined by the connection set S.

This gives: β_m = Σ_{λ^n=1} β_m^(λ)

Key advantages:
1. Each eigenspace is much smaller than the full chain complex
2. For large n, most eigenspaces contribute 0
3. Makes dependence on n and S transparent
"""
import numpy as np
from itertools import combinations
from collections import Counter
import sys
sys.path.insert(0, '04-computation')
from path_homology_v2 import path_betti_numbers, circulant_digraph

def fourier_symbol_matrix_1(S, n, lam):
    """
    Symbol matrix M_1(λ) for the boundary ∂_1 on eigenspace λ.

    For C_n^S with S = {s_0, ..., s_{d-1}}, the allowed 1-paths
    in eigenspace λ correspond to generators e_{s_j}.

    M_1(λ) maps from 1-chains to 0-chains.
    ∂(v, v+s) = (v+s) - (v) = λ^s * e_0 - e_0 = (λ^s - 1) * e_0

    So M_1(λ) = row vector [λ^{s_0} - 1, λ^{s_1} - 1, ..., λ^{s_{d-1}} - 1]

    Returns: numpy matrix (m x d) where m = dim of 0-eigenspace = 1
    """
    d = len(S)
    M = np.array([[lam**s - 1 for s in S]])
    return M

def omega_2_generators(S, n):
    """
    Find which 2-paths (v, v+s1, v+s1+s2) are ∂-invariant.

    A 2-path (v, v+s1, v+s1+s2) has boundary:
      (v+s1, v+s1+s2) - (v, v+s1+s2) + (v, v+s1)

    The first and third terms are allowed (they're edges with steps s2, s1).
    The middle term (v, v+s1+s2) is allowed iff (s1+s2) mod n ∈ S.

    So the 2-path with steps (s1, s2) is ∂-invariant iff (s1+s2) mod n ∈ S.
    """
    omega_pairs = []
    for s1 in S:
        for s2 in S:
            if s1 == s2:
                continue  # vertices must be distinct
            total = (s1 + s2) % n
            if total in S and total != s1 and total != 0:
                # Also need v+s1+s2 ≠ v, i.e., s1+s2 ≠ 0 mod n
                # And v+s1+s2 ≠ v+s1, i.e., s2 ≠ 0 (already guaranteed since s2 ∈ S, S ⊂ {1,...,n-1})
                # And v+s1+s2 ≠ v, i.e., total ≠ 0
                omega_pairs.append((s1, s2, total))
    return omega_pairs

def fourier_betti(S, n, max_dim=4):
    """
    Compute Betti numbers via Fourier decomposition.
    This is a simplified version — we compute β_0 and β_1 via symbol matrices,
    and compare with the full computation.
    """
    S_set = set(S)
    roots = [np.exp(2j * np.pi * k / n) for k in range(n)]

    # β_0: contributed by eigenspaces where M_1(λ) has rank 0
    # M_1(λ) = [λ^s - 1 for s in S]
    # This is the zero vector iff λ^s = 1 for all s ∈ S
    # i.e., ord(λ) divides gcd(s for s in S, n)
    beta_0 = 0
    beta_1_fourier = 0

    for k in range(n):
        lam = roots[k]
        M1 = fourier_symbol_matrix_1(S, n, lam)
        rank1 = np.linalg.matrix_rank(M1, tol=1e-8)

        dim_omega_1 = len(S)  # one generator per step s ∈ S
        dim_omega_0 = 1  # one generator (the eigenspace of 0-chains)

        # β_0^(λ) = dim(Ω_0^(λ)) - rank(M_1(λ)) - rank(M_1_from_above) = dim_omega_0 - rank1
        # but Ω_{-1} = 0, so there's no im from above
        beta_0_lam = dim_omega_0 - rank1  # = ker(∂_0) (∂_0 = 0)... actually β_0 = ker(∂_0)/im(∂_1) = dim_omega_0 - rank1
        beta_0 += max(0, beta_0_lam)

        # For β_1: need dim(ker ∂_1) - rank(∂_2)
        # ker(∂_1) on eigenspace λ: dim = dim_omega_1 - rank1... no
        # ker(∂_1) = nullity of M_1(λ) = dim_omega_1 - rank1
        ker_1 = dim_omega_1 - rank1

        # For ∂_2, need Ω_2 generators
        omega_2 = omega_2_generators(S, n)
        dim_omega_2 = len(omega_2)

        if dim_omega_2 > 0 and ker_1 > 0:
            # Build M_2(λ): boundary ∂ on 2-chains in eigenspace λ
            # ∂(v, v+s1, v+s1+s2) = (v+s1, v+s1+s2) - (v, v+s1+s2) + (v, v+s1)
            # In Fourier: the eigenspace generator for step s is e_s with τ(e_s) = λ * e_s
            # Actually, the 2-chain with steps (s1, s2) maps to:
            #   λ^{s1} * e_{s2} - e_{s1+s2} + e_{s1}  (in terms of 1-chain generators)

            # Build the matrix: columns indexed by omega_2 generators, rows by S (1-chain gens)
            S_list = list(S)
            S_idx = {s: i for i, s in enumerate(S_list)}
            M2 = np.zeros((len(S_list), dim_omega_2), dtype=complex)
            for j, (s1, s2, total) in enumerate(omega_2):
                # ∂(steps s1, s2) = λ^{s1} * e_{s2} - e_{total} + e_{s1}
                if s1 in S_idx:
                    M2[S_idx[s1], j] += 1
                if s2 in S_idx:
                    M2[S_idx[s2], j] += lam**s1
                if total in S_idx:
                    M2[S_idx[total], j] -= 1
            rank2 = np.linalg.matrix_rank(M2, tol=1e-8)
        else:
            rank2 = 0

        beta_1_lam = ker_1 - rank2
        beta_1_fourier += max(0, beta_1_lam)

    return beta_0, beta_1_fourier

# ===== TEST: Compare Fourier vs full computation =====
print("=" * 70)
print("FOURIER vs FULL COMPUTATION")
print("=" * 70)

for n in [5, 6, 7, 8, 9, 10]:
    print(f"\nn={n}:")
    for size in range(1, min(n, 4)):
        for S in combinations(range(1, n), size):
            S_list = list(S)
            S_set = set(S_list)

            # Full computation
            A = circulant_digraph(n, S_list)
            betti_full = path_betti_numbers(A, n, max_dim=min(n-1, 4))

            # Fourier computation (β_0 and β_1 only)
            b0_f, b1_f = fourier_betti(S_set, n)

            match = (betti_full[0] == b0_f and betti_full[1] == b1_f)
            if not match:
                print(f"  MISMATCH C_{n}^{S_set}: full β={betti_full}, Fourier β_0={b0_f}, β_1={b1_f}")
            elif any(b > 0 for b in betti_full[2:]):
                print(f"  C_{n}^{S_set}: β={betti_full}, Fourier β_0={b0_f}✓, β_1={b1_f}✓")

# ===== Omega_2 structure analysis =====
print("\n\n" + "=" * 70)
print("Ω_2 STRUCTURE: WHEN DO 2-PATHS SURVIVE?")
print("=" * 70)

print("\nThe key condition: 2-path with steps (s1,s2) is in Ω_2 iff (s1+s2) mod n ∈ S")
print("This defines a 'closure' property of S under addition mod n.\n")

for n in [5, 7, 11]:
    print(f"n={n}:")
    for size in range(1, min(n, 5)):
        for S in combinations(range(1, n), size):
            S_set = set(S)
            omega_2 = omega_2_generators(S_set, n)
            if omega_2:
                A = circulant_digraph(n, list(S))
                betti = path_betti_numbers(A, n, max_dim=min(n-1, 4))
                print(f"  S={S_set}: |Ω_2 generators|={len(omega_2)}, β={betti}")
                if len(omega_2) <= 6:
                    for s1, s2, total in omega_2:
                        print(f"    steps ({s1},{s2}), sum={total}")

# ===== Large n exploration =====
print("\n\n" + "=" * 70)
print("STABILITY: BETTI NUMBERS FOR LARGE n")
print("=" * 70)

print("\nTang-Yau predict: β stabilizes for large primes")
print("Testing S={1,3} across different n:\n")

for n in [5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 19, 23]:
    if 3 < n:
        A = circulant_digraph(n, [1, 3])
        betti = path_betti_numbers(A, n, max_dim=min(n-1, 4))
        print(f"  C_{n}^{{1,3}}: β={betti}")

print("\nTesting S={1,2,4} across different n:\n")
for n in [5, 6, 7, 8, 9, 10, 11, 13, 17]:
    A = circulant_digraph(n, [1, 2, 4])
    betti = path_betti_numbers(A, n, max_dim=min(n-1, 4))
    print(f"  C_{n}^{{1,2,4}}: β={betti}")

print("\n\nDone.")
