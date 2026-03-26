#!/usr/bin/env python3
"""
SEESAW ALGEBRAIC PROOF ATTEMPT: Compute the invariant S that controls β₁+β₃.

From the chain complex with β₂=0 (proved THM-108):
  β₁ + β₃ = S - dim(im ∂₄)
  where S = dim(ker ∂₁) - dim(Ω₂) + dim(Ω₃)

If S ≤ 1 for all tournaments, then β₁+β₃ ≤ 1 → seesaw β₁·β₃=0.

This script computes S for random tournaments at n=6,7 to check.

kind-pasteur-2026-03-25-S26
"""
import sys, time, random
import numpy as np
from itertools import permutations
from collections import Counter

sys.stdout.reconfigure(line_buffering=True)
random.seed(2026)

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

def enum_allowed(A, n, p):
    """Enumerate allowed p-paths (directed paths of length p)."""
    paths = []
    for perm in permutations(range(n), p+1):
        ok = True
        for i in range(p):
            if A[perm[i]][perm[i+1]] != 1:
                ok = False; break
        if ok: paths.append(perm)
    return paths

def compute_omega_dim(A, n, p):
    """Compute dim(Ω_p) = dimension of ∂-invariant subspace of A_p."""
    ap = enum_allowed(A, n, p)
    if not ap: return 0

    if p == 0: return len(ap)  # Ω₀ = A₀ (all vertices)

    ap_m1 = enum_allowed(A, n, p-1)
    ap_m1_set = set(ap_m1)

    # Find non-allowed faces
    non_allowed = {}
    na_idx = 0
    for path in ap:
        for i in range(p+1):
            face = path[:i] + path[i+1:]
            if face not in ap_m1_set and face not in non_allowed:
                non_allowed[face] = na_idx
                na_idx += 1

    if not non_allowed:
        return len(ap)  # All faces allowed → Ω_p = A_p

    # Build non-allowed face projection matrix
    P = np.zeros((len(non_allowed), len(ap)))
    for j, path in enumerate(ap):
        for i in range(p+1):
            face = path[:i] + path[i+1:]
            sign = (-1)**i
            if face in non_allowed:
                P[non_allowed[face], j] += sign

    # Ω_p = null space of P
    _, S, _ = np.linalg.svd(P)
    rank_P = sum(s > 1e-8 for s in S)
    return len(ap) - rank_P

def compute_boundary_rank(A, n, p, omega_p_dim=None):
    """Compute rank of ∂_p: Ω_p → Ω_{p-1}."""
    if p == 0: return 0

    ap = enum_allowed(A, n, p)
    ap_m1 = enum_allowed(A, n, p-1)
    if not ap or not ap_m1: return 0

    # Build boundary matrix ∂_p: A_p → A_{p-1}
    idx_m1 = {path: i for i, path in enumerate(ap_m1)}
    D = np.zeros((len(ap_m1), len(ap)))
    for j, path in enumerate(ap):
        for i in range(p+1):
            face = path[:i] + path[i+1:]
            sign = (-1)**i
            if face in idx_m1:
                D[idx_m1[face], j] += sign

    # For exact rank on Ω_p, we'd need to restrict D to Ω_p.
    # But for ker(∂₁), we can use the full D since Ω₁ = A₁ for tournaments
    # (every 1-path = directed arc is trivially ∂-invariant: ∂(a→b) = b - a ∈ A₀)

    _, S, _ = np.linalg.svd(D)
    return sum(s > 1e-8 for s in S)

# ============================================================
# MAIN COMPUTATION
# ============================================================

print("=" * 80)
print("  SEESAW ALGEBRAIC ANALYSIS: Computing S = dim(ker ∂₁) - dim(Ω₂) + dim(Ω₃)")
print("  kind-pasteur-2026-03-25-S26")
print("=" * 80)

for n in [6, 7]:
    print(f"\n--- n = {n} ---")

    # dim(ker ∂₁) is always C(n,2) - (n-1) for connected tournaments
    ker_d1 = n*(n-1)//2 - (n-1)
    print(f"  dim(ker ∂₁) = C({n},2) - ({n}-1) = {ker_d1}")

    n_samples = 50 if n <= 6 else 20
    S_values = []
    beta_pairs = []

    t0 = time.time()
    for trial in range(n_samples):
        A = random_tournament(n)

        # Compute Ω dimensions
        omega2 = compute_omega_dim(A, n, 2)
        omega3 = compute_omega_dim(A, n, 3)

        S = ker_d1 - omega2 + omega3
        S_values.append(S)

        # Also compute β₁ and β₃ to verify
        # β₁ = dim(ker ∂₁) - rank(∂₂|_{Ω₂})
        # For tournaments: Ω₁ = A₁ (all directed arcs are in Ω₁)
        # rank(∂₁) = n-1 (incidence matrix rank)
        # So dim(ker ∂₁) = |A₁| - rank(∂₁) = C(n,2) - (n-1) = ker_d1

        # To compute β₁ properly, need rank(∂₂ restricted to Ω₂)
        # Approximation: compute rank of full ∂₂ (A₂ → A₁) matrix
        rank_d2_full = compute_boundary_rank(A, n, 2)
        # This is an upper bound on rank(∂₂|_{Ω₂})
        # β₁ ≥ ker_d1 - rank_d2_full (but this can be negative, meaning β₁ = 0)

        if (trial + 1) % 10 == 0:
            print(f"  [{trial+1}/{n_samples}] {time.time()-t0:.1f}s")

    S_counter = Counter(S_values)
    print(f"\n  S distribution: {dict(sorted(S_counter.items()))}")
    print(f"  S range: [{min(S_values)}, {max(S_values)}]")
    print(f"  S mean: {sum(S_values)/len(S_values):.2f}")

    if max(S_values) <= 1:
        print(f"  *** S ≤ 1 for all samples → SEESAW PROVED (if this holds for all tournaments)")
    else:
        print(f"  S > 1 observed → seesaw requires dim(im ∂₄) ≥ {max(S_values) - 1}")

print(f"""
  INTERPRETATION:
  If S = dim(ker ∂₁) - dim(Ω₂) + dim(Ω₃) ≤ 1 for ALL tournaments,
  then β₁ + β₃ = S - dim(im ∂₄) ≤ 1 (since dim(im ∂₄) ≥ 0),
  which gives β₁·β₃ = 0 (the seesaw).

  If S > 1 for some tournaments, the seesaw requires dim(im ∂₄) ≥ S-1.
  This would need a separate argument about ∂₄.
""")
