#!/usr/bin/env python3
"""
FULL SEESAW ANALYSIS: Compute complete chain complex data for specific tournaments.

For each tournament: compute dim(Ω_p), rank(∂_p), β_p for all p.
Check: is dim(Ω_p) palindromic? Is β₁·β₃ = 0? What determines β₄?

kind-pasteur-2026-03-25-S26
"""
import sys, time, random
import numpy as np
from itertools import permutations
from collections import Counter

sys.stdout.reconfigure(line_buffering=True)
random.seed(314)

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

def transitive_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            A[i][j] = 1
    return A

def paley_tournament(p):
    """Paley tournament on p vertices (p prime, p ≡ 3 mod 4)."""
    qr = set()
    for i in range(1, p):
        qr.add((i*i) % p)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and ((j - i) % p) in qr:
                A[i][j] = 1
    return A

def enum_allowed(A, n, p):
    paths = []
    for perm in permutations(range(n), p+1):
        ok = True
        for i in range(p):
            if A[perm[i]][perm[i+1]] != 1:
                ok = False; break
        if ok: paths.append(perm)
    return paths

def compute_omega_dim(A, n, p):
    ap = enum_allowed(A, n, p)
    if not ap: return 0, np.eye(0)
    if p == 0: return len(ap), np.eye(len(ap))

    ap_m1 = enum_allowed(A, n, p-1)
    ap_m1_set = set(ap_m1)

    non_allowed = {}; na_idx = 0
    for path in ap:
        for i in range(p+1):
            face = path[:i] + path[i+1:]
            if face not in ap_m1_set and face not in non_allowed:
                non_allowed[face] = na_idx; na_idx += 1

    if not non_allowed:
        return len(ap), np.eye(len(ap))

    P = np.zeros((len(non_allowed), len(ap)))
    for j, path in enumerate(ap):
        for i in range(p+1):
            face = path[:i] + path[i+1:]
            sign = (-1)**i
            if face in non_allowed:
                P[non_allowed[face], j] += sign

    _, S, Vt = np.linalg.svd(P)
    rk = sum(s > 1e-8 for s in S)
    dim = len(ap) - rk
    return dim, Vt[rk:].T if dim > 0 else np.zeros((len(ap), 0))

def compute_betti_full(A, n, max_dim=None):
    """Compute ALL Betti numbers and chain complex dimensions."""
    if max_dim is None:
        max_dim = n - 1

    # Compute Ω dimensions and bases
    omega_dims = []
    omega_bases = []
    allowed_paths = []

    for p in range(max_dim + 2):
        ap = enum_allowed(A, n, p)
        allowed_paths.append(ap)
        dim, basis = compute_omega_dim(A, n, p)
        omega_dims.append(dim)
        omega_bases.append(basis)

    # Compute boundary ranks and Betti numbers
    betti = []
    ranks = []

    for p in range(max_dim + 1):
        # rank of ∂_p on Ω_p
        if p == 0 or omega_dims[p] == 0:
            ranks.append(0)
        else:
            ap = allowed_paths[p]
            ap_m1 = allowed_paths[p-1]
            idx = {path: i for i, path in enumerate(ap_m1)}
            D = np.zeros((len(ap_m1), len(ap)))
            for j, path in enumerate(ap):
                for i in range(p+1):
                    face = path[:i] + path[i+1:]
                    sign = (-1)**i
                    if face in idx:
                        D[idx[face], j] += sign
            # Restrict to Ω_p
            D_omega = D @ omega_bases[p]
            _, S, _ = np.linalg.svd(D_omega)
            rk = sum(s > 1e-8 for s in S)
            ranks.append(rk)

        # rank of ∂_{p+1} on Ω_{p+1} (image into Ω_p)
        if p + 1 > max_dim or omega_dims[p+1] == 0:
            im_rank = 0
        else:
            ap1 = allowed_paths[p+1]
            ap = allowed_paths[p]
            idx = {path: i for i, path in enumerate(ap)}
            D1 = np.zeros((len(ap), len(ap1)))
            for j, path in enumerate(ap1):
                for i in range(p+2):
                    face = path[:i] + path[i+1:]
                    sign = (-1)**i
                    if face in idx:
                        D1[idx[face], j] += sign
            D1_omega = D1 @ omega_bases[p+1]
            # Project into Ω_p coordinates
            im_mat = omega_bases[p].T @ D1_omega if omega_dims[p] > 0 else np.zeros((0, D1_omega.shape[1]))
            _, S1, _ = np.linalg.svd(im_mat)
            im_rank = sum(s > 1e-8 for s in S1)

        ker_rank = omega_dims[p] - ranks[p]
        beta = max(0, ker_rank - im_rank)
        betti.append(beta)

    return omega_dims[:max_dim+1], ranks, betti

# ============================================================
# MAIN: Analyze specific tournaments
# ============================================================

print("=" * 80)
print("  FULL CHAIN COMPLEX ANALYSIS")
print("  kind-pasteur-2026-03-25-S26")
print("=" * 80)

n = 7
print(f"\n  n = {n}")
print(f"  dim(ker ∂₁) = C({n},2) - ({n}-1) = {n*(n-1)//2 - n + 1}")

# Specific tournaments
tournaments = {
    "Transitive": transitive_tournament(n),
    "Paley T_7": paley_tournament(n),
}
# Add some random ones
for i in range(5):
    tournaments[f"Random_{i}"] = random_tournament(n)

for name, A in tournaments.items():
    print(f"\n  --- {name} ---")
    t0 = time.time()
    omega_dims, ranks, betti = compute_betti_full(A, n, max_dim=5)
    elapsed = time.time() - t0

    print(f"  dim(Ω_p): {omega_dims}")
    print(f"  rank(∂_p): {ranks}")
    print(f"  β_p:       {betti}")

    # Check palindromy of Ω dimensions
    is_palindrome = all(omega_dims[i] == omega_dims[-(i+1)] for i in range(len(omega_dims)//2))
    print(f"  Ω palindromic: {is_palindrome}")

    # S = dim(ker ∂₁) - dim(Ω₂) + dim(Ω₃)
    ker_d1 = n*(n-1)//2 - n + 1
    S = ker_d1 - omega_dims[2] + omega_dims[3] if len(omega_dims) > 3 else 0
    print(f"  S = {ker_d1} - {omega_dims[2]} + {omega_dims[3] if len(omega_dims) > 3 else '?'} = {S}")
    print(f"  β₁+β₃ = {betti[1]+betti[3] if len(betti) > 3 else '?'}")
    print(f"  β₁·β₃ = {betti[1]*betti[3] if len(betti) > 3 else '?'}")
    print(f"  seesaw: {'HOLDS' if len(betti) > 3 and betti[1]*betti[3] == 0 else 'UNKNOWN'}")

    # Euler characteristic
    chi = sum((-1)**p * d for p, d in enumerate(omega_dims))
    chi_betti = sum((-1)**p * b for p, b in enumerate(betti))
    print(f"  χ = {chi} (from Ω), {chi_betti} (from β)")
    print(f"  [{elapsed:.1f}s]")

print(f"""
  KEY QUESTIONS:
  1. Is Ω palindromic for all tournaments at n=7? (Not just Paley)
  2. Does S grow with tournament complexity?
  3. What determines β₄ > 0?
  4. Is there a pattern in the β profiles?
""")
