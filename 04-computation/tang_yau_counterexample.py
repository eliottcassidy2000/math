#!/usr/bin/env python3
"""
tang_yau_counterexample.py - Test Tang-Yau Conjecture 4.8

Tang-Yau (arXiv:2602.04140) conjecture:
  For circulant digraph C_n^S with S cap (-S) = empty (no-wrap-around),
  H_m(C_n^S) = 0 for all m >= 3.

We test this conjecture for:
1. Small |S| (their regime) - should hold
2. Tournament-size |S| = (n-1)/2 - Paley P_7 is a counterexample!

Also compute beta_2 for various circulant digraphs to check if
beta_2 = 0 extends beyond tournaments.

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, os
import numpy as np
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved

def build_circulant_digraph(n, S):
    """Build adjacency matrix for circulant digraph C_n^S."""
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for s in S:
            j = (i + s) % n
            if j != i:
                A[i][j] = 1
    return A

def compute_betti(A, n, max_dim=None):
    """Compute all Betti numbers up to max_dim."""
    if max_dim is None:
        max_dim = n - 1

    allowed = {}
    for p in range(-1, max_dim + 2):
        if p < 0:
            allowed[p] = []
        else:
            allowed[p] = enumerate_allowed_paths(A, n, p)
            if not allowed[p]:
                max_dim = min(max_dim, p - 1)
                break

    omega_basis = {}
    for p in range(max_dim + 2):
        if p not in allowed or not allowed[p]:
            omega_basis[p] = np.zeros((0, 0))
            continue
        basis = compute_omega_basis(A, n, p, allowed[p],
                                     allowed[p-1] if p-1 in allowed else [])
        omega_basis[p] = basis

    bd_omega = {}
    for p in range(1, max_dim + 2):
        dim_p = omega_basis[p].shape[1] if omega_basis[p].ndim == 2 else 0
        if dim_p == 0:
            continue
        bd = build_full_boundary_matrix(allowed[p], allowed[p-1] if p-1 in allowed else [])
        bd_omega[p] = bd @ omega_basis[p]

    betti = []
    for p in range(max_dim + 1):
        dim_p = omega_basis[p].shape[1] if omega_basis[p].ndim == 2 else 0
        if dim_p == 0:
            betti.append(0)
            continue
        if p not in bd_omega:
            ker = dim_p
        else:
            S_vals = np.linalg.svd(bd_omega[p], compute_uv=False)
            ker = dim_p - sum(s > 1e-8 for s in S_vals)
        if p+1 not in bd_omega:
            im = 0
        else:
            S_vals = np.linalg.svd(bd_omega[p+1], compute_uv=False)
            im = sum(s > 1e-8 for s in S_vals)
        betti.append(ker - im)

    return betti


def neg_S(S, n):
    """Compute -S mod n."""
    return frozenset((n - s) % n for s in S)


print("=" * 70)
print("TANG-YAU CONJECTURE 4.8 TEST")
print("=" * 70)
print("Conjecture: H_m = 0 for m >= 3 when S cap (-S) = empty")
print()

# Part 1: Small |S| cases (their regime)
print("=" * 70)
print("PART 1: Small |S| (Tang-Yau's regime)")
print("=" * 70)

for n in range(3, 12):
    for s_size in [1, 2]:
        if s_size >= n:
            continue
        # Generate all S with |S| = s_size and S cap (-S) = empty
        from itertools import combinations
        tested = set()
        for S_tuple in combinations(range(1, n), s_size):
            S = frozenset(S_tuple)
            nS = neg_S(S, n)
            if S & nS:
                continue  # wrap-around
            canon = min(S, nS)  # avoid duplicates from reversal
            if canon in tested:
                continue
            tested.add(canon)

            A = build_circulant_digraph(n, S)
            try:
                betti = compute_betti(A, n, max_dim=min(5, n-1))
            except Exception as e:
                betti = [f"ERR: {e}"]
            has_high = any(b > 0 for b in betti[3:]) if len(betti) > 3 else False
            marker = " *** COUNTEREXAMPLE ***" if has_high else ""
            print(f"  C_{n}^{sorted(S)}: b={betti}{marker}")

# Part 2: Tournament-size |S| = (n-1)/2
print()
print("=" * 70)
print("PART 2: Tournaments (|S| = (n-1)/2)")
print("=" * 70)

# Paley tournaments
paley_sets = {
    3: {1},       # QR(3) = {1}
    5: {1, 4},    # QR(5) = {1, 4}
    7: {1, 2, 4}, # QR(7) = {1, 2, 4}
}

for p, S in sorted(paley_sets.items()):
    nS = neg_S(S, p)
    wrap = S & nS
    A = build_circulant_digraph(p, S)
    t0 = time.time()
    betti = compute_betti(A, p, max_dim=p-1)
    dt = time.time() - t0
    has_high = any(b > 0 for b in betti[3:]) if len(betti) > 3 else False
    marker = " *** COUNTEREXAMPLE ***" if has_high else ""
    print(f"  P_{p} = C_{p}^{sorted(S)}: wrap={bool(wrap)}, b={betti} ({dt:.1f}s){marker}")

# Non-Paley circulant tournaments at n=9
print()
print("  --- Non-Paley circulant tournaments at n=9 ---")
n = 9
from itertools import combinations
count = 0
for S_tuple in combinations(range(1, n), (n-1)//2):
    S = frozenset(S_tuple)
    nS = neg_S(S, n)
    if S & nS:
        continue  # not a tournament (S must be disjoint from -S)
    count += 1
    A = build_circulant_digraph(n, S)
    t0 = time.time()
    # Only compute up to dim 5 (beta_5 is the interesting one)
    betti = compute_betti(A, n, max_dim=5)
    dt = time.time() - t0
    has_high = any(b > 0 for b in betti[3:]) if len(betti) > 3 else False
    marker = " *** COUNTEREXAMPLE ***" if has_high else ""
    print(f"  C_9^{sorted(S)}: b={betti} ({dt:.1f}s){marker}")

print(f"\n  Total circulant tournaments at n=9: {count}")

# Part 3: Check beta_2 for NON-tournament circulant digraphs
print()
print("=" * 70)
print("PART 3: beta_2 for non-tournament circulant digraphs (no-wrap-around)")
print("=" * 70)
print("Do non-tournament circulants with S cap (-S) = empty also have beta_2=0?")
print()

b2_nonzero = 0
b2_total = 0
for n in range(3, 10):
    for s_size in range(1, n):
        if s_size == (n-1)//2 and n % 2 == 1:
            continue  # skip tournaments (already tested)
        for S_tuple in combinations(range(1, n), s_size):
            S = frozenset(S_tuple)
            nS = neg_S(S, n)
            if S & nS:
                continue
            A = build_circulant_digraph(n, S)
            betti = compute_betti(A, n, max_dim=min(4, n-1))
            b2_total += 1
            if len(betti) > 2 and betti[2] > 0:
                b2_nonzero += 1
                print(f"  beta_2 > 0: C_{n}^{sorted(S)}, b={betti}")

print(f"\n  Tested {b2_total} non-tournament circulant digraphs (no-wrap-around)")
print(f"  beta_2 > 0 in {b2_nonzero} cases")

print(f"\nDone.")
