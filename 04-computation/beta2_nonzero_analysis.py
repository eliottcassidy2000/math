#!/usr/bin/env python3
"""
beta2_nonzero_analysis.py - Characterize when beta_2 > 0 for circulant digraphs

Key finding from tang_yau_counterexample.py: many non-tournament circulant
digraphs with S cap (-S) = empty have beta_2 > 0. Tournament completeness
is essential for beta_2 = 0.

This script analyzes the structure: WHEN does beta_2 > 0 for |S|=2?

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, os
import numpy as np
from itertools import combinations
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
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for s in S:
            A[i][(i + s) % n] = 1
    return A

def compute_betti_quick(A, n, max_dim=3):
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
            S_v = np.linalg.svd(bd_omega[p], compute_uv=False)
            ker = dim_p - sum(s > 1e-8 for s in S_v)
        if p+1 not in bd_omega:
            im = 0
        else:
            S_v = np.linalg.svd(bd_omega[p+1], compute_uv=False)
            im = sum(s > 1e-8 for s in S_v)
        betti.append(ker - im)
    return betti

def neg_S(S, n):
    return frozenset((n - s) % n for s in S)

print("=" * 70)
print("BETA_2 CHARACTERIZATION FOR |S|=2 CIRCULANT DIGRAPHS")
print("=" * 70)

for n in [5, 7, 8, 9, 11, 13]:
    print(f"\n--- n = {n} ---")
    tested = set()
    b2_zero = []
    b2_pos = []

    for S_tuple in combinations(range(1, n), 2):
        S = frozenset(S_tuple)
        nS = neg_S(S, n)
        if S & nS:
            continue
        canon = min(tuple(sorted(S)), tuple(sorted(nS)))
        if canon in tested:
            continue
        tested.add(canon)

        A = build_circulant_digraph(n, S)
        betti = compute_betti_quick(A, n, max_dim=3)
        b2 = betti[2] if len(betti) > 2 else 0

        s1, s2 = sorted(S)
        # Check: is one element double the other mod n?
        double_closed = ((2*s1) % n == s2) or ((2*s2) % n == s1)
        # Check: s1 + s2 = n?
        sum_n = (s1 + s2 == n)
        # Check: s2 - s1 divides n?
        diff_div = (n % (s2 - s1) == 0) if s2 > s1 else False

        info = f"S={{{s1},{s2}}}: b2={b2}, 2*closed={double_closed}, sum_n={sum_n}, diff|n={diff_div}"

        if b2 == 0:
            b2_zero.append((s1, s2, info))
        else:
            b2_pos.append((s1, s2, b2, info))

    print(f"  beta_2 = 0 ({len(b2_zero)} cases):")
    for s1, s2, info in b2_zero:
        print(f"    {info}")
    print(f"  beta_2 > 0 ({len(b2_pos)} cases):")
    for s1, s2, b2, info in b2_pos:
        print(f"    {info}")

    # Summary pattern
    print(f"\n  Pattern check:")
    zero_double = sum(1 for s1,s2,_ in b2_zero if (2*s1)%n==s2 or (2*s2)%n==s1)
    pos_double = sum(1 for s1,s2,_,_ in b2_pos if (2*s1)%n==s2 or (2*s2)%n==s1)
    print(f"    2*-closure in beta_2=0: {zero_double}/{len(b2_zero)}")
    print(f"    2*-closure in beta_2>0: {pos_double}/{len(b2_pos)}")

    # Check consecutive
    zero_consec = sum(1 for s1,s2,_ in b2_zero if s2-s1==1 or s1+1==s2)
    pos_consec = sum(1 for s1,s2,_,_ in b2_pos if s2-s1==1 or s1+1==s2)
    print(f"    Consecutive in beta_2=0: {zero_consec}/{len(b2_zero)}")
    print(f"    Consecutive in beta_2>0: {pos_consec}/{len(b2_pos)}")

print("\nDone.")
