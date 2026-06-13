#!/usr/bin/env python3
"""
beta2_threshold_analysis.py - Find the exact |S| threshold for beta_2=0

Key finding from beta2_doubling_closure_general.py:
- |S|=2: beta_2=0 iff has-doubling-pair
- |S|=3: mixed (some beta_2>0)
- |S|>=4: ALL beta_2=0 at n=9,11

Question: What is the exact threshold? Does it depend on n?

Also analyzes the exceptions: S=(1,4,7) at n=9 has beta_2=0 without doubling pair.

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

def compute_betti_2(A, n):
    allowed = {}
    for p in range(5):
        allowed[p] = enumerate_allowed_paths(A, n, p)
        if not allowed[p]:
            break
    if 2 not in allowed or not allowed[2]:
        return 0
    omega_basis = {}
    for p in range(4):
        if p not in allowed or not allowed[p]:
            omega_basis[p] = np.zeros((0, 0))
            continue
        basis = compute_omega_basis(A, n, p, allowed[p],
                                     allowed[p-1] if p-1 in allowed and p >= 1 else [])
        omega_basis[p] = basis
    bd_omega = {}
    for p in [2, 3]:
        dim_p = omega_basis[p].shape[1] if omega_basis[p].ndim == 2 else 0
        if dim_p == 0:
            continue
        bd = build_full_boundary_matrix(allowed[p], allowed[p-1] if p-1 in allowed else [])
        bd_omega[p] = bd @ omega_basis[p]
    dim_2 = omega_basis[2].shape[1] if omega_basis[2].ndim == 2 else 0
    if dim_2 == 0:
        return 0
    if 2 not in bd_omega:
        ker = dim_2
    else:
        S_v = np.linalg.svd(bd_omega[2], compute_uv=False)
        ker = dim_2 - sum(s > 1e-8 for s in S_v)
    if 3 not in bd_omega:
        im = 0
    else:
        S_v = np.linalg.svd(bd_omega[3], compute_uv=False)
        im = sum(s > 1e-8 for s in S_v)
    return ker - im

def neg_S(S, n):
    return frozenset((n - s) % n for s in S)

def has_doubling_pair(S, n):
    for s in S:
        if (2*s) % n in S:
            return True
    return False

def gcd(a, b):
    while b:
        a, b = b, a % b
    return a

def is_arithmetic_prog(S, n):
    """Check if S is an arithmetic progression mod n."""
    S_sorted = sorted(S)
    if len(S_sorted) < 2:
        return True
    d = (S_sorted[1] - S_sorted[0]) % n
    for i in range(2, len(S_sorted)):
        if (S_sorted[i] - S_sorted[i-1]) % n != d:
            return False
    return True

def S_generates_Zn(S, n):
    """Check if S generates Z_n under addition."""
    g = 0
    for s in S:
        g = gcd(g, s)
    return gcd(g, n) == 1


print("=" * 70)
print("BETA_2 THRESHOLD ANALYSIS")
print("=" * 70)

# Focus on n=13 to see if |S|>=4 still gives all beta_2=0
for n in [13, 15]:
    print(f"\n{'='*70}")
    print(f"n = {n}")
    print(f"{'='*70}")

    for s_size in range(2, min(n//2, 6)):
        print(f"\n  --- |S| = {s_size} ---")
        tested = set()
        b2_zero = 0
        b2_pos = 0
        b2_pos_examples = []
        b2_zero_no_dp = []

        count = 0
        for S_tuple in combinations(range(1, n), s_size):
            S = frozenset(S_tuple)
            nS = neg_S(S, n)
            if S & nS:
                continue
            canon = min(tuple(sorted(S)), tuple(sorted(nS)))
            if canon in tested:
                continue
            tested.add(canon)
            count += 1

            # Limit computation for large cases
            if count > 200:
                break

            A = build_circulant_digraph(n, S)
            b2 = compute_betti_2(A, n)

            dp = has_doubling_pair(S, n)
            ap = is_arithmetic_prog(sorted(S), n)
            gen = S_generates_Zn(S, n)

            if b2 == 0:
                b2_zero += 1
                if not dp:
                    b2_zero_no_dp.append((tuple(sorted(S)), ap, gen))
            else:
                b2_pos += 1
                if len(b2_pos_examples) < 5:
                    b2_pos_examples.append((tuple(sorted(S)), b2, dp, ap, gen))

        print(f"  Tested: {count}, beta_2=0: {b2_zero}, beta_2>0: {b2_pos}")
        if b2_pos_examples:
            print(f"  Examples with beta_2 > 0:")
            for S_t, b2, dp, ap, gen in b2_pos_examples:
                print(f"    S={S_t}: b2={b2}, dp={dp}, AP={ap}, gen={gen}")
        if b2_zero_no_dp:
            print(f"  beta_2=0 WITHOUT doubling pair ({len(b2_zero_no_dp)} cases):")
            for S_t, ap, gen in b2_zero_no_dp[:10]:
                print(f"    S={S_t}: AP={ap}, gen={gen}")

# Specific analysis of S=(1,4,7) at n=9
print(f"\n{'='*70}")
print(f"ANALYSIS OF S=(1,4,7) AT n=9")
print(f"{'='*70}")
n = 9
S = {1, 4, 7}
print(f"  S = {sorted(S)}")
print(f"  -S = {sorted(neg_S(S, n))}")
print(f"  S ∩ (-S) = {sorted(S & neg_S(S, n))}")
print(f"  Has doubling pair: {has_doubling_pair(S, n)}")
print(f"  Is arithmetic progression: {is_arithmetic_prog(sorted(S), n)}")
print(f"  Common diff: {(4-1) % 9} = {3} = n/3")
print(f"  Generates Z_9: {S_generates_Zn(S, n)}")
# Check: 2*1=2 not in S, 2*4=8 not in S, 2*7=14%9=5 not in S
print(f"  Doubles: 2*1={2%9}, 2*4={8%9}, 2*7={14%9} — none in S")
# S is a coset of the subgroup {0,3,6} of Z_9
print(f"  S = 1 + {{0, 3, 6}} = coset of subgroup <3> in Z_9")
print(f"  Subgroup <3> = {{0, 3, 6}}, order 3")
print(f"  This is a COSET structure: S = 1 + H where H is a subgroup")

# Check: does coset structure always give beta_2=0?
print(f"\n  Coset structure check at n=9:")
for a in range(1, 9):
    S_coset = frozenset((a + h) % 9 for h in [0, 3, 6])
    if 0 in S_coset:
        continue
    nS = neg_S(S_coset, 9)
    if S_coset & nS:
        continue
    A_c = build_circulant_digraph(9, S_coset)
    b2_c = compute_betti_2(A_c, 9)
    dp_c = has_doubling_pair(S_coset, 9)
    print(f"    a={a}: S={sorted(S_coset)}, -S={sorted(nS)}, b2={b2_c}, dp={dp_c}")

print("\nDone.")
