#!/usr/bin/env python3
"""
beta2_deletion_implies_global.py — Test the KEY claim:

    β₁(T\v) > 0 for ALL v  ==>  β₁(T) > 0

If this holds, then β₁(T) = 0 implies ∃v with β₁(T\v) = 0,
which gives h₂_rel(v) = 0, hence β₂(T) = 0 by induction.

This is WEAKER than Σ ≤ 3 but SUFFICIENT for the β₂ proof.

Also test the contrapositive: 
    β₁(T) = 0  ==>  ∃v with β₁(T\v) = 0

Verified exhaustive n=5,6; sampled n=7,8,9.

Author: opus-2026-03-08-S49
"""
import sys, time, random
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = __import__('os').fdopen(__import__('os').open(__import__('os').devnull, __import__('os').O_WRONLY), 'w')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved

def dim_om(om):
    return om.shape[1] if om.ndim == 2 and om.shape[0] > 0 else 0

def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A

def build_random_adj(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

def compute_beta1(A, n):
    ap0 = enumerate_allowed_paths(A, n, 0)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
    d1 = dim_om(om1)
    if d1 == 0: return 0
    bd1 = build_full_boundary_matrix(ap1, ap0)
    rk1 = np.linalg.matrix_rank(bd1 @ om1, tol=1e-8)
    d2 = dim_om(om2)
    if d2 > 0:
        bd2 = build_full_boundary_matrix(ap2, ap1)
        bd2om = np.linalg.lstsq(om1, bd2 @ om2, rcond=None)[0]
        b1 = np.linalg.matrix_rank(bd2om, tol=1e-8)
    else:
        b1 = 0
    return (d1 - rk1) - b1


# ============================================================
# Exhaustive test
# ============================================================
print("=" * 70)
print("KEY CLAIM: β₁(T\v)>0 ∀v ==> β₁(T)>0")
print("Contrapositive: β₁(T)=0 ==> ∃v with β₁(T\v)=0")
print("=" * 70)

for n in [4, 5, 6]:
    m = n*(n-1)//2
    total = 1 << m
    t0 = time.time()
    
    claim_holds = True
    count_all_del_positive = 0
    count_b1_zero = 0
    
    for bits in range(total):
        A = build_adj(n, bits)
        beta1_T = compute_beta1(A, n)
        
        all_del_positive = True
        has_zero_del = False
        for v in range(n):
            others = [i for i in range(n) if i != v]
            A_sub = [[A[others[i]][others[j]] for j in range(n-1)] for i in range(n-1)]
            b1v = compute_beta1(A_sub, n-1)
            if b1v == 0:
                all_del_positive = False
                has_zero_del = True
                break
            
        if not has_zero_del:
            # Check remaining vertices
            for v2 in range(v+1, n):
                others2 = [i for i in range(n) if i != v2]
                A_sub2 = [[A[others2[i]][others2[j]] for j in range(n-1)] for i in range(n-1)]
                b1v2 = compute_beta1(A_sub2, n-1)
                if b1v2 == 0:
                    has_zero_del = True
                    break
            all_del_positive = not has_zero_del
        
        if all_del_positive:
            count_all_del_positive += 1
            if beta1_T == 0:
                claim_holds = False
                scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]
                print(f"  COUNTEREXAMPLE at n={n}: bits={bits}, β₁(T)=0, all β₁(T\v)>0, scores={scores}")
        
        if beta1_T == 0:
            count_b1_zero += 1
            if not has_zero_del:
                scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]
                print(f"  VIOLATION at n={n}: β₁(T)=0 but no v with β₁(T\v)=0, scores={scores}")
        
        if (bits+1) % 5000 == 0:
            elapsed = time.time() - t0
            print(f"  n={n}: {bits+1}/{total} ({elapsed:.0f}s)")
    
    elapsed = time.time() - t0
    print(f"\nn={n}: {total} tournaments in {elapsed:.0f}s")
    print(f"  ##all_del_pos = {count_all_del_positive}")
    print(f"  ##b1_zero = {count_b1_zero}")
    print(f"  Claim holds: {'✓' if claim_holds else '✗'}")


# ============================================================
# Sampled test for larger n
# ============================================================
for n, samples in [(7, 5000), (8, 2000), (9, 500)]:
    t0 = time.time()
    claim_holds = True
    count_all_del_positive = 0
    count_b1_zero_no_del = 0
    
    for trial in range(samples):
        A = build_random_adj(n)
        beta1_T = compute_beta1(A, n)
        
        has_zero_del = False
        for v in range(n):
            others = [i for i in range(n) if i != v]
            A_sub = [[A[others[i]][others[j]] for j in range(n-1)] for i in range(n-1)]
            b1v = compute_beta1(A_sub, n-1)
            if b1v == 0:
                has_zero_del = True
                break
        
        if not has_zero_del:
            count_all_del_positive += 1
            if beta1_T == 0:
                claim_holds = False
                scores = sorted([sum(A[i][j] for j in range(n) if j!=i) for i in range(n)])
                print(f"  COUNTEREXAMPLE at n={n}: β₁(T)=0 but all β₁(T\v)>0, scores={scores}")
        
        if beta1_T == 0 and not has_zero_del:
            count_b1_zero_no_del += 1
        
        if (trial+1) % 1000 == 0:
            elapsed = time.time() - t0
            print(f"  n={n}: {trial+1}/{samples} ({elapsed:.0f}s)")
    
    elapsed = time.time() - t0
    print(f"\nn={n}: {samples} tournaments in {elapsed:.0f}s")
    print(f"  ##all_del_pos = {count_all_del_positive}")
    print(f"  ##b1z_no_del = {count_b1_zero_no_del}")
    print(f"  Claim holds: {'✓' if claim_holds else '✗'}")

print("\nDone.")
