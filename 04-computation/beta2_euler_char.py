#!/usr/bin/env python3
"""
beta2_euler_char.py — Euler characteristic and dimension analysis for Ω chain complex

Key question: Does χ(Ω) = 1 for all tournaments? If so, this gives a relation
between Betti numbers: β₁ = β₂ - β₃ + β₄ - ...

If β₂ = 0 universally, then β₁ + β₃ = β₄ + ..., and mutual exclusivity
of β₁ and β₃ would give further constraints.

Author: opus-2026-03-08-S49
"""
import sys, time
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

def compute_omega_dims(A, n, max_p=None):
    """Compute dim(Ω_p) for all p."""
    if max_p is None:
        max_p = n - 1
    dims = [n]  # Ω_0 = vertices
    
    ap_prev = list(range(n))  # A_0 = vertices
    for p in range(1, max_p + 1):
        ap = enumerate_allowed_paths(A, n, p)
        if not ap:
            dims.append(0)
            break
        ap_prev_p = enumerate_allowed_paths(A, n, p-1)
        om = compute_omega_basis(A, n, p, ap, ap_prev_p)
        dims.append(dim_om(om))
        if dim_om(om) == 0:
            break
    
    # Pad with zeros
    while len(dims) <= max_p:
        dims.append(0)
    
    return dims

# ============================================================
# Euler characteristic at n=5
# ============================================================
print("=" * 70)
print("EULER CHARACTERISTIC χ(Ω) FOR TOURNAMENTS")
print("=" * 70)

for n in [4, 5, 6]:
    m = n*(n-1)//2
    total = 1 << m
    if n >= 7:
        break
    
    chi_dist = Counter()
    dim_examples = {}
    t0 = time.time()
    
    for bits in range(total):
        A = build_adj(n, bits)
        dims = compute_omega_dims(A, n)
        chi = sum((-1)**p * d for p, d in enumerate(dims))
        chi_dist[chi] += 1
        
        # Store a few examples
        key = tuple(dims)
        if key not in dim_examples:
            dim_examples[key] = bits
        
        if (bits+1) % 5000 == 0:
            print(f"  n={n}: {bits+1}/{total}")
    
    elapsed = time.time() - t0
    print(f"\nn={n}: {total} tournaments, {elapsed:.0f}s")
    print(f"  χ distribution: {dict(sorted(chi_dist.items()))}")
    print(f"  Ω_* dimension profiles ({len(dim_examples)} distinct):")
    for dims, bits in sorted(dim_examples.items()):
        chi = sum((-1)**p * d for p, d in enumerate(dims))
        print(f"    {dims} (χ={chi}), e.g. bits={bits}")

# ============================================================
# At n=7: sample to check χ and get dim profiles
# ============================================================
print(f"\n{'='*70}")
print("n=7: SAMPLED χ AND DIMENSION PROFILES")
print("=" * 70)

import random
n = 7
samples = 500
chi_dist = Counter()
dim_profiles = Counter()

for _ in range(samples):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    
    dims = compute_omega_dims(A, n)
    chi = sum((-1)**p * d for p, d in enumerate(dims))
    chi_dist[chi] += 1
    dim_profiles[tuple(dims)] += 1

print(f"  χ distribution: {dict(sorted(chi_dist.items()))}")
print(f"  Most common dim profiles:")
for dims, cnt in sorted(dim_profiles.items(), key=lambda x: -x[1])[:15]:
    chi = sum((-1)**p * d for p, d in enumerate(dims))
    print(f"    {dims} (χ={chi}): {cnt}")


# ============================================================
# CRITICAL: Rank analysis of ∂₂ and ∂₃ in Ω
# ============================================================
print(f"\n{'='*70}")
print("RANK ANALYSIS: rk(∂₂) + rk(∂₃) vs dim(Ω₂)")
print("=" * 70)

n = 5
m = n*(n-1)//2
total = 1 << m

defect_dist = Counter()  # dim(Ω₂) - rk(∂₂) - rk(∂₃) = β₂

for bits in range(total):
    A = build_adj(n, bits)
    
    ap = {}
    om = {}
    for p in range(5):
        ap[p] = enumerate_allowed_paths(A, n, p)
    om[0] = np.eye(n)
    for p in range(1, 5):
        if ap[p]:
            om[p] = compute_omega_basis(A, n, p, ap[p], ap[p-1])
        else:
            om[p] = np.zeros((0,0))
    
    d2 = dim_om(om[2])
    if d2 == 0:
        defect_dist[0] += 1
        continue
    
    # rk(∂₂: Ω₂ → Ω₁)
    bd2 = build_full_boundary_matrix(ap[2], ap[1])
    d2_mat = np.linalg.lstsq(om[1], bd2 @ om[2], rcond=None)[0]
    rk2 = np.linalg.matrix_rank(d2_mat, tol=1e-8)
    
    # rk(∂₃: Ω₃ → Ω₂)
    d3 = dim_om(om[3])
    if d3 > 0:
        bd3 = build_full_boundary_matrix(ap[3], ap[2])
        d3_mat = np.linalg.lstsq(om[2], bd3 @ om[3], rcond=None)[0]
        rk3 = np.linalg.matrix_rank(d3_mat, tol=1e-8)
    else:
        rk3 = 0
    
    defect = d2 - rk2 - rk3  # = β₂
    defect_dist[defect] += 1

print(f"\nn=5: dim(Ω₂) - rk(∂₂) - rk(∂₃) distribution: {dict(sorted(defect_dist.items()))}")
print(f"  β₂ = 0 for all: {'✓' if all(d == 0 for d in defect_dist) else '✗'}")


# ============================================================
# FORMULA HUNT: What determines dim(Ω₂), dim(Ω₃)?
# ============================================================
print(f"\n{'='*70}")
print("FORMULA: dim(Ω_p) vs tournament invariants")
print("=" * 70)

n = 5

# Count transitive triples and 3-cycles
from collections import defaultdict
t3_to_dims = defaultdict(list)

for bits in range(total):
    A = build_adj(n, bits)
    
    # Count 3-cycles
    t3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (A[i][j] and A[j][k] and A[k][i]) or (A[i][k] and A[k][j] and A[j][i]):
                    t3 += 1
    
    # Count transitive triples
    tt = 0
    for a in range(n):
        for b in range(n):
            if b == a or not A[a][b]: continue
            for c in range(n):
                if c in (a,b): continue
                if A[b][c] and A[a][c]:
                    tt += 1
    
    dims = compute_omega_dims(A, n, max_p=4)
    t3_to_dims[t3].append((tt, dims))

print(f"\nn=5: dim(Ω_p) by t₃:")
for t3 in sorted(t3_to_dims.keys()):
    examples = t3_to_dims[t3]
    tt_vals = set(e[0] for e in examples)
    dim_vals = Counter(tuple(e[1]) for e in examples)
    print(f"  t₃={t3}: {len(examples)} tournaments, tt={tt_vals}")
    for dims, cnt in sorted(dim_vals.items()):
        print(f"    Ω dims = {dims}: {cnt}")

# ============================================================
# KEY: Relation between dim(Ω₂) and transitive triples
# ============================================================
print(f"\n{'='*70}")
print("dim(Ω₂) = number of transitive triples?")
print("=" * 70)

all_match = True
for bits in range(total):
    A = build_adj(n, bits)
    tt = 0
    for a in range(n):
        for b in range(n):
            if b == a or not A[a][b]: continue
            for c in range(n):
                if c in (a,b): continue
                if A[b][c] and A[a][c]:
                    tt += 1
    
    dims = compute_omega_dims(A, n, max_p=2)
    if dims[2] != tt:
        print(f"  MISMATCH: bits={bits}, dim(Ω₂)={dims[2]}, tt={tt}")
        all_match = False

print(f"  dim(Ω₂) = #transitive_triples: {'✓' if all_match else '✗'}")

# Ω₃ = doubly transitive 4-paths?
print(f"\n  dim(Ω₃) = #doubly_transitive_4paths?")
all_match_3 = True
for bits in range(total):
    A = build_adj(n, bits)
    dt4 = 0
    for a in range(n):
        for b in range(n):
            if b == a or not A[a][b]: continue
            for c in range(n):
                if c in (a,b) or not A[b][c] or not A[a][c]: continue
                for d in range(n):
                    if d in (a,b,c) or not A[c][d] or not A[b][d]: continue
                    dt4 += 1
    
    dims = compute_omega_dims(A, n, max_p=3)
    if dims[3] != dt4:
        print(f"  MISMATCH: bits={bits}, dim(Ω₃)={dims[3]}, dt4={dt4}")
        all_match_3 = False

print(f"  dim(Ω₃) = #(a→b→c→d, a→c, b→d): {'✓' if all_match_3 else '✗'}")

# Formula for #transitive triples
print(f"\n  #transitive triples = C(n,3) + 2*C(n,3) - 3*t₃ ?")
print(f"  Actually: #tt = 3*C(n,3) - 3*t₃ = n(n-1)(n-2)/2 - 3*t₃")
print(f"  (Each 3-set contributes either 3 transitive triples (if transitive) or 0 (if 3-cycle))")
print(f"  Wait: a 3-set with a 3-cycle has 0 transitive triples. A transitive 3-set has 3.")
print(f"  So #tt = 3*(C(n,3) - t₃)")

all_match_formula = True
for bits in range(total):
    A = build_adj(n, bits)
    t3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (A[i][j] and A[j][k] and A[k][i]) or (A[i][k] and A[k][j] and A[j][i]):
                    t3 += 1
    tt_formula = 3 * (n*(n-1)*(n-2)//6 - t3)
    
    tt_actual = 0
    for a in range(n):
        for b in range(n):
            if b == a or not A[a][b]: continue
            for c in range(n):
                if c in (a,b): continue
                if A[b][c] and A[a][c]:
                    tt_actual += 1
    
    if tt_formula != tt_actual:
        print(f"  WRONG: bits={bits}, formula={tt_formula}, actual={tt_actual}")
        all_match_formula = False

print(f"  #tt = 3*(C(n,3) - t₃): {'✓' if all_match_formula else '✗'}")
print(f"  So dim(Ω₂) = 3*(C(n,3) - t₃)")

print("\nDone.")
