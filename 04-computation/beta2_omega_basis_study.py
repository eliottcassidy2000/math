#!/usr/bin/env python3
"""
beta2_omega_basis_study.py — What ARE the Ω₂ and Ω₃ basis elements?

The previous script showed dim(Ω₂) ≠ #transitive_triples and 
dim(Ω₃) ≠ #doubly_transitive_4paths. So Ω contains COMBINATIONS
of allowed paths that aren't individually in Ω.

Study the actual basis vectors to understand the structure.

Author: opus-2026-03-08-S49
"""
import sys, numpy as np
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

# Pick a tournament where dim(Ω₂) ≠ #transitive triples
# bits=4 at n=5 had dim(Ω₃)=4 but dt4=3
n = 5
bits = 4
A = build_adj(n, bits)
print(f"Tournament bits={bits}, n={n}")
print(f"Adjacency:")
for i in range(n):
    row = [A[i][j] for j in range(n)]
    print(f"  {i}: {row}")
scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]
print(f"Scores: {scores}")

# Compute Ω₂
ap1 = enumerate_allowed_paths(A, n, 1)
ap2 = enumerate_allowed_paths(A, n, 2)
ap3 = enumerate_allowed_paths(A, n, 3)

om1 = compute_omega_basis(A, n, 1, ap1, list(range(n)))
om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
om3 = compute_omega_basis(A, n, 3, ap3, ap2) if ap3 else np.zeros((0,0))

print(f"\n|A₁|={len(ap1)}, dim(Ω₁)={dim_om(om1)}")
print(f"|A₂|={len(ap2)}, dim(Ω₂)={dim_om(om2)}")
print(f"|A₃|={len(ap3)}, dim(Ω₃)={dim_om(om3)}")

# Show A₂ paths
print(f"\nA₂ (allowed 2-paths):")
tt_count = 0
for p in ap2:
    a,b,c = p
    is_trans = bool(A[a][c])
    tag = "TT" if is_trans else "IT"
    if is_trans: tt_count += 1
    print(f"  ({a},{b},{c}) [{tag}]")
print(f"  Transitive: {tt_count}, Intransitive: {len(ap2)-tt_count}")

# Show Ω₂ basis
print(f"\nΩ₂ basis ({dim_om(om2)} vectors in ℝ^{len(ap2)}):")
for j in range(dim_om(om2)):
    col = om2[:, j]
    nonzero = [(ap2[i], round(col[i], 4)) for i in range(len(col)) if abs(col[i]) > 1e-8]
    print(f"  ω₂[{j}]: {nonzero}")

# Show A₃ paths  
print(f"\nA₃ (allowed 3-paths):")
for p in ap3:
    a,b,c,d = p
    ac = A[a][c]
    bd = A[b][d]
    ad = A[a][d]
    tag = f"a→c={ac}, b→d={bd}, a→d={ad}"
    print(f"  ({a},{b},{c},{d}) [{tag}]")

# Show Ω₃ basis
print(f"\nΩ₃ basis ({dim_om(om3)} vectors in ℝ^{len(ap3)}):")
for j in range(dim_om(om3)):
    col = om3[:, j]
    nonzero = [(ap3[i], round(col[i], 4)) for i in range(len(col)) if abs(col[i]) > 1e-8]
    print(f"  ω₃[{j}]: {nonzero}")

# Verify: which individual A₂ paths are in Ω₂?
print(f"\nIndividual A₂ paths in Ω₂?")
for i, p in enumerate(ap2):
    e_i = np.zeros(len(ap2))
    e_i[i] = 1.0
    # Check if e_i is in column space of om2
    if dim_om(om2) > 0:
        coords = np.linalg.lstsq(om2, e_i, rcond=None)[0]
        resid = np.linalg.norm(om2 @ coords - e_i)
        in_omega = resid < 1e-6
    else:
        in_omega = False
    a,b,c = p
    tag = "TT" if A[a][c] else "IT"
    print(f"  ({a},{b},{c}) [{tag}]: {'✓ in Ω₂' if in_omega else '✗ NOT in Ω₂'}")

print("\nDone.")
