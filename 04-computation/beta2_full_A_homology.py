#!/usr/bin/env python3
"""
beta2_full_A_homology.py — Homology of the FULL allowed path complex A_*

The key question: is H₂(A_*) = 0 for all tournaments?

If YES, then Z₂(A) = B₂(A) = im(∂₃|A₃).
Combined with im(∂₃|A₃) ∩ Ω₂ = im(∂₃|Ω₃), this gives β₂(Ω) = 0.

The A_* complex:
- A_p = allowed p-paths in the tournament (directed, all edges present)
- ∂_p: alternating sum of face maps (delete vertex i from position i)
- H_p(A) = ker(∂_p)/im(∂_{p+1})

Note: A_* is NOT the Ω complex (which has the ∂∘∂=0 constraint).
The A complex is the "regular path complex" of the digraph.

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
    enumerate_allowed_paths,
    build_full_boundary_matrix
)
sys.stdout = _saved

def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A


print("=" * 70)
print("HOMOLOGY OF FULL ALLOWED PATH COMPLEX A_*")
print("=" * 70)

for n in [3, 4, 5]:
    m = n*(n-1)//2
    total = 1 << m
    t0 = time.time()

    h2_dist = Counter()
    h1_dist = Counter()
    h3_dist = Counter()

    for bits in range(total):
        A = build_adj(n, bits)

        # Compute all allowed paths
        ap = {}
        for p in range(n):
            ap[p] = enumerate_allowed_paths(A, n, p)

        # Build boundary matrices
        bds = {}
        for p in range(1, n):
            if ap[p] and ap[p-1]:
                bds[p] = build_full_boundary_matrix(ap[p], ap[p-1])
            else:
                bds[p] = np.zeros((len(ap[p-1]) if ap.get(p-1) else 0,
                                   len(ap[p]) if ap.get(p) else 0))

        # H₁(A): ker(∂₁)/im(∂₂)
        if len(ap[1]) > 0:
            rk_d1 = np.linalg.matrix_rank(bds[1], tol=1e-8) if bds[1].size > 0 else 0
            ker_d1 = len(ap[1]) - rk_d1
        else:
            ker_d1 = 0

        if len(ap[2]) > 0 and len(ap[1]) > 0:
            rk_d2 = np.linalg.matrix_rank(bds[2], tol=1e-8) if bds[2].size > 0 else 0
        else:
            rk_d2 = 0

        h1 = ker_d1 - rk_d2
        h1_dist[h1] += 1

        # H₂(A): ker(∂₂)/im(∂₃)
        ker_d2 = len(ap[2]) - rk_d2 if len(ap[2]) > 0 else 0

        if len(ap.get(3, [])) > 0 and len(ap[2]) > 0:
            rk_d3 = np.linalg.matrix_rank(bds[3], tol=1e-8) if bds[3].size > 0 else 0
        else:
            rk_d3 = 0

        h2 = ker_d2 - rk_d3
        h2_dist[h2] += 1

        # H₃(A): ker(∂₃)/im(∂₄) (if applicable)
        if n >= 5:
            ker_d3 = len(ap.get(3, [])) - rk_d3 if len(ap.get(3, [])) > 0 else 0
            if len(ap.get(4, [])) > 0 and len(ap.get(3, [])) > 0:
                rk_d4 = np.linalg.matrix_rank(bds[4], tol=1e-8) if bds[4].size > 0 else 0
            else:
                rk_d4 = 0
            h3 = ker_d3 - rk_d4
            h3_dist[h3] += 1

    elapsed = time.time() - t0
    print(f"\nn={n}: {total} tournaments in {elapsed:.1f}s")
    print(f"  H₁(A): {dict(sorted(h1_dist.items()))}")
    print(f"  H₂(A): {dict(sorted(h2_dist.items()))}")
    if h3_dist:
        print(f"  H₃(A): {dict(sorted(h3_dist.items()))}")

    if all(v == 0 for v in h2_dist.keys() if h2_dist[v] > 0):
        print(f"  *** H₂(A_*) = 0 for ALL n={n} tournaments! ***")

# n=6 exhaustive
print(f"\n{'='*70}")
print("n=6 EXHAUSTIVE")
print("=" * 70)

n = 6
m = n*(n-1)//2
total = 1 << m
t0 = time.time()

h2_dist = Counter()
h2_nonzero_examples = []

for bits in range(total):
    A = build_adj(n, bits)

    ap = {}
    for p in range(n):
        ap[p] = enumerate_allowed_paths(A, n, p)

    bds = {}
    for p in [2, 3]:
        if ap.get(p) and ap.get(p-1):
            bds[p] = build_full_boundary_matrix(ap[p], ap[p-1])

    # H₂(A)
    rk_d2 = np.linalg.matrix_rank(bds[2], tol=1e-8) if 2 in bds and bds[2].size > 0 else 0
    ker_d2 = len(ap[2]) - rk_d2

    rk_d3 = np.linalg.matrix_rank(bds[3], tol=1e-8) if 3 in bds and bds[3].size > 0 else 0
    h2 = ker_d2 - rk_d3

    h2_dist[h2] += 1
    if h2 > 0 and len(h2_nonzero_examples) < 3:
        scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]
        h2_nonzero_examples.append((bits, scores, h2))

    if (bits+1) % 10000 == 0:
        elapsed = time.time() - t0
        print(f"  {bits+1}/{total} ({elapsed:.0f}s)")

elapsed = time.time() - t0
print(f"\nn={n}: {total} tournaments in {elapsed:.0f}s")
print(f"  H₂(A): {dict(sorted(h2_dist.items()))}")

if h2_nonzero_examples:
    print(f"  Examples with H₂(A)>0:")
    for b, s, h in h2_nonzero_examples:
        print(f"    bits={b}, scores={s}, H₂={h}")
else:
    print(f"  *** H₂(A_*) = 0 for ALL n={n} tournaments! ***")

print("\nDone.")
