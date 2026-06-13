#!/usr/bin/env python3
"""
beta2_n7_verify.py - Verify beta2=0 at n=7 (large sample)

Also: test if the "contractibility" / "quasi-contractibility" of
tournaments forces exactness at dimension 2.

Key idea: A contractible chain complex has all betti = 0.
A "quasi-contractible at dim 2" complex has beta2 = 0.

Author: kind-pasteur-2026-03-08-S41
"""
import sys, os, time
import numpy as np
import random
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import path_betti_numbers
sys.stdout = _saved

random.seed(42)


def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A


# ============================================================
# Verify beta2=0 at n=7 (large sample)
# ============================================================
print("=" * 70)
print("BETA2=0 VERIFICATION AT n=7")
print("=" * 70)

n = 7
total = 1 << (n*(n-1)//2)  # 2^21 = 2097152
sample_size = 2000

beta2_nonzero = 0
betti_dist = Counter()
t0 = time.time()

for trial in range(sample_size):
    bits = random.randint(0, total - 1)
    A = build_adj(n, bits)
    betti = path_betti_numbers(A, n, max_dim=2)
    b2 = betti[2]

    if b2 > 0:
        beta2_nonzero += 1
        scores = tuple(sorted([sum(row) for row in A]))
        print(f"  COUNTEREXAMPLE at trial {trial}: bits={bits}, scores={scores}, betti={betti}")

    betti_dist[tuple(betti[:3])] += 1

    if (trial + 1) % 500 == 0:
        elapsed = time.time() - t0
        print(f"  {trial+1}/{sample_size} done ({elapsed:.1f}s), beta2>0: {beta2_nonzero}")

elapsed = time.time() - t0
print(f"\nn=7: beta2>0 in {beta2_nonzero}/{sample_size} random tournaments ({elapsed:.1f}s)")
print(f"Betti (b0,b1,b2) distribution: {dict(sorted(betti_dist.items()))}")


# ============================================================
# Also check n=8 (smaller sample)
# ============================================================
print(f"\n{'='*70}")
print("BETA2=0 VERIFICATION AT n=8")
print("=" * 70)

n = 8
total = 1 << (n*(n-1)//2)  # 2^28
sample_size = 500

beta2_nonzero = 0
betti_dist = Counter()
t0 = time.time()

for trial in range(sample_size):
    bits = random.randint(0, total - 1)
    A = build_adj(n, bits)
    betti = path_betti_numbers(A, n, max_dim=2)
    b2 = betti[2]

    if b2 > 0:
        beta2_nonzero += 1
        scores = tuple(sorted([sum(row) for row in A]))
        print(f"  COUNTEREXAMPLE: bits={bits}, scores={scores}, betti={betti}")

    betti_dist[tuple(betti[:3])] += 1

    if (trial + 1) % 100 == 0:
        elapsed = time.time() - t0
        print(f"  {trial+1}/{sample_size} done ({elapsed:.1f}s), beta2>0: {beta2_nonzero}")

elapsed = time.time() - t0
print(f"\nn=8: beta2>0 in {beta2_nonzero}/{sample_size} random tournaments ({elapsed:.1f}s)")
print(f"Betti (b0,b1,b2) distribution: {dict(sorted(betti_dist.items()))}")


# ============================================================
# Check n=9 (very small sample)
# ============================================================
print(f"\n{'='*70}")
print("BETA2=0 VERIFICATION AT n=9")
print("=" * 70)

n = 9
total = 1 << (n*(n-1)//2)  # 2^36
sample_size = 100

beta2_nonzero = 0
t0 = time.time()

for trial in range(sample_size):
    bits = random.randint(0, total - 1)
    A = build_adj(n, bits)
    betti = path_betti_numbers(A, n, max_dim=2)
    b2 = betti[2]

    if b2 > 0:
        beta2_nonzero += 1
        scores = tuple(sorted([sum(row) for row in A]))
        print(f"  COUNTEREXAMPLE: bits={bits}, scores={scores}")

elapsed = time.time() - t0
print(f"\nn=9: beta2>0 in {beta2_nonzero}/{sample_size} random tournaments ({elapsed:.1f}s)")


print("\n\nDone.")
