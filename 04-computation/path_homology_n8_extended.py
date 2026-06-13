#!/usr/bin/env python3
"""Extended n=8 tournament sampling - looking for pattern in which β_p appear"""
import numpy as np
from collections import Counter
import random
import sys
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import path_betti_numbers, count_3cycles

random.seed(123)  # different seed for fresh samples

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

n = 8
betti_dist = Counter()
b_by_dim = {p: Counter() for p in range(8)}
co_occur = Counter()  # which (p,q) pairs both nonzero
total = 0

for trial in range(50):
    A = random_tournament(n)
    betti = path_betti_numbers(A, n, max_dim=7)
    bt = tuple(betti)
    betti_dist[bt] += 1
    total += 1

    nonzero_dims = [p for p in range(len(betti)) if betti[p] > 0 and p > 0]
    for p in nonzero_dims:
        b_by_dim[p][betti[p]] += 1
    for p in nonzero_dims:
        for q in nonzero_dims:
            if p < q:
                co_occur[(p,q)] += 1

    if trial % 10 == 9:
        print(f"  ... {trial+1} done", flush=True)

print(f"\nn=8 Betti distribution ({total} samples):")
for bt in sorted(betti_dist.keys()):
    print(f"  β={list(bt)}: {betti_dist[bt]}")

print(f"\nNonzero β_p by dimension:")
for p in range(1, 8):
    if b_by_dim[p]:
        print(f"  β_{p}: {dict(b_by_dim[p])}")
    else:
        print(f"  β_{p}: never nonzero")

print(f"\nCo-occurrences (β_p>0 AND β_q>0 simultaneously):")
for (p,q), cnt in sorted(co_occur.items()):
    print(f"  β_{p} and β_{q}: {cnt}")

# Also check: mutual exclusion patterns
print(f"\nMutual exclusion check:")
for bt in betti_dist:
    b = list(bt)
    nz = [p for p in range(1, len(b)) if b[p] > 0]
    if len(nz) > 1:
        print(f"  MULTIPLE nonzero: β={b}, nonzero at p={nz}")

print("\nDone.")
