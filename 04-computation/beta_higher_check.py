#!/usr/bin/env python3
"""
beta_higher_check.py - Check if ALL beta_p = 0 for p >= 2 in tournaments

We know beta2=0 for all tested n. Check beta3, beta4 too.
Also develop a formula for beta1.

Author: kind-pasteur-2026-03-08-S41
"""
import sys, os, time
import numpy as np
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix, path_betti_numbers
)
sys.stdout = _saved


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


def count_c3(A, n):
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[k][j] and A[i][k]):
                    c3 += 1
    return c3


# ============================================================
# n=5: FULL Betti numbers (all dimensions)
# ============================================================
print("=" * 70)
print("FULL BETTI NUMBERS: n=5 (exhaustive)")
print("=" * 70)

n = 5
total = 1 << (n*(n-1)//2)

betti_all = Counter()
any_nonzero = False

for bits in range(total):
    A = build_adj(n, bits)
    betti = path_betti_numbers(A, n, max_dim=n-1)
    bt = tuple(betti)
    betti_all[bt] += 1
    if any(b > 0 for b in betti[2:]):
        any_nonzero = True
        print(f"  NONZERO higher: bits={bits}, betti={betti}")

print(f"\nAll Betti vectors:")
for bt, cnt in sorted(betti_all.items()):
    print(f"  {list(bt)}: {cnt}")
print(f"\nAny beta_p > 0 for p>=2: {any_nonzero}")


# ============================================================
# n=6: Check beta3, beta4, beta5 (sampled, since exhaustive is expensive)
# ============================================================
print(f"\n{'='*70}")
print("HIGHER BETTI NUMBERS: n=6 (sampled)")
print("=" * 70)

import random
random.seed(42)
n = 6
n_samples = 2000

higher_nonzero = 0
betti_dist6 = Counter()
t0 = time.time()

for trial in range(n_samples):
    bits = random.randint(0, (1 << (n*(n-1)//2)) - 1)
    A = build_adj(n, bits)
    betti = path_betti_numbers(A, n, max_dim=n-1)
    bt = tuple(betti)
    betti_dist6[bt] += 1
    if any(b > 0 for b in betti[2:]):
        higher_nonzero += 1
        if higher_nonzero <= 3:
            print(f"  NONZERO: bits={bits}, betti={betti}")

    if trial % 500 == 499:
        dt = time.time() - t0
        print(f"  ... {trial+1}/{n_samples} ({dt:.0f}s)")

dt = time.time() - t0
print(f"\nn=6: {n_samples} samples in {dt:.0f}s")
print(f"Any beta_p > 0 for p>=2: {higher_nonzero}/{n_samples}")

print(f"\nBetti vector distribution:")
for bt, cnt in sorted(betti_dist6.items(), key=lambda x: -x[1]):
    print(f"  {list(bt)}: {cnt}")


# ============================================================
# n=7: Just check beta2 and beta3 (sampled)
# ============================================================
print(f"\n{'='*70}")
print("BETTI NUMBERS: n=7 (sampled, max_dim=3)")
print("=" * 70)

n = 7
n_samples = 1000

higher_nonzero7 = 0
t0 = time.time()

for trial in range(n_samples):
    bits = random.randint(0, (1 << (n*(n-1)//2)) - 1)
    A = build_adj(n, bits)
    betti = path_betti_numbers(A, n, max_dim=3)
    if any(b > 0 for b in betti[2:]):
        higher_nonzero7 += 1
        if higher_nonzero7 <= 3:
            print(f"  NONZERO: bits={bits}, betti={betti}")

    if trial % 500 == 499:
        dt = time.time() - t0
        print(f"  ... {trial+1}/{n_samples} ({dt:.0f}s)")

dt = time.time() - t0
print(f"\nn=7: {n_samples} samples in {dt:.0f}s")
print(f"Any beta_p > 0 for p>=2,3: {higher_nonzero7}/{n_samples}")


# ============================================================
# BETA1 FORMULA: deeper analysis
# ============================================================
print(f"\n{'='*70}")
print("BETA1 FORMULA ANALYSIS")
print("=" * 70)

# Key pattern from n=4,5:
# beta1 relates to the 3-cycle STRUCTURE, not just count
# At n=4: beta1=1 iff (1,1,2,2) iff the tournament is a 4-cycle plus diagonal
# At n=5: beta1=0 or 1

# Hypothesis: beta1 = max(0, c3 - something)
# At n=4: c3=2 -> b1=1, so subtract 1? But c3=1 -> b1=0, so subtract >= 1
# At n=5: c3=5 -> b1=1, c3=4 -> b1=0 or 1, c3=3 -> b1=0 or 1

# Alternative: beta1 = c3 - rk(bd2|_Om2) - dim(something)
# We know: beta1 = dim(Z1) - rk(bd2) = C(n-1,2) - rk(bd2)
# So beta1 is determined by rk(bd2|_Om2).

# What determines rk(bd2|_Om2)?
# bd2 maps Om2 -> A1, and rk = dim(Om2) - dim(Z2)
# So rk(bd2) = dim(Om2) - dim(Z2)

# From our data at n=5:
# rk(bd2) in {5, 6}
# rk(bd2) = 6 = C(4,2): means bd2|_Om2 is full rank (image = Z1)
# rk(bd2) = 5: beta1 = 1

# When is bd2|_Om2 full rank (= C(n-1,2))?
# This means im(bd2|_Om2) = Z1, i.e., EVERY Z1 element is a boundary

# At n=5, rk=6 iff beta1=0, rk=5 iff beta1=1
# beta1 counts independent 1-cycles NOT filled by 2-boundaries

# Let's check: is beta1 related to the clique cover number?
# Or to the directed cycle structure?

n = 5
total = 1 << (n*(n-1)//2)

# For each tournament, compute directed cycle structure
for bits in [0, 2, 31, 100, 341, 682, 1023]:
    A = build_adj(n, bits)
    betti = path_betti_numbers(A, n, max_dim=4)
    c3 = count_c3(A, n)
    scores = tuple(sorted([sum(row) for row in A]))

    # Count directed 5-cycles
    a4 = enumerate_allowed_paths(A, n, 4)
    c5 = 0
    for p in a4:
        if A[p[4]][p[0]]:  # closes the cycle
            c5 += 1
    c5 //= 5  # each cycle counted 5 times

    print(f"  bits={bits}, scores={scores}, c3={c3}, c5={c5}, betti={betti}")

# At n=6, what are the possible beta1 values?
print(f"\nbeta1 distribution at n=6:")
n = 6
b1_dist = Counter()
for bits in range(1 << (n*(n-1)//2)):
    A = build_adj(n, bits)
    a0 = enumerate_allowed_paths(A, n, 0)
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)

    om1 = compute_omega_basis(A, n, 1, a1, a0)
    om2 = compute_omega_basis(A, n, 2, a2, a1)
    d_om1 = om1.shape[1] if om1.ndim == 2 else 0
    d_om2 = om2.shape[1] if om2.ndim == 2 else 0

    bd1 = build_full_boundary_matrix(a1, a0)
    bd2 = build_full_boundary_matrix(a2, a1)

    rk1 = 0
    if d_om1 > 0:
        S = np.linalg.svd(bd1 @ om1, compute_uv=False)
        rk1 = sum(s > 1e-8 for s in S)

    rk2 = 0
    if d_om2 > 0:
        S = np.linalg.svd(bd2 @ om2, compute_uv=False)
        rk2 = sum(s > 1e-8 for s in S)

    b1 = (d_om1 - rk1) - rk2
    c3 = count_c3(A, n)
    scores = tuple(sorted([sum(row) for row in A]))
    b1_dist[(scores, c3, b1)] += 1

print(f"\n{'scores':<25} {'c3':>3} {'b1':>3} {'count':>5}")
for (sc, c3, b1), cnt in sorted(b1_dist.items()):
    if b1 > 0:
        print(f"{str(sc):<25} {c3:>3} {b1:>3} {cnt:>5} ***")

# Count total beta1 > 0
b1_pos = sum(cnt for (sc, c3, b1), cnt in b1_dist.items() if b1 > 0)
print(f"\nbeta1 > 0: {b1_pos}/{1 << (n*(n-1)//2)}")

# Max beta1
max_b1 = max(b1 for (sc, c3, b1) in b1_dist.keys())
print(f"Max beta1 at n=6: {max_b1}")


print("\n\nDone.")
