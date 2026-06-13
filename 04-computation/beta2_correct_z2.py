#!/usr/bin/env python3
"""
beta2_correct_z2.py - Correct Z2 analysis given beta1 != c3

Key finding: beta1 is NOT c3 for GLMY path homology of tournaments!
At n=5: beta1 in {0, 1}, NOT c3 in {0,...,5}.

So Z2 = dim(Om2) - rk(bd2|_Om2) must be computed directly.

The question for beta2=0 is: rk(bd3|_Om3) = Z2 = dim(Om2) - rk(bd2|_Om2).

We know beta2=0 holds computationally. What determines rk(bd2|_Om2)?

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


def compute_J2(A, n):
    J2 = 0
    for a in range(n):
        for c in range(n):
            if a == c or not A[c][a]:
                continue
            for b in range(n):
                if b != a and b != c and A[a][b] and A[b][c]:
                    J2 += 1
                    break
    return J2


# ============================================================
# Full Betti analysis at n=3,4,5,6
# ============================================================
for n in [3, 4, 5]:
    print("=" * 70)
    print(f"FULL BETTI NUMBERS: n={n}")
    print("=" * 70)

    total = 1 << (n*(n-1)//2)
    betti_dist = Counter()

    for bits in range(total):
        A = build_adj(n, bits)
        betti = path_betti_numbers(A, n, max_dim=n-1)
        c3 = count_c3(A, n)
        scores = tuple(sorted([sum(row) for row in A]))
        betti_dist[(scores, c3, tuple(betti))] += 1

    print(f"\n{'scores':<20} {'c3':>3} {'betti':<25} {'count':>5}")
    for (sc, c3, bt), cnt in sorted(betti_dist.items()):
        print(f"{str(sc):<20} {c3:>3} {str(list(bt)):<25} {cnt:>5}")


# ============================================================
# n=6: exhaustive betti numbers (just beta0, beta1, beta2)
# ============================================================
print(f"\n{'='*70}")
print("BETTI NUMBERS: n=6 (beta0, beta1, beta2)")
print("=" * 70)

n = 6
total = 1 << (n*(n-1)//2)

betti_dist6 = Counter()
t0 = time.time()
for bits in range(total):
    if bits % 10000 == 0 and bits > 0:
        dt = time.time() - t0
        print(f"  ... {bits}/{total} ({dt:.0f}s)")

    A = build_adj(n, bits)
    betti = path_betti_numbers(A, n, max_dim=2)
    c3 = count_c3(A, n)
    scores = tuple(sorted([sum(row) for row in A]))
    betti_dist6[(scores, c3, tuple(betti))] += 1

dt = time.time() - t0
print(f"\nDone in {dt:.0f}s")

print(f"\n{'scores':<25} {'c3':>3} {'betti':<15} {'count':>5}")
for (sc, c3, bt), cnt in sorted(betti_dist6.items()):
    print(f"{str(sc):<25} {c3:>3} {str(list(bt)):<15} {cnt:>5}")


# ============================================================
# ANALYSIS: What determines beta1?
# ============================================================
print(f"\n{'='*70}")
print("BETA1 ANALYSIS")
print("=" * 70)

# At n=5: beta1 in {0,1}
# Let's understand when beta1=1 vs beta1=0
n = 5
total = 1 << (n*(n-1)//2)

print(f"\nn={n}:")
b1_examples = defaultdict(list)
for bits in range(total):
    A = build_adj(n, bits)
    betti = path_betti_numbers(A, n, max_dim=2)
    c3 = count_c3(A, n)
    scores = tuple(sorted([sum(row) for row in A]))
    b1_examples[(scores, c3, betti[1])].append(bits)

print(f"{'scores':<20} {'c3':>3} {'b1':>3} {'count':>5}")
for (sc, c3, b1), examples in sorted(b1_examples.items()):
    print(f"{str(sc):<20} {c3:>3} {b1:>3} {len(examples):>5}")

# beta1=1 only occurs at specific score sequences
# Let's check: is there a pattern?
print(f"\nbeta1=1 scores at n=5: (1,1,2,3,3)[half], (1,2,2,2,3)[2/3], (2,2,2,2,2)[all]")
print("These all have c3 >= 3.")
print("But c3 >= 3 alone is not sufficient (scores (0,2,2,3,3) has c3=2, b1=0)")

# At n=4:
n = 4
total = 1 << (n*(n-1)//2)
print(f"\nn={n}:")
for bits in range(total):
    A = build_adj(n, bits)
    betti = path_betti_numbers(A, n, max_dim=3)
    c3 = count_c3(A, n)
    scores = tuple(sorted([sum(row) for row in A]))
    if bits < 64:
        pass

b1_dist4 = Counter()
for bits in range(total):
    A = build_adj(n, bits)
    betti = path_betti_numbers(A, n, max_dim=3)
    c3 = count_c3(A, n)
    scores = tuple(sorted([sum(row) for row in A]))
    b1_dist4[(scores, c3, tuple(betti))] += 1

print(f"{'scores':<15} {'c3':>3} {'betti':<25} {'count':>5}")
for (sc, c3, bt), cnt in sorted(b1_dist4.items()):
    print(f"{str(sc):<15} {c3:>3} {str(list(bt)):<25} {cnt:>5}")


# ============================================================
# Check if beta1 = c3 - (n-2) for some interpretation
# ============================================================
print(f"\n{'='*70}")
print("BETA1 FORMULA SEARCH")
print("=" * 70)

n = 5
total = 1 << (n*(n-1)//2)

print(f"\nn={n}:")
print(f"{'scores':<20} {'c3':>3} {'b1':>3} {'rk_bd2':>6} {'Om2':>4} {'Z2':>3}")
seen = set()
for bits in range(total):
    A = build_adj(n, bits)
    betti = path_betti_numbers(A, n, max_dim=2)
    c3 = count_c3(A, n)
    j2 = compute_J2(A, n)
    scores = tuple(sorted([sum(row) for row in A]))

    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    om2 = compute_omega_basis(A, n, 2, a2, a1)
    d_om2 = om2.shape[1] if om2.ndim == 2 else 0

    bd2 = build_full_boundary_matrix(a2, a1)
    if d_om2 > 0:
        S = np.linalg.svd(bd2 @ om2, compute_uv=False)
        rk_bd2 = sum(s > 1e-8 for s in S)
    else:
        rk_bd2 = 0

    dim_z2 = d_om2 - rk_bd2
    key = (scores, c3, betti[1], rk_bd2, d_om2, dim_z2)
    if key not in seen:
        seen.add(key)
        print(f"{str(scores):<20} {c3:>3} {betti[1]:>3} {rk_bd2:>6} {d_om2:>4} {dim_z2:>3}")


print("\n\nDone.")
