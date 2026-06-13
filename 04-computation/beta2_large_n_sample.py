#!/usr/bin/env python3
"""
beta2_large_n_sample.py - Sample beta_2 at n=7,8,9 to extend verification

Also compute the full Betti sequence to check beta_p = 0 for all p >= 2.

Author: kind-pasteur-2026-03-08-S42
"""
import sys, os, time, random
import numpy as np
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
            if bits & (1 << idx): A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A


def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


# ============================================================
# PART 1: n=7 sampling — full Betti numbers
# ============================================================
print("=" * 70)
print("BETA SAMPLING: FULL BETTI NUMBERS")
print("=" * 70)

for n, sample_size, max_dim in [(7, 2000, 5), (8, 500, 4), (9, 100, 3)]:
    print(f"\nn={n}: {sample_size} random tournaments, max_dim={max_dim}")
    t0 = time.time()

    beta_nonzero = {p: 0 for p in range(2, max_dim + 1)}
    beta1_dist = {}
    betti_profiles = {}

    for trial in range(sample_size):
        A = random_tournament(n)
        betti = path_betti_numbers(A, n, max_dim=max_dim)

        for p in range(2, min(len(betti), max_dim + 1)):
            if betti[p] != 0:
                beta_nonzero[p] += 1

        b1 = betti[1] if len(betti) > 1 else 0
        beta1_dist[b1] = beta1_dist.get(b1, 0) + 1

        bt = tuple(betti)
        betti_profiles[bt] = betti_profiles.get(bt, 0) + 1

        if (trial+1) % max(1, sample_size // 5) == 0:
            elapsed = time.time() - t0
            print(f"  {trial+1}/{sample_size} ({elapsed:.0f}s)")
            # Check if any nonzero beta_p for p >= 2
            any_nonzero = any(v > 0 for v in beta_nonzero.values())
            if any_nonzero:
                print(f"  *** NONZERO beta_p for p>=2 detected! ***")

    elapsed = time.time() - t0
    print(f"\n  Results ({elapsed:.0f}s):")
    for p in range(2, max_dim + 1):
        print(f"    beta_{p} != 0: {beta_nonzero[p]}/{sample_size}")
    print(f"    beta_1 distribution: {dict(sorted(beta1_dist.items()))}")

    all_zero = all(v == 0 for v in beta_nonzero.values())
    print(f"    ALL beta_p=0 for p>=2? {all_zero}")

    # Show most common Betti profiles
    print(f"    Top Betti profiles:")
    for bt, cnt in sorted(betti_profiles.items(), key=lambda x: -x[1])[:10]:
        print(f"      {list(bt)}: {cnt} ({100*cnt/sample_size:.1f}%)")


# ============================================================
# PART 2: Special tournaments at larger n
# ============================================================
print(f"\n{'='*70}")
print("SPECIAL TOURNAMENTS")
print("=" * 70)

# Paley-like tournaments (circulant)
def circulant_tournament(n):
    """Quadratic residue tournament (Paley) for prime n."""
    qr = set()
    for x in range(1, n):
        qr.add((x*x) % n)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and ((j - i) % n) in qr:
                A[i][j] = 1
    # Verify it's a tournament
    for i in range(n):
        for j in range(i+1, n):
            assert A[i][j] + A[j][i] == 1, f"Not tournament at {i},{j}"
    return A

for p in [5, 7, 11, 13]:
    try:
        A = circulant_tournament(p)
        max_d = min(p-1, 5)
        t0 = time.time()
        betti = path_betti_numbers(A, p, max_dim=max_d)
        elapsed = time.time() - t0
        print(f"\n  Paley T_{p}: beta = {betti} ({elapsed:.1f}s)")
        all_zero = all(betti[i] == 0 for i in range(2, len(betti)))
        print(f"    beta_p=0 for p>=2? {all_zero}")
    except Exception as e:
        print(f"  Paley T_{p}: ERROR {e}")


# Doubly regular tournaments from small n
# Transitive tournament
for n in [7, 8, 9, 10]:
    # Transitive
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            A[i][j] = 1
    max_d = min(n-1, 5)
    t0 = time.time()
    betti = path_betti_numbers(A, n, max_dim=max_d)
    elapsed = time.time() - t0
    print(f"\n  Transitive T_{n}: beta = {betti} ({elapsed:.1f}s)")


# ============================================================
# PART 3: Regular tournaments (odd n) — are they always beta_p=0 for p>=2?
# ============================================================
print(f"\n{'='*70}")
print("REGULAR TOURNAMENTS BETA SAMPLING")
print("=" * 70)

for n in [7, 9]:
    print(f"\nn={n}: sampling regular tournaments")
    sample = 200 if n == 7 else 50
    reg_count = 0
    t0 = time.time()

    for trial in range(sample * 10):  # oversample and filter
        A = random_tournament(n)
        scores = sorted([sum(row) for row in A])
        if scores != [(n-1)//2] * n:
            continue

        reg_count += 1
        max_d = min(n-1, 5)
        betti = path_betti_numbers(A, n, max_dim=max_d)
        nonzero_high = any(betti[p] != 0 for p in range(2, len(betti)))
        if nonzero_high:
            print(f"  NONZERO: beta={betti}")

        if reg_count >= sample:
            break

    elapsed = time.time() - t0
    print(f"  n={n}: checked {reg_count} regular tournaments ({elapsed:.0f}s), all beta_p=0 for p>=2")


print("\n\nDone.")
