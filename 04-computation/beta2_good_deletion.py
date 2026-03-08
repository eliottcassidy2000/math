#!/usr/bin/env python3
"""
beta2_good_deletion.py - THE KEY QUESTION for beta2=0 proof:

Does every tournament T on n vertices have a vertex v such that
beta1(T\v) = 0?

If YES, then beta2=0 follows by induction:
  - Base: beta2=0 at n<=3 (trivially)
  - Inductive step: Pick v with beta1(T\v) = 0.
    By induction: beta2(T\v) = 0.
    H2(T,T\v) = 0 (since H2_rel > 0 requires beta1(T\v) = 1).
    LES: 0 = H2(T\v) -> H2(T) -> H2(T,T\v) = 0
    => H2(T) = 0. QED!

This is the SIMPLEST possible proof structure for beta2=0.

Author: kind-pasteur-2026-03-08-S42
"""
import sys, os, time, random
import numpy as np
from collections import Counter
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


# ============================================================
# CHECK: Does every tournament have a "good" vertex v with
#        beta1(T\v) = 0?
# ============================================================
print("=" * 70)
print("GOOD DELETION TEST: exists v with beta1(T\\v) = 0?")
print("=" * 70)

for n in [4, 5, 6]:
    m = n*(n-1)//2
    total = 1 << m

    no_good_v = 0
    good_count_dist = Counter()
    b1_all_dist = Counter()
    t0 = time.time()

    for bits in range(total):
        A = build_adj(n, bits)
        has_good = False
        good_count = 0

        for v in range(n):
            others = [i for i in range(n) if i != v]
            A_sub = [[A[others[i]][others[j]] for j in range(n-1)] for i in range(n-1)]
            b = path_betti_numbers(A_sub, n-1, max_dim=1)
            b1 = b[1]
            if b1 == 0:
                has_good = True
                good_count += 1

        if not has_good:
            no_good_v += 1
            scores = tuple(sorted([sum(row) for row in A]))
            if no_good_v <= 10:
                print(f"  NO GOOD V! T#{bits} scores={scores}")
                # Show all deletion betti numbers
                for v in range(n):
                    others = [i for i in range(n) if i != v]
                    A_sub = [[A[others[i]][others[j]] for j in range(n-1)] for i in range(n-1)]
                    b = path_betti_numbers(A_sub, n-1, max_dim=2)
                    dv = sum(A[v])
                    print(f"    v={v} (d+={dv}): betti={b}")

        good_count_dist[good_count] += 1

    elapsed = time.time() - t0
    print(f"\nn={n}: {total} tournaments ({elapsed:.0f}s)")
    print(f"  Tournaments with NO good vertex: {no_good_v}/{total}")
    print(f"  Good vertex count distribution: {dict(sorted(good_count_dist.items()))}")

    if no_good_v == 0:
        print(f"  ==> EVERY tournament at n={n} has a good vertex!")


# ============================================================
# Sample at n=7
# ============================================================
print(f"\n{'='*70}")
print("GOOD DELETION: SAMPLING AT n=7")
print("=" * 70)

n = 7
total = 1 << (n*(n-1)//2)
sample_size = 5000

no_good_v = 0
good_count_dist = Counter()
t0 = time.time()

for trial in range(sample_size):
    bits = random.randint(0, total - 1)
    A = build_adj(n, bits)
    good_count = 0
    has_good = False

    for v in range(n):
        others = [i for i in range(n) if i != v]
        A_sub = [[A[others[i]][others[j]] for j in range(6)] for i in range(6)]
        b = path_betti_numbers(A_sub, 6, max_dim=1)
        if b[1] == 0:
            has_good = True
            good_count += 1

    if not has_good:
        no_good_v += 1

    good_count_dist[good_count] += 1

    if (trial+1) % 1000 == 0:
        elapsed = time.time() - t0
        print(f"  {trial+1}/{sample_size} ({elapsed:.0f}s), no_good={no_good_v}")

elapsed = time.time() - t0
print(f"\nn=7: {sample_size} samples ({elapsed:.0f}s)")
print(f"  No good vertex: {no_good_v}/{sample_size}")
print(f"  Good vertex count dist: {dict(sorted(good_count_dist.items()))}")


# ============================================================
# Sample at n=8
# ============================================================
print(f"\n{'='*70}")
print("GOOD DELETION: SAMPLING AT n=8")
print("=" * 70)

n = 8
total = 1 << (n*(n-1)//2)
sample_size = 2000

no_good_v = 0
good_count_dist = Counter()
t0 = time.time()

for trial in range(sample_size):
    bits = random.randint(0, total - 1)
    A = build_adj(n, bits)
    good_count = 0
    has_good = False

    for v in range(n):
        others = [i for i in range(n) if i != v]
        A_sub = [[A[others[i]][others[j]] for j in range(7)] for i in range(7)]
        b = path_betti_numbers(A_sub, 7, max_dim=1)
        if b[1] == 0:
            has_good = True
            good_count += 1
            break  # only need one

    if not has_good:
        no_good_v += 1
        scores = tuple(sorted([sum(row) for row in A]))
        print(f"  NO GOOD V! bits={bits} scores={scores}")

    good_count_dist[good_count] += 1

    if (trial+1) % 500 == 0:
        elapsed = time.time() - t0
        print(f"  {trial+1}/{sample_size} ({elapsed:.0f}s), no_good={no_good_v}")

elapsed = time.time() - t0
print(f"\nn=8: {sample_size} samples ({elapsed:.0f}s)")
print(f"  No good vertex: {no_good_v}/{sample_size}")


# ============================================================
# STRUCTURAL QUESTION: Which tournaments have ALL deletions
#                      with beta1=1?
# ============================================================
print(f"\n{'='*70}")
print("WHICH TOURNAMENTS HAVE ALL beta1(T\\v) = 1?")
print("=" * 70)

# At n-1, beta1=1 happens when the tournament has a specific structure.
# beta1=0 iff the path complex is "contractible at dim 1".
# beta1=1 iff there's a nontrivial 1-cycle.
# At n=4: beta1=1 for score (1,1,2,2) only (24/64 = 37.5%)
# At n=5: beta1=1 for 304/1024 = 29.7%
# At n=6: beta1=1 for ?

# If beta1(T')=1 has probability p at n-1, and deletions were independent,
# then Prob(all bad) = p^n. At n=6, p~0.3, so p^6 ~ 0.0007.
# But deletions are NOT independent, so the actual count could be higher.

n = 5
count_all_b1 = 0
for bits in range(1 << (n*(n-1)//2)):
    A = build_adj(n, bits)
    all_b1 = True
    for v in range(n):
        others = [i for i in range(n) if i != v]
        A_sub = [[A[others[i]][others[j]] for j in range(4)] for i in range(4)]
        b = path_betti_numbers(A_sub, 4, max_dim=1)
        if b[1] == 0:
            all_b1 = False
            break
    if all_b1:
        count_all_b1 += 1
        scores = tuple(sorted([sum(row) for row in A]))
        print(f"  n={n}: T#{bits} scores={scores} has ALL beta1(T\\v)=1")

print(f"\nn={n}: {count_all_b1}/{1 << (n*(n-1)//2)} have all beta1(T\\v)=1")

# Check n=6
n = 6
count_all_b1 = 0
t0 = time.time()
for bits in range(1 << (n*(n-1)//2)):
    A = build_adj(n, bits)
    all_b1 = True
    for v in range(n):
        others = [i for i in range(n) if i != v]
        A_sub = [[A[others[i]][others[j]] for j in range(5)] for i in range(5)]
        b = path_betti_numbers(A_sub, 5, max_dim=1)
        if b[1] == 0:
            all_b1 = False
            break
    if all_b1:
        count_all_b1 += 1
        scores = tuple(sorted([sum(row) for row in A]))
        if count_all_b1 <= 5:
            print(f"  n={n}: T#{bits} scores={scores}")

    if (bits+1) % 10000 == 0:
        elapsed = time.time() - t0
        print(f"  n={n}: {bits+1}/{1<<(n*(n-1)//2)} ({elapsed:.0f}s), all_b1={count_all_b1}")

elapsed = time.time() - t0
print(f"\nn={n}: {count_all_b1}/{1 << (n*(n-1)//2)} have all beta1(T\\v)=1 ({elapsed:.0f}s)")


print("\n\nDone.")
