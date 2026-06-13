#!/usr/bin/env python3
"""
maximizer_hereditary_topology.py — Topological hereditary property

Discovery: β_{k+1}(T) > 0 when ALL vertex-deletions have β_k > 0.

Test this hypothesis at n=6 (exhaustive):
1. For each S-phase tournament (β₃=1) at n=6, do ALL vertex-deletions have β₂>0?
   (This would mean the topological hereditary goes further down)
2. For each C-phase tournament (β₁=1) at n=6, do ALL vertex-deletions have β₀=1?
   (Trivially yes, since β₀=1 always)

Also test: at n=5, which tournaments have β₁=1?
Do ALL their vertex-deletions have β₀=1? (Yes trivially)
Do THEIR vertex-deletions (n=4) have any pattern?

The hereditary chain should be:
  n=7 BIBD: β₄=6 ← n=6 deletions all β₃=1
  n=6 S-phase: β₃=1 ← n=5 deletions all β₂=?
  n=5 C-phase: β₁=1 ← n=4 deletions β₀=1 (trivial)

Wait, at n=5 the maximizers have β₁=1, and at n=6 the S-phase maximizers have β₃=1.
The dimension JUMPS from 1 to 3, skipping β₂. So the chain might be:
  β₄ ← β₃ ← β₁ ← β₀ (with a gap at dim 2)

Let's check if β₂ ever appears for tournaments.

Author: kind-pasteur-2026-03-08-S40
"""
import sys, os, time
from collections import Counter
from itertools import combinations
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w')
from path_homology_v2 import path_betti_numbers
sys.stdout = _saved

def H_tournament(A, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if dp[mask][v] == 0: continue
            for u in range(n):
                if mask & (1 << u): continue
                if A[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full])

# ===== Check: does beta_2 EVER appear for tournaments? =====
print("=" * 70)
print("BETA_2 EXISTENCE CHECK at n=5,6")
print("=" * 70)

# n=5 exhaustive
n = 5
m = n * (n-1) // 2
total = 1 << m
beta2_count = 0

for bits in range(total):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    try:
        beta = path_betti_numbers(A, n, max_dim=3)
        if len(beta) > 2 and int(beta[2]) > 0:
            beta2_count += 1
            H = H_tournament(A, n)
            beta_list = [int(beta[k]) for k in range(min(4, len(beta)))]
            print(f"  n=5, H={H}, beta={beta_list}")
    except:
        pass

print(f"n=5: {beta2_count} tournaments with beta_2 > 0 (out of {total})")

# n=6 exhaustive
n = 6
m = n * (n-1) // 2
total = 1 << m
beta2_count_6 = 0
t0 = time.time()

for bits in range(total):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    try:
        beta = path_betti_numbers(A, n, max_dim=3)
        if len(beta) > 2 and int(beta[2]) > 0:
            beta2_count_6 += 1
            H = H_tournament(A, n)
            beta_list = [int(beta[k]) for k in range(min(4, len(beta)))]
            if beta2_count_6 <= 10:
                print(f"  n=6, H={H}, beta={beta_list}")
    except:
        pass

    if (bits + 1) % 10000 == 0 and beta2_count_6 > 0:
        print(f"  ... {bits+1}/{total}, {beta2_count_6} with beta_2>0 so far")

print(f"n=6: {beta2_count_6} tournaments with beta_2 > 0 (out of {total}) in {time.time()-t0:.1f}s")

# ===== Deletion analysis for S-phase at n=6 =====
print("\n" + "=" * 70)
print("S-PHASE DELETIONS at n=6")
print("=" * 70)

n = 6
m = n * (n-1) // 2
total = 1 << m
s_phase_del_data = []

for bits in range(total):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    try:
        beta = path_betti_numbers(A, n, max_dim=4)
    except:
        continue

    if len(beta) > 3 and int(beta[3]) > 0:
        # S-phase — check all vertex deletions
        del_betti = []
        for v in range(n):
            # Delete vertex v
            remaining = [i for i in range(n) if i != v]
            sub_A = [[0]*5 for _ in range(5)]
            for i in range(5):
                for j in range(5):
                    sub_A[i][j] = A[remaining[i]][remaining[j]]
            try:
                sub_beta = path_betti_numbers(sub_A, 5, max_dim=3)
                sub_list = [int(sub_beta[k]) if k < len(sub_beta) else 0 for k in range(4)]
            except:
                sub_list = [1, 0, 0, 0]
            del_betti.append(tuple(sub_list))

        s_phase_del_data.append({
            'beta': [int(beta[k]) for k in range(min(5, len(beta)))],
            'del_betti': del_betti,
            'H': H_tournament(A, n)
        })

print(f"S-phase tournaments: {len(s_phase_del_data)}")

# Check: do all deletions have beta_2 > 0?
all_del_b2 = 0
any_del_b2 = 0
for d in s_phase_del_data:
    if all(b[2] > 0 for b in d['del_betti']):
        all_del_b2 += 1
    if any(b[2] > 0 for b in d['del_betti']):
        any_del_b2 += 1

print(f"ALL deletions have beta_2>0: {all_del_b2}")
print(f"ANY deletion has beta_2>0: {any_del_b2}")

# What DO the deletions look like?
del_patterns = Counter()
for d in s_phase_del_data:
    pattern = tuple(sorted(d['del_betti']))
    del_patterns[pattern] += 1

print(f"\nDeletion Betti patterns:")
for pat, cnt in del_patterns.most_common():
    print(f"  {list(pat)}: {cnt}")

# Show examples
for d in s_phase_del_data[:3]:
    print(f"\n  H={d['H']}, beta={d['beta']}")
    for i, db in enumerate(d['del_betti']):
        print(f"    del v{i}: beta={list(db)}")

# ===== C-phase deletions at n=6 =====
print("\n" + "=" * 70)
print("C-PHASE DELETIONS at n=6 (sample)")
print("=" * 70)

c_phase_del_data = []

for bits in range(total):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    H = H_tournament(A, n)
    if H != 45:  # Only H=45 C-phase maximizers
        continue

    try:
        beta = path_betti_numbers(A, n, max_dim=4)
    except:
        continue

    if len(beta) > 1 and int(beta[1]) > 0 and (len(beta) <= 3 or int(beta[3]) == 0):
        # C-phase maximizer
        del_betti = []
        for v in range(n):
            remaining = [i for i in range(n) if i != v]
            sub_A = [[0]*5 for _ in range(5)]
            for i in range(5):
                for j in range(5):
                    sub_A[i][j] = A[remaining[i]][remaining[j]]
            try:
                sub_beta = path_betti_numbers(sub_A, 5, max_dim=3)
                sub_list = [int(sub_beta[k]) if k < len(sub_beta) else 0 for k in range(4)]
            except:
                sub_list = [1, 0, 0, 0]
            del_betti.append(tuple(sub_list))

        c_phase_del_data.append({
            'beta': [int(beta[k]) for k in range(min(5, len(beta)))],
            'del_betti': del_betti,
            'H': H
        })

print(f"C-phase H=45 maximizers: {len(c_phase_del_data)}")

del_patterns_c = Counter()
for d in c_phase_del_data:
    pattern = tuple(sorted(d['del_betti']))
    del_patterns_c[pattern] += 1

print(f"\nDeletion Betti patterns:")
for pat, cnt in del_patterns_c.most_common():
    print(f"  {list(pat)}: {cnt}")

# Show examples
for d in c_phase_del_data[:2]:
    print(f"\n  H={d['H']}, beta={d['beta']}")
    for i, db in enumerate(d['del_betti']):
        print(f"    del v{i}: beta={list(db)}")
