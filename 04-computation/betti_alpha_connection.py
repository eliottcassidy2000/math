#!/usr/bin/env python3
"""
betti_alpha_connection.py — opus-2026-03-13-S71c

QUESTION: How do Betti numbers relate to the α-hierarchy?
- β_1 ∈ {0,1}: β_1=1 iff not strongly connected
- β_3 ∈ {0,1}: β_3=1 for ~8% at n=7
- β_1 · β_3 = 0 (seesaw)
- Both are lambda-determined

Does β_3 correlate with α₁, α₂, or specific cycle counts?
Does the "phase" (β_1=0,β_3=0 vs β_1=0,β_3=1) correspond to
a threshold in the independence polynomial?
"""

import sys, time
import numpy as np
from itertools import combinations
from collections import defaultdict, deque
sys.stdout.reconfigure(line_buffering=True)

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def is_strongly_connected(A, n):
    visited = {0}
    queue = deque([0])
    while queue:
        v = queue.popleft()
        for u in range(n):
            if u not in visited and A[v][u]:
                visited.add(u); queue.append(u)
    if len(visited) < n: return False
    visited = {0}
    queue = deque([0])
    while queue:
        v = queue.popleft()
        for u in range(n):
            if u not in visited and A[u][v]:
                visited.add(u); queue.append(u)
    return len(visited) == n

def compute_betti_3(A, n):
    """Compute β_3 (upper bound, ignoring im(∂_4) which is 0 at n≤7)."""
    # Enumerate 3-paths
    paths3 = []
    for v0 in range(n):
        for v1 in range(n):
            if v1==v0 or not A[v0][v1]: continue
            for v2 in range(n):
                if v2==v0 or v2==v1 or not A[v1][v2]: continue
                for v3 in range(n):
                    if v3==v0 or v3==v1 or v3==v2 or not A[v2][v3]: continue
                    paths3.append((v0,v1,v2,v3))

    paths2 = []
    for v0 in range(n):
        for v1 in range(n):
            if v1==v0 or not A[v0][v1]: continue
            for v2 in range(n):
                if v2==v0 or v2==v1 or not A[v1][v2]: continue
                paths2.append((v0,v1,v2))

    p2_set = set(paths2)
    p2_idx = {p: i for i, p in enumerate(paths2)}

    junk_set = set()
    face_junk = []
    face_allowed = []

    for p in paths3:
        jf = []; af = []
        for fi in range(4):
            face = p[:fi] + p[fi+1:]
            sign = 1 if fi%2==0 else -1
            if face in p2_set:
                af.append((face, sign))
            else:
                junk_set.add(face)
                jf.append((face, sign))
        face_junk.append(jf)
        face_allowed.append(af)

    n3 = len(paths3)
    n2 = len(paths2)
    junk_list = sorted(junk_set)
    n_junk = len(junk_list)
    junk_idx = {j: i for i, j in enumerate(junk_list)}

    C = np.zeros((n_junk, n3), dtype=np.float64)
    for j, jf in enumerate(face_junk):
        for face, sign in jf:
            C[junk_idx[face], j] += sign

    rank_c = int(np.linalg.matrix_rank(C)) if n_junk > 0 else 0
    omega3 = n3 - rank_c

    CB = np.zeros((n_junk + n2, n3), dtype=np.float64)
    CB[:n_junk, :] = C
    for j, af in enumerate(face_allowed):
        for face, sign in af:
            row = n_junk + p2_idx[face]
            CB[row, j] += sign

    rank_cb = int(np.linalg.matrix_rank(CB))
    bd3_rank = rank_cb - rank_c

    beta3 = omega3 - bd3_rank
    return beta3

def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for ms in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != ms:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                pm = mask ^ (1 << v)
                t = 0
                for u in range(n):
                    if (pm & (1 << u)) and A[u][v]:
                        t += dp.get((pm, u), 0)
                if t:
                    dp[(mask, v)] = t
    return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))

n = 7
tb = n*(n-1)//2
np.random.seed(42)

print(f"n={n}: Computing β_3, H, c3, c5, α₁ for 500 samples...")
results = []
t0 = time.time()

for trial in range(500):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)

    sc = is_strongly_connected(A, n)
    beta1 = 0 if sc else 1
    H = count_ham_paths(A, n)
    c3 = int(np.trace(A @ A @ A)) // 3
    c5 = int(np.trace(np.linalg.matrix_power(A, 5))) // 5

    if beta1 == 0:
        beta3 = compute_betti_3(A, n)
    else:
        beta3 = 0  # seesaw: β₁β₃=0

    # Find α₂ quickly
    cycles3 = []
    for a in range(n):
        for b in range(a+1, n):
            for c in range(b+1, n):
                if (A[a][b] and A[b][c] and A[c][a]) or (A[a][c] and A[c][b] and A[b][a]):
                    cycles3.append(frozenset([a,b,c]))
    cycles3 = list(set(cycles3))
    a2 = sum(1 for i in range(len(cycles3)) for j in range(i+1, len(cycles3)) if not (cycles3[i] & cycles3[j]))
    a1 = c3 + c5 + (H - 1 - 4*a2) // 2 - c3 - c5  # This gives c7
    a1_total = c3 + c5 + a1  # a1 = c7 here, so a1_total = c3+c5+c7 = total α₁

    results.append((beta1, beta3, H, c3, c5, a1, a2, a1_total))

    if trial % 100 == 0:
        dt = time.time() - t0
        print(f"  trial {trial}: {dt:.1f}s")

dt = time.time() - t0
print(f"Done: {dt:.1f}s")

# Analyze β₃ vs various quantities
b3_0 = [r for r in results if r[0]==0 and r[1]==0]  # SC, β₃=0
b3_1 = [r for r in results if r[0]==0 and r[1]==1]  # SC, β₃=1
b1_1 = [r for r in results if r[0]==1]               # NSC, β₁=1

print(f"\nPhase distribution:")
print(f"  SC, β₃=0: {len(b3_0)} ({100*len(b3_0)/len(results):.1f}%)")
print(f"  SC, β₃=1: {len(b3_1)} ({100*len(b3_1)/len(results):.1f}%)")
print(f"  NSC, β₁=1: {len(b1_1)} ({100*len(b1_1)/len(results):.1f}%)")

print(f"\nH statistics by phase:")
for label, group in [("SC β₃=0", b3_0), ("SC β₃=1", b3_1), ("NSC β₁=1", b1_1)]:
    if group:
        H_vals = [r[2] for r in group]
        print(f"  {label}: H mean={np.mean(H_vals):.1f}, range=[{min(H_vals)}, {max(H_vals)}]")

print(f"\nc3 statistics by phase:")
for label, group in [("SC β₃=0", b3_0), ("SC β₃=1", b3_1), ("NSC β₁=1", b1_1)]:
    if group:
        c3_vals = [r[3] for r in group]
        print(f"  {label}: c3 mean={np.mean(c3_vals):.1f}, range=[{min(c3_vals)}, {max(c3_vals)}]")

print(f"\nc5 statistics by phase:")
for label, group in [("SC β₃=0", b3_0), ("SC β₃=1", b3_1), ("NSC β₁=1", b1_1)]:
    if group:
        c5_vals = [r[4] for r in group]
        print(f"  {label}: c5 mean={np.mean(c5_vals):.1f}, range=[{min(c5_vals)}, {max(c5_vals)}]")

print(f"\nc7 (= α₁ - c3 - c5) statistics by phase:")
for label, group in [("SC β₃=0", b3_0), ("SC β₃=1", b3_1), ("NSC β₁=1", b1_1)]:
    if group:
        c7_vals = [r[5] for r in group]
        print(f"  {label}: c7 mean={np.mean(c7_vals):.1f}, range=[{min(c7_vals)}, {max(c7_vals)}]")

print(f"\nα₂ (disjoint 3-cycle pairs) statistics by phase:")
for label, group in [("SC β₃=0", b3_0), ("SC β₃=1", b3_1), ("NSC β₁=1", b1_1)]:
    if group:
        a2_vals = [r[6] for r in group]
        print(f"  {label}: α₂ mean={np.mean(a2_vals):.1f}, range=[{min(a2_vals)}, {max(a2_vals)}]")

# Is β₃ a threshold function of α₁ or c3?
print(f"\nβ₃=1 iff c3 ≥ threshold?")
c3_b3_0 = set(r[3] for r in b3_0)
c3_b3_1 = set(r[3] for r in b3_1)
print(f"  c3 values with β₃=0: {sorted(c3_b3_0)}")
print(f"  c3 values with β₃=1: {sorted(c3_b3_1)}")
overlap = c3_b3_0 & c3_b3_1
print(f"  Overlap: {sorted(overlap) if overlap else 'NONE'}")
if overlap:
    print(f"  β₃ is NOT a threshold in c3 (overlap exists)")
else:
    print(f"  β₃ MIGHT be a threshold in c3!")

# H threshold?
H_b3_0 = set(r[2] for r in b3_0)
H_b3_1 = set(r[2] for r in b3_1)
overlap_H = H_b3_0 & H_b3_1
print(f"\nβ₃=1 iff H ≥ threshold?")
print(f"  H overlap: {len(overlap_H)} values")
if overlap_H:
    print(f"  NOT a simple threshold (overlap: min={min(overlap_H)}, max={max(overlap_H)})")

print("\nDone.")
