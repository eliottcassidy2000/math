#!/usr/bin/env python3
"""
betti_alpha_v2.py — opus-2026-03-13-S71c

Simplified: check β₃ vs c3/c5/H/α₁ at n=7.
Use the existing betti computation from betti_sigma_hierarchy.py.
"""

import sys, time
import numpy as np
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

def compute_betti(A, n, max_d=3):
    """Compute Betti numbers through max_d."""
    paths = {0: [(v,) for v in range(n)]}
    for d in range(1, max_d + 2):
        pd = []
        for p in paths[d-1]:
            last = p[-1]
            used = set(p)
            for v in range(n):
                if v not in used and A[last][v] == 1:
                    pd.append(p + (v,))
        paths[d] = pd
        if not pd:
            break

    actual_max = max(d for d in paths if paths[d])

    path_idx = {}
    for d in range(actual_max + 1):
        path_idx[d] = {p: i for i, p in enumerate(paths[d])}

    allowed_set = {d: set(paths[d]) for d in range(actual_max + 1)}

    omega_dims = [n]
    bd_ranks = [0]

    for d in range(1, min(max_d + 1, actual_max + 1)):
        n_d = len(paths[d])
        n_dm1 = len(paths[d-1])
        if n_d == 0:
            omega_dims.append(0)
            bd_ranks.append(0)
            continue

        junk_set = set()
        face_junk = []
        face_allowed = []

        for p in paths[d]:
            jf = []; af = []
            for fi in range(d + 1):
                face = p[:fi] + p[fi+1:]
                sign = 1 if fi % 2 == 0 else -1
                if face in allowed_set[d-1]:
                    af.append((face, sign))
                else:
                    junk_set.add(face)
                    jf.append((face, sign))
            face_junk.append(jf)
            face_allowed.append(af)

        junk_list = sorted(junk_set)
        n_junk = len(junk_list)
        junk_idx = {j: i for i, j in enumerate(junk_list)}

        C = np.zeros((n_junk, n_d), dtype=float)
        for j, jf in enumerate(face_junk):
            for face, sign in jf:
                C[junk_idx[face], j] += sign

        rank_c = int(np.linalg.matrix_rank(C)) if n_junk > 0 else 0
        omega_d = n_d - rank_c
        omega_dims.append(omega_d)

        CB = np.zeros((n_junk + n_dm1, n_d), dtype=float)
        CB[:n_junk, :] = C
        for j, af in enumerate(face_allowed):
            for face, sign in af:
                row = n_junk + path_idx[d-1][face]
                CB[row, j] += sign

        rank_cb = int(np.linalg.matrix_rank(CB))
        bd_rank = rank_cb - rank_c
        bd_ranks.append(bd_rank)

    betti = []
    for d in range(min(max_d + 1, len(omega_dims))):
        od = omega_dims[d]
        rd = bd_ranks[d] if d < len(bd_ranks) else 0
        rd1 = bd_ranks[d+1] if d+1 < len(bd_ranks) else 0
        betti.append(od - rd - rd1)

    return betti

n = 7
tb = n*(n-1)//2
np.random.seed(42)

print(f"n={n}: Betti numbers vs cycle/H structure, 300 samples...")
results = []
t0 = time.time()

for trial in range(300):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)

    betti = compute_betti(A, n, max_d=3)
    b0, b1 = betti[0], betti[1]
    b2 = betti[2] if len(betti) > 2 else 0
    b3 = betti[3] if len(betti) > 3 else 0

    H = sum(1 for v in range(n)
            for mask_end in range(n)
            if True)  # placeholder
    # Compute H properly
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for ms in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != ms: continue
            for v in range(n):
                if not (mask & (1 << v)): continue
                pm = mask ^ (1 << v)
                t = 0
                for u in range(n):
                    if (pm & (1 << u)) and A[u][v]:
                        t += dp.get((pm, u), 0)
                if t:
                    dp[(mask, v)] = t
    H = sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))

    c3 = int(np.trace(A @ A @ A)) // 3
    c5 = int(np.trace(np.linalg.matrix_power(A, 5))) // 5

    results.append((b0, b1, b2, b3, H, c3, c5))

    if trial % 100 == 0:
        dt = time.time() - t0
        print(f"  trial {trial}: {dt:.1f}s, β=[{b0},{b1},{b2},{b3}]")

dt = time.time() - t0
print(f"Done: {dt:.1f}s")

# Analyze
phase_counts = defaultdict(int)
for r in results:
    phase_counts[(r[1], r[3])] += 1

print(f"\nPhase (β₁,β₃) distribution:")
for (b1,b3), cnt in sorted(phase_counts.items()):
    print(f"  (β₁={b1}, β₃={b3}): {cnt} ({100*cnt/len(results):.1f}%)")

# β₃ statistics
b3_0 = [r for r in results if r[1]==0 and r[3]==0]
b3_1 = [r for r in results if r[1]==0 and r[3]==1]
b1_1 = [r for r in results if r[1]==1]

for label, group in [("β₁=0,β₃=0", b3_0), ("β₁=0,β₃=1", b3_1), ("β₁=1", b1_1)]:
    if group:
        H_vals = [r[4] for r in group]
        c3_vals = [r[5] for r in group]
        c5_vals = [r[6] for r in group]
        print(f"\n{label} (n={len(group)}):")
        print(f"  H: mean={np.mean(H_vals):.1f}, range=[{min(H_vals)},{max(H_vals)}]")
        print(f"  c3: mean={np.mean(c3_vals):.1f}, range=[{min(c3_vals)},{max(c3_vals)}]")
        print(f"  c5: mean={np.mean(c5_vals):.1f}, range=[{min(c5_vals)},{max(c5_vals)}]")

# Check c3 overlap
if b3_0 and b3_1:
    c3_set_0 = set(r[5] for r in b3_0)
    c3_set_1 = set(r[5] for r in b3_1)
    overlap = c3_set_0 & c3_set_1
    print(f"\nc3 overlap between β₃=0 and β₃=1: {sorted(overlap) if overlap else 'NONE'}")

    H_set_0 = set(r[4] for r in b3_0)
    H_set_1 = set(r[4] for r in b3_1)
    overlap_H = H_set_0 & H_set_1
    print(f"H overlap: {len(overlap_H)} values, min={min(overlap_H) if overlap_H else 'N/A'}, max={max(overlap_H) if overlap_H else 'N/A'}")

print("\nDone.")
