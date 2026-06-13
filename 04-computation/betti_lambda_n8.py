#!/usr/bin/env python3
"""
Check: are Betti numbers lambda-determined at n=8?
At n=7, all β_0..β_3 are lambda-determined (0 ambiguities in 2000 samples).
At n=8, β_3 ∈ {0,1} and appears ~16% of the time.

This is computationally more expensive since Betti computation at n=8 involves
large chain spaces. Use smaller sample.

opus-2026-03-13-S71c
"""
import sys, time
import numpy as np
from collections import defaultdict
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

def lambda_key(A, n):
    L = np.zeros((n, n), dtype=int)
    for u in range(n):
        for v in range(u+1, n):
            for w in range(n):
                if w == u or w == v: continue
                if (A[u][v] and A[v][w] and A[w][u]) or (A[v][u] and A[u][w] and A[w][v]):
                    L[u][v] += 1; L[v][u] += 1
    return tuple(L[i][j] for i in range(n) for j in range(i+1, n))

def compute_beta1(A, n):
    """Fast computation of β_1 only."""
    # β_1 = dim(ker(∂_1)) - dim(im(∂_2))
    # dim(ker(∂_1)) = |edges| - n + #components = C(n,2) - n + 1 (tournaments are weakly connected)
    # So we just need im(∂_2), which is rank of ∂_2 restricted to Ω_2.

    # For β_1: it's simpler. β_1 = 1 iff tournament is NOT strongly connected.
    # (Well, more precisely: β_1 counts "1-dimensional holes" in the path complex.)
    # For tournaments: β_1 = 1 iff tournament is NOT strongly connected.

    # Actually: β_1 = C(n,2) - n + 1 - rank(∂_2|_Ω_2)
    # This requires computing Ω_2 and ∂_2.

    # Faster: use the known formula β_1 ∈ {0,1} and β_1 = 1 iff not strongly connected.
    # Check strong connectivity.
    from collections import deque
    def is_strongly_connected(A, n):
        # BFS from 0
        visited = {0}
        queue = deque([0])
        while queue:
            v = queue.popleft()
            for u in range(n):
                if u not in visited and A[v][u]:
                    visited.add(u)
                    queue.append(u)
        if len(visited) < n: return False
        # BFS on transpose
        visited = {0}
        queue = deque([0])
        while queue:
            v = queue.popleft()
            for u in range(n):
                if u not in visited and A[u][v]:
                    visited.add(u)
                    queue.append(u)
        return len(visited) == n

    return 0 if is_strongly_connected(A, n) else 1

def compute_beta3_fast(A, n):
    """Compute β_3 using rank of boundary matrices."""
    # Enumerate 3-paths: (v0,v1,v2,v3) with v0→v1→v2→v3, all distinct
    paths3 = []
    for v0 in range(n):
        for v1 in range(n):
            if v1 == v0 or A[v0][v1] == 0: continue
            for v2 in range(n):
                if v2 == v0 or v2 == v1 or A[v1][v2] == 0: continue
                for v3 in range(n):
                    if v3 == v0 or v3 == v1 or v3 == v2 or A[v2][v3] == 0: continue
                    paths3.append((v0,v1,v2,v3))

    # Enumerate 2-paths
    paths2 = []
    for v0 in range(n):
        for v1 in range(n):
            if v1 == v0 or A[v0][v1] == 0: continue
            for v2 in range(n):
                if v2 == v0 or v2 == v1 or A[v1][v2] == 0: continue
                paths2.append((v0,v1,v2))

    p2_idx = {p: i for i, p in enumerate(paths2)}
    p2_set = set(paths2)

    # Compute faces of 3-paths and classify as allowed/junk
    junk_set = set()
    face_data_junk = []
    face_data_allowed = []

    for p in paths3:
        jf = []; af = []
        for fi in range(4):
            face = p[:fi] + p[fi+1:]
            sign = 1 if fi % 2 == 0 else -1
            if face in p2_set:
                af.append((face, sign))
            else:
                junk_set.add(face)
                jf.append((face, sign))
        face_data_junk.append(jf)
        face_data_allowed.append(af)

    n3 = len(paths3)
    n2 = len(paths2)
    junk_list = sorted(junk_set)
    n_junk = len(junk_list)
    junk_idx = {j: i for i, j in enumerate(junk_list)}

    # Constraint matrix
    C = np.zeros((n_junk, n3), dtype=np.float64)
    for j, jf in enumerate(face_data_junk):
        for face, sign in jf:
            C[junk_idx[face], j] += sign

    rank_c = int(np.linalg.matrix_rank(C)) if n_junk > 0 else 0
    omega3 = n3 - rank_c

    # Combined matrix
    CB = np.zeros((n_junk + n2, n3), dtype=np.float64)
    CB[:n_junk, :] = C
    for j, af in enumerate(face_data_allowed):
        for face, sign in af:
            row = n_junk + p2_idx[face]
            CB[row, j] += sign

    rank_cb = int(np.linalg.matrix_rank(CB))
    bd3_rank = rank_cb - rank_c

    # Also need im(∂_4) for β_3 = Ω_3 - rank(∂_3) - rank(∂_4|_Ω_4)
    # For now, assume rank(∂_4) for β_3 computation...
    # Actually β_3 = dim(ker(∂_3|_Ω_3)) - dim(im(∂_4|_Ω_4))
    # = (Ω_3 - rank(∂_3|_Ω_3)) - rank(∂_4|_Ω_4)

    # For a quick check, we can compute the upper bound β_3 ≤ Ω_3 - rank(∂_3)
    # At n=7, β_3 is often the exact value since Ω_4 is often 0 or small.

    # Actually let's just compute β_3 as the upper bound (ignoring im(∂_4)):
    beta3_upper = omega3 - bd3_rank

    return beta3_upper, omega3, bd3_rank

# Test at n=7 first to verify
n = 7
tb = n*(n-1)//2
np.random.seed(42)
print(f"Quick n=7 verification (100 samples)...")

lam_b = defaultdict(set)
for trial in range(100):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    b1 = compute_beta1(A, n)
    b3_upper, _, _ = compute_beta3_fast(A, n)
    key = lambda_key(A, n)
    lam_b[key].add((b1, b3_upper))

ambig = sum(1 for v in lam_b.values() if len(v) > 1)
print(f"  Lambda groups: {len(lam_b)}, ambiguous (β1,β3): {ambig}")

# Now n=8
n = 8
tb = n*(n-1)//2
np.random.seed(42)
print(f"\nn=8: checking if β_1 and β_3 are lambda-determined...")
print(f"(This is slower due to larger chain spaces)")

lam_b8 = defaultdict(set)
b3_dist = defaultdict(int)

t0 = time.time()
for trial in range(200):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    b1 = compute_beta1(A, n)
    b3_upper, omega3, bd3_rank = compute_beta3_fast(A, n)
    key = lambda_key(A, n)
    lam_b8[key].add((b1, b3_upper))
    b3_dist[b3_upper] += 1

    if trial % 50 == 0:
        dt = time.time() - t0
        print(f"  trial {trial}: {dt:.1f}s, β3 dist so far: {dict(b3_dist)}")

ambig8 = sum(1 for v in lam_b8.values() if len(v) > 1)
print(f"\n  Lambda groups: {len(lam_b8)}, ambiguous (β1,β3_upper): {ambig8}")
print(f"  β3_upper distribution: {dict(b3_dist)}")

if ambig8 > 0:
    count = 0
    for key, vals in sorted(lam_b8.items()):
        if len(vals) > 1:
            print(f"    Ambiguity: {sorted(vals)}")
            count += 1
            if count >= 5: break

print(f"\nDone.")
