#!/usr/bin/env python3
"""
n9_maximizer_betti.py — Find the actual n=9 maximizer (H=3357) and compute its Betti

OEIS A038375: max H(T) for n vertices = 1,1,3,5,15,45,189,661,3357,15745,95095
At n=9, H_max = 3357. The beta5_n9_check only found H=3255 in 10k samples.
We need many more samples, or use the known maximizer structure.

Strategy: Sample more aggressively, and also try SC tournaments with near-regular scores.

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, random
import numpy as np
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

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

def random_regular_tournament(n):
    """Generate a random regular tournament on n vertices (n must be odd).
    Uses the Wilson method: start with any tournament, improve by arc flips."""
    assert n % 2 == 1
    target = (n - 1) // 2

    # Start with random tournament
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    # Iteratively improve scores toward regular
    for _ in range(1000):
        scores = [sum(A[i]) for i in range(n)]
        max_score_diff = max(abs(s - target) for s in scores)
        if max_score_diff == 0:
            break
        # Find vertex with too high score and one with too low
        high_v = max(range(n), key=lambda v: scores[v])
        low_v = min(range(n), key=lambda v: scores[v])
        if high_v == low_v:
            break
        # If high_v beats low_v, flip the arc
        if A[high_v][low_v]:
            A[high_v][low_v] = 0
            A[low_v][high_v] = 1
        else:
            # Find w that high_v beats and low_v loses to, flip both
            candidates = [w for w in range(n) if w != high_v and w != low_v
                         and A[high_v][w] and A[w][low_v]]
            if candidates:
                w = random.choice(candidates)
                A[high_v][w] = 0
                A[w][high_v] = 1
                A[w][low_v] = 0
                A[low_v][w] = 1

    return A

# ===== Phase 1: Massive sampling for H=3357 =====
print("=" * 70)
print("SEARCHING FOR n=9 MAXIMIZER (H=3357)")
print("=" * 70)

n = 9
m = n * (n-1) // 2  # 36
t0 = time.time()
best_H = 0
best_tours = []

# Try random regular tournaments (maximizer is always regular at odd n)
print("\nPhase 1a: Random regular tournaments...")
for trial in range(100000):
    A = random_regular_tournament(n)
    scores = [sum(A[i]) for i in range(n)]
    if max(scores) - min(scores) > 0:
        continue  # Not actually regular
    H = H_tournament(A, n)
    if H > best_H:
        best_H = H
        best_tours = [(H, [row[:] for row in A])]
        print(f"  trial {trial}: NEW MAX H={H}")
    elif H == best_H:
        best_tours.append((H, [row[:] for row in A]))

    if trial % 10000 == 0:
        elapsed = time.time() - t0
        print(f"  {trial}/100000 ({elapsed:.1f}s), best H={best_H}, found {len(best_tours)} maximizers")

print(f"\nPhase 1a done: best H={best_H}, {len(best_tours)} maximizers found")

# Phase 1b: Random general tournaments
print("\nPhase 1b: Random general tournaments...")
for trial in range(200000):
    bits = random.randint(0, (1 << m) - 1)
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
    if H > best_H:
        best_H = H
        best_tours = [(H, [row[:] for row in A])]
        print(f"  trial {trial}: NEW MAX H={H}")
    elif H == best_H:
        best_tours.append((H, [row[:] for row in A]))

    if trial % 50000 == 0:
        elapsed = time.time() - t0
        print(f"  {trial}/200000 ({elapsed:.1f}s), best H={best_H}, found {len(best_tours)} maximizers")

print(f"\nPhase 1b done: best H={best_H}, {len(best_tours)} maximizers found")

# ===== Report best tournaments =====
print("\n" + "=" * 70)
print(f"TOP TOURNAMENTS (H={best_H})")
print("=" * 70)

# Deduplicate and show first 5
seen = set()
unique_tours = []
for H, A in best_tours:
    key = tuple(tuple(row) for row in A)
    if key not in seen:
        seen.add(key)
        unique_tours.append((H, A))

print(f"Unique maximizers: {len(unique_tours)}")
for idx, (H, A) in enumerate(unique_tours[:5]):
    scores = sorted([sum(A[i]) for i in range(n)])
    # Compute c3
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (A[i][j] and A[j][k] and A[k][i]) or (A[i][k] and A[k][j] and A[j][i]):
                    c3 += 1

    # Skew eigenvalues
    S = np.array([[A[i][j] - A[j][i] for j in range(n)] for i in range(n)], dtype=float)
    evals_raw = np.linalg.eigvals(S)
    pos_imag = sorted([e.imag for e in evals_raw if e.imag > 0.01])
    gap = max(pos_imag) - min(pos_imag) if len(pos_imag) >= 2 else 0

    print(f"\n  [{idx+1}] H={H}, score={tuple(scores)}, c3={c3}, gap={gap:.4f}")
    print(f"      eigenvalues: [{', '.join(f'{e:.4f}' for e in pos_imag)}]")
    print(f"      Adjacency row sums: {[sum(A[i]) for i in range(n)]}")

print(f"\nTotal time: {time.time()-t0:.1f}s")
