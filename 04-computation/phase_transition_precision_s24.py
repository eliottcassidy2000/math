#!/usr/bin/env python3
"""
phase_transition_precision_s24.py — opus-2026-04-05-S24

Precision analysis of the order-disorder phase transition under random arc flips.

Key question: Starting from the transitive tournament, how many random arc
flips does it take for H(T) to reach [some fraction] of its equilibrium value?

The initial experiment found the midpoint at ~0.33*m for n=5,6,7.
Is this exactly 1/3? Does it depend on n? What's the functional form?

We also measure:
- The EXACT equilibrium mean H (= n!/2^{n-1})
- The relaxation curve H(step) / E[H_random]
- The autocorrelation time
- The variance of H along the trajectory
"""

import numpy as np
from collections import defaultdict
import time
import random
import math
import sys

sys.stdout.reconfigure(encoding='utf-8')

def bits_to_adj_fast(bits, n):
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

def count_ham_paths_fast(A, n):
    m = 1 << n
    dp = [[0]*n for _ in range(m)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, m):
        for v in range(n):
            if dp[mask][v] == 0:
                continue
            if not (mask & (1 << v)):
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v][u] == 1:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = m - 1
    return sum(dp[full][v] for v in range(n))

def count_3cycles_fast(A, n):
    c = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j] and A[j][k] and A[k][i]:
                    c += 1
                if A[i][k] and A[k][j] and A[j][i]:
                    c += 1
    return c

def transitive_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            A[i][j] = 1
    return A

def flip_arc_inplace(A, i, j):
    A[i][j] = 1 - A[i][j]
    A[j][i] = 1 - A[j][i]

print("=" * 70)
print("PHASE TRANSITION PRECISION ANALYSIS")
print("=" * 70)

results = {}

for n in [5, 6, 7, 8]:
    m_arcs = n * (n - 1) // 2
    mean_H_theory = math.factorial(n) / 2**(n-1)
    max_flips = 5 * m_arcs  # generous

    num_trials = {5: 2000, 6: 1000, 7: 200, 8: 50}[n]

    print(f"\n{'='*60}")
    print(f"n = {n}, m = {m_arcs}, E[H_random] = {mean_H_theory:.2f}, trials = {num_trials}")
    print(f"{'='*60}")

    # Record H at every step for precise analysis
    H_by_step = defaultdict(list)

    t0 = time.time()
    for trial in range(num_trials):
        A = transitive_tournament(n)
        for step in range(max_flips + 1):
            # Record at key steps
            if step <= 2 * m_arcs or step % max(1, m_arcs // 5) == 0:
                H = count_ham_paths_fast(A, n)
                H_by_step[step].append(H)

            if step < max_flips:
                i = random.randint(0, n-2)
                j = random.randint(i+1, n-1)
                flip_arc_inplace(A, i, j)

    elapsed = time.time() - t0
    print(f"  Time: {elapsed:.1f}s")

    # Compute mean H at each step
    steps_sorted = sorted(H_by_step.keys())
    mean_H = {s: np.mean(H_by_step[s]) for s in steps_sorted}

    # Find various crossing points
    crossings = {}
    for frac in [0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99]:
        target = 1 + frac * (mean_H_theory - 1)
        for i, s in enumerate(steps_sorted):
            if mean_H[s] >= target:
                crossings[frac] = s
                break

    print(f"\n  Crossing analysis (fraction of equilibrium reached):")
    print(f"  {'Fraction':>8} | {'Step':>5} | {'Step/m':>7} | {'Step/m^2':>8}")
    print(f"  {'-'*8}-|-{'-'*5}-|-{'-'*7}-|-{'-'*8}")
    for frac in sorted(crossings.keys()):
        s = crossings[frac]
        print(f"  {frac:8.2f} | {s:5d} | {s/m_arcs:7.3f} | {s/m_arcs**2:8.4f}")

    # The 50% crossing in terms of m
    if 0.5 in crossings:
        t_half = crossings[0.5]
        print(f"\n  MIDPOINT: step {t_half} = {t_half/m_arcs:.4f} * m")

    # Exponential fit: H(t) ≈ H_eq * (1 - exp(-t/tau))
    # Find tau by fitting
    data_steps = [s for s in steps_sorted if s > 0 and mean_H[s] > 1.1]
    if len(data_steps) > 5:
        # Linear regression on log(H_eq - H(t)) vs t
        log_data = []
        for s in data_steps:
            diff = mean_H_theory - mean_H[s]
            if diff > 0.1:
                log_data.append((s, np.log(diff)))

        if len(log_data) > 5:
            x = np.array([d[0] for d in log_data[:len(log_data)//2]])  # use first half
            y = np.array([d[1] for d in log_data[:len(log_data)//2]])
            slope, intercept = np.polyfit(x, y, 1)
            tau = -1.0 / slope
            print(f"\n  Exponential relaxation: tau = {tau:.3f} = {tau/m_arcs:.4f} * m")
            print(f"  Relaxation time in units of m: {tau/m_arcs:.4f}")

    # Print detailed trajectory
    print(f"\n  Detailed trajectory (step, mean_H, std_H, normalized):")
    for s in steps_sorted[:40]:
        vals = H_by_step[s]
        mH = np.mean(vals)
        sH = np.std(vals)
        norm = (mH - 1) / (mean_H_theory - 1) if mean_H_theory > 1 else 0
        print(f"    step {s:4d}: mean_H={mH:8.1f}, std={sH:7.1f}, normalized={norm:.4f}")

    results[n] = {
        'crossings': crossings,
        'm': m_arcs,
        'mean_H_theory': mean_H_theory,
    }

# Cross-n comparison
print("\n\n" + "=" * 70)
print("CROSS-n COMPARISON")
print("=" * 70)
print(f"\n  {'n':>3} | {'m':>4} | {'E[H]':>8} | {'t_50%':>6} | {'t_50%/m':>8} | {'tau/m':>7}")
print(f"  {'-'*3}-|-{'-'*4}-|-{'-'*8}-|-{'-'*6}-|-{'-'*8}-|-{'-'*7}")
for n in sorted(results.keys()):
    r = results[n]
    t50 = r['crossings'].get(0.5, -1)
    t50m = t50 / r['m'] if t50 > 0 else -1
    print(f"  {n:3d} | {r['m']:4d} | {r['mean_H_theory']:8.1f} | {t50:6d} | {t50m:8.4f} | ")

# Also compare the "10% crossing" and "90% crossing" for universality
print(f"\n  Crossing at 10%:")
for n in sorted(results.keys()):
    r = results[n]
    t10 = r['crossings'].get(0.1, -1)
    if t10 > 0:
        print(f"    n={n}: step {t10} = {t10/r['m']:.4f} * m")

print(f"\n  Crossing at 90%:")
for n in sorted(results.keys()):
    r = results[n]
    t90 = r['crossings'].get(0.9, -1)
    if t90 > 0:
        print(f"    n={n}: step {t90} = {t90/r['m']:.4f} * m")
