#!/usr/bin/env python3
"""
DEEP EXPLORATION: The boundary between real-rooted and non-real-rooted
I(Omega(T), x) at n=9.

Key questions:
1. What graph-theoretic property of Omega(T) determines failure?
2. Is there a spectral signature in the TOURNAMENT adjacency matrix?
3. Can we interpolate from the counterexample to a "generic" tournament?

Strategy: Start from the counterexample, flip arcs one at a time
toward a regular tournament, and track when real-rootedness is restored.

Author: opus-2026-03-06-S19
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import random_tournament
from itertools import combinations, permutations
import numpy as np
from collections import Counter

T_ce = [
    [0, 1, 0, 1, 0, 0, 1, 1, 0],
    [0, 0, 0, 1, 0, 0, 0, 0, 0],
    [1, 1, 0, 0, 1, 1, 1, 1, 0],
    [0, 0, 1, 0, 0, 1, 0, 1, 0],
    [1, 1, 0, 1, 0, 0, 0, 1, 0],
    [1, 1, 0, 0, 1, 0, 1, 1, 1],
    [0, 1, 0, 1, 1, 0, 0, 1, 0],
    [0, 1, 0, 0, 0, 0, 0, 0, 0],
    [1, 1, 1, 1, 1, 0, 1, 1, 0],
]
n = 9

def find_all_directed_odd_cycles(T, n):
    cycles = set()
    for length in range(3, n + 1, 2):
        for subset in combinations(range(n), length):
            for perm in permutations(subset[1:]):
                cycle = (subset[0],) + perm
                if all(T[cycle[i]][cycle[(i+1) % length]] for i in range(length)):
                    min_idx = cycle.index(min(cycle))
                    rotated = cycle[min_idx:] + cycle[:min_idx]
                    cycles.add(rotated)
    return list(cycles)

def independence_poly_dc(adj_dict):
    memo = {}
    def solve(verts):
        if verts in memo:
            return memo[verts]
        if not verts:
            return [1]
        v = max(verts, key=lambda u: len(adj_dict[u] & verts))
        p1 = solve(verts - {v})
        p2 = solve(verts - (adj_dict[v] & verts) - {v})
        maxlen = max(len(p1), len(p2) + 1)
        result = [0] * maxlen
        for i in range(len(p1)):
            result[i] += p1[i]
        for i in range(len(p2)):
            result[i + 1] += p2[i]
        memo[verts] = result
        return result
    return solve(frozenset(adj_dict.keys()))

def analyze_tournament(T, n):
    """Return (coeffs, newton_fail, disc, scores, spectral_data)."""
    cycles = find_all_directed_odd_cycles(T, n)
    m = len(cycles)
    if m == 0:
        return [1], False, None, tuple(sorted(sum(row) for row in T)), {}

    cvsets = [frozenset(c) for c in cycles]
    adj = {}
    for i in range(m):
        adj[i] = frozenset(j for j in range(m) if j != i and cvsets[i] & cvsets[j])
    coeffs = independence_poly_dc(adj)
    deg = len(coeffs) - 1
    while deg > 0 and coeffs[deg] == 0:
        deg -= 1

    newton_fail = False
    if deg >= 2:
        for k in range(1, deg):
            if coeffs[k]**2 < coeffs[k-1] * coeffs[k+1] * (k+1) / k:
                newton_fail = True
                break

    # Score and spectral data
    scores = tuple(sorted(sum(row) for row in T))
    A = np.array(T, dtype=float)
    eigs = np.linalg.eigvals(A)
    spectral_radius = max(abs(e) for e in eigs)
    aat_eigs = sorted(np.linalg.eigvals(A @ A.T).real, reverse=True)

    # Cycle length distribution
    len_dist = Counter(len(c) for c in cycles)

    # Omega density
    if m > 1:
        edges = sum(len(adj[i]) for i in range(m)) // 2
        density = 2 * edges / (m * (m - 1))
    else:
        density = 0.0

    # Independence number (from poly degree)
    alpha = deg

    spectral_data = {
        'spectral_radius': spectral_radius,
        'aat_lambda1': aat_eigs[0],
        'num_cycles': m,
        'len_dist': dict(len_dist),
        'omega_density': density,
        'alpha': alpha,
    }

    return coeffs[:deg+1], newton_fail, None, scores, spectral_data

print("=" * 70)
print("FAILURE BOUNDARY ANALYSIS AT n=9")
print("=" * 70)

# Part 1: Analyze the counterexample
print("\n--- Part 1: Counterexample analysis ---")
ce_coeffs, ce_fail, _, ce_scores, ce_spec = analyze_tournament(T_ce, n)
print(f"  Scores: {ce_scores}")
print(f"  I(Omega,x) coeffs: {ce_coeffs}")
print(f"  Newton fails: {ce_fail}")
print(f"  Spectral radius: {ce_spec['spectral_radius']:.4f}")
print(f"  AAT lambda_1: {ce_spec['aat_lambda1']:.4f}")
print(f"  Odd cycles: {ce_spec['num_cycles']}")
print(f"  Length dist: {ce_spec['len_dist']}")
print(f"  Omega density: {ce_spec['omega_density']:.4f}")
print(f"  Alpha: {ce_spec['alpha']}")

# Score variance as measure of irregularity
mean_score = sum(ce_scores) / n
score_var = sum((s - mean_score)**2 for s in ce_scores) / n
print(f"  Score variance: {score_var:.4f}")

# Part 2: Walk from counterexample toward regularity
print("\n--- Part 2: Arc-flip walk toward regularity ---")
# At each step, flip the arc that most reduces score variance
T_current = [row[:] for row in T_ce]
path = []

for step in range(20):
    coeffs, fail, _, scores, spec = analyze_tournament(T_current, n)
    mean_s = sum(scores) / n
    var = sum((s - mean_s)**2 for s in scores) / n

    status = "FAIL" if fail else "ok"
    cyc_str = f"{spec['num_cycles']}cyc"
    print(f"  Step {step:2d}: var={var:.3f}, {cyc_str}, coeffs={coeffs}, {status}")

    path.append((var, fail, coeffs, spec['num_cycles']))

    if not fail and step > 0:
        print(f"  -> Restored real-rootedness at step {step}!")
        break

    # Find the arc flip that most reduces variance
    best_flip = None
    best_var = var
    for u in range(n):
        for v in range(u+1, n):
            T_try = [row[:] for row in T_current]
            T_try[u][v], T_try[v][u] = T_try[v][u], T_try[u][v]
            new_scores = tuple(sorted(sum(row) for row in T_try))
            new_mean = sum(new_scores) / n
            new_var = sum((s - new_mean)**2 for s in new_scores) / n
            if new_var < best_var:
                best_var = new_var
                best_flip = (u, v)

    if best_flip is None:
        print(f"  -> No variance-reducing flip found!")
        break

    u, v = best_flip
    T_current[u][v], T_current[v][u] = T_current[v][u], T_current[u][v]

# Part 3: Random walk from counterexample
print("\n--- Part 3: Random arc-flip walk (10 trials) ---")
import random
random.seed(123)

for trial in range(10):
    T_walk = [row[:] for row in T_ce]
    steps_to_fix = -1

    for step in range(50):
        # Random arc flip
        u, v = random.choice([(i,j) for i in range(n) for j in range(i+1,n)])
        T_walk[u][v], T_walk[v][u] = T_walk[v][u], T_walk[u][v]

        cycles = find_all_directed_odd_cycles(T_walk, n)
        m = len(cycles)
        if m == 0:
            continue
        cvsets = [frozenset(c) for c in cycles]
        adj = {}
        for i in range(m):
            adj[i] = frozenset(j for j in range(m) if j != i and cvsets[i] & cvsets[j])
        coeffs = independence_poly_dc(adj)
        deg = len(coeffs) - 1
        while deg > 0 and coeffs[deg] == 0:
            deg -= 1

        newton_fail = False
        if deg >= 2:
            for k in range(1, deg):
                if coeffs[k]**2 < coeffs[k-1] * coeffs[k+1] * (k+1) / k:
                    newton_fail = True
                    break

        if not newton_fail:
            steps_to_fix = step + 1
            break

    print(f"  Trial {trial}: fixed after {steps_to_fix} random flips" if steps_to_fix > 0 else f"  Trial {trial}: still failing after 50 flips")

# Part 4: What makes a2 so small in the counterexample?
print("\n--- Part 4: Why is alpha_2 so small? ---")
cycles_ce = find_all_directed_odd_cycles(T_ce, n)
cvsets_ce = [frozenset(c) for c in cycles_ce]
len_dist = Counter(len(c) for c in cycles_ce)
print(f"  {len(cycles_ce)} cycles: {dict(len_dist)}")

# Count independent pairs by length combination
pair_types = Counter()
for i in range(len(cycles_ce)):
    for j in range(i+1, len(cycles_ce)):
        if not (cvsets_ce[i] & cvsets_ce[j]):
            key = tuple(sorted([len(cycles_ce[i]), len(cycles_ce[j])]))
            pair_types[key] += 1

print(f"  Independent pairs by type:")
for key, count in sorted(pair_types.items()):
    print(f"    ({key[0]},{key[1]})-pairs: {count}")

# Which vertices appear in the MOST cycles?
vertex_freq = Counter()
for c in cycles_ce:
    for v in frozenset(c):
        vertex_freq[v] += 1

print(f"\n  Vertex frequency in odd cycles:")
for v in range(n):
    print(f"    v={v} (out-deg={sum(T_ce[v])}): appears in {vertex_freq[v]} cycles ({vertex_freq[v]/len(cycles_ce)*100:.1f}%)")

# The "hub" effect: if one vertex appears in almost all cycles,
# then most pairs share that vertex and can't be independent.
max_v = max(vertex_freq, key=vertex_freq.get)
print(f"\n  Hub vertex: {max_v} (in {vertex_freq[max_v]}/{len(cycles_ce)} = {vertex_freq[max_v]/len(cycles_ce)*100:.1f}% of cycles)")
print(f"  Cycles NOT containing hub: {len(cycles_ce) - vertex_freq[max_v]}")
# Independent pairs can only come from cycles both avoiding the hub,
# or one avoiding hub and the other...
cycles_without_hub = [i for i in range(len(cycles_ce)) if max_v not in cvsets_ce[i]]
print(f"  Cycles not containing vertex {max_v}: {len(cycles_without_hub)}")
for i in cycles_without_hub:
    print(f"    {cycles_ce[i]} (len={len(cycles_ce[i])})")

# Part 5: Score sequence and cycle distribution
print("\n--- Part 5: High-variance score → hub structure? ---")
print("  Hypothesis: Very irregular score sequences create 'hub' vertices")
print("  that appear in almost all odd cycles, suppressing independent pairs.")
print()
# Test: for random n=9 tournaments, compute hub fraction and check Newton
np.random.seed(42)
hub_data = []
for trial in range(100):
    T2 = random_tournament(n)
    scores = sorted(sum(row) for row in T2)
    var2 = sum((s - 4)**2 for s in scores) / n

    cycles2 = find_all_directed_odd_cycles(T2, n)
    m2 = len(cycles2)
    if m2 < 5:
        continue

    vf2 = Counter()
    for c in cycles2:
        for v in frozenset(c):
            vf2[v] += 1
    max_frac = max(vf2.values()) / m2

    cvsets2 = [frozenset(c) for c in cycles2]
    adj2 = {}
    for i in range(m2):
        adj2[i] = frozenset(j for j in range(m2) if j != i and cvsets2[i] & cvsets2[j])
    coeffs2 = independence_poly_dc(adj2)
    deg2 = len(coeffs2) - 1
    while deg2 > 0 and coeffs2[deg2] == 0:
        deg2 -= 1

    newton_fail2 = False
    if deg2 >= 2:
        for k in range(1, deg2):
            if coeffs2[k]**2 < coeffs2[k-1] * coeffs2[k+1] * (k+1) / k:
                newton_fail2 = True
                break

    hub_data.append((var2, max_frac, m2, newton_fail2))

# Sort by hub fraction
hub_data.sort(key=lambda x: -x[1])
print(f"  Tested {len(hub_data)} tournaments. Top 15 by hub fraction:")
print(f"  {'var':>6s} {'hub%':>6s} {'#cyc':>5s} {'Newton':>7s}")
for var2, frac, m2, fail in hub_data[:15]:
    print(f"  {var2:6.3f} {frac*100:5.1f}% {m2:5d} {'FAIL' if fail else 'ok':>7s}")

# Any failures?
failures = [d for d in hub_data if d[3]]
print(f"\n  Newton failures among 100 random: {len(failures)}")

# Counterexample hub fraction
print(f"\n  Counterexample: var={score_var:.3f}, hub={vertex_freq[max_v]/len(cycles_ce)*100:.1f}%, #cyc={len(cycles_ce)}")

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
