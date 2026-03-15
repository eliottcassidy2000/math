#!/usr/bin/env python3
"""
simplicial_H_n6_s112.py — Simplicial H at n=6 + characterization of sim_H=1
kind-pasteur-2026-03-15-S112

KEY FINDINGS SO FAR:
  sim_H in {0, 1} at n=4,5 (BINARY invariant)
  sim_H = 1 iff transitive core is a total order (DAG with unique LE)
  sim_H = 0 iff transitive core has cycles (contradictory constraints)

THIS SCRIPT:
  1. Check n=6 (2^15 = 32768 tournaments)
  2. Characterize sim_H=1 tournaments (score sequences, structure)
  3. What fraction have sim_H=1 as n grows?
  4. Is there a formula for P(sim_H=1)?
"""

from itertools import permutations, combinations
from collections import Counter, defaultdict
from math import factorial, comb

def analyze_fast(T, n):
    """Fast analysis: compute H, sim_H, score sequence."""
    # Build transitive triple constraints
    constraints = set()  # set of (a,b) meaning a must come before b
    n_3cycles = 0

    for i, j, k in combinations(range(n), 3):
        # Check if transitive
        s = T[i][j] + T[j][k] + T[k][i]
        if s == 3 or s == 0:
            n_3cycles += 1
            continue
        # Find the transitive ordering
        for a, b, c in [(i,j,k),(i,k,j),(j,i,k),(j,k,i),(k,i,j),(k,j,i)]:
            if T[a][b] and T[b][c] and T[a][c]:
                constraints.add((a, b))
                constraints.add((b, c))
                constraints.add((a, c))
                break

    # Check if constraints form a DAG (no cycles)
    # Build adjacency list
    adj = defaultdict(set)
    for a, b in constraints:
        adj[a].add(b)

    # Topological sort attempt
    in_degree = defaultdict(int)
    for a, b in constraints:
        in_degree[b] += 1
        if a not in in_degree:
            in_degree[a] += 0  # ensure all nodes present

    all_nodes = set()
    for a, b in constraints:
        all_nodes.add(a)
        all_nodes.add(b)

    # Kahn's algorithm
    queue = [v for v in range(n) if in_degree.get(v, 0) == 0]
    topo_order = []
    remaining = dict(in_degree)

    while queue:
        if len(queue) > 1:
            # Multiple choices: not a unique total order -> sim_H > 1 (or 0 if cycle)
            # Actually if DAG but not chain: multiple LEs
            pass
        v = queue.pop(0)
        topo_order.append(v)
        for u in adj.get(v, []):
            remaining[u] -= 1
            if remaining[u] == 0:
                queue.append(u)

    is_dag = len(topo_order) == len(all_nodes) or len(all_nodes) == 0

    # If DAG: count linear extensions
    # For efficiency: check if it's a CHAIN (total order)
    is_chain = is_dag and len(constraints) == n*(n-1)//2
    # Actually: a chain on n vertices has C(n,2) edges.
    # But constraints might not cover all vertices (isolated vertices from 3-cycles)

    # More precisely: constraints form a total order iff for every pair (i,j),
    # either (i,j) or (j,i) is in constraints.
    is_total = True
    for i in range(n):
        for j in range(i+1, n):
            if (i,j) not in constraints and (j,i) not in constraints:
                is_total = False
                break
        if not is_total:
            break

    if not is_dag:
        sim_H = 0
    elif is_total:
        sim_H = 1
    else:
        # DAG but not total: count linear extensions (could be > 1)
        # For small n, brute force
        sim_H = 0
        for perm in permutations(range(n)):
            pos = {perm[i]: i for i in range(n)}
            ok = True
            for a, b in constraints:
                if pos[a] > pos[b]:
                    ok = False
                    break
            if ok:
                sim_H += 1

    # Count H (tournament Hamiltonian paths)
    H = 0
    for perm in permutations(range(n)):
        ok = True
        for p in range(n-1):
            if T[perm[p]][perm[p+1]] != 1:
                ok = False
                break
        if ok:
            H += 1

    # Score sequence
    scores = tuple(sorted([sum(T[i][j] for j in range(n)) for i in range(n)]))

    return H, sim_H, n_3cycles, scores

# ============================================================
# n=6 analysis
# ============================================================
print("="*70)
print("n=6: EXHAUSTIVE ANALYSIS (32768 tournaments)")
print("="*70)

n = 6
total = 2**(n*(n-1)//2)
sH_dist = Counter()
H_sH_dist = Counter()
score_sH = defaultdict(lambda: Counter())
sim_H_values = set()

count = 0
for bits in range(total):
    T = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                T[i][j] = 1
            else:
                T[j][i] = 1
            idx += 1

    H, sH, c3, scores = analyze_fast(T, n)
    sH_dist[sH] += 1
    H_sH_dist[(H, sH)] += 1
    score_sH[scores][sH] += 1
    sim_H_values.add(sH)

    count += 1
    if count % 5000 == 0:
        print(f"  {count}/{total} processed...")

print(f"\nDone! {total} tournaments analyzed.")

print(f"\nsim_H values at n=6: {sorted(sim_H_values)}")
print(f"\nsim_H distribution:")
for sH in sorted(sH_dist.keys()):
    pct = sH_dist[sH] / total * 100
    print(f"  sim_H={sH}: {sH_dist[sH]} ({pct:.1f}%)")

# Is sim_H still in {0, 1}?
binary = all(sH in [0, 1] for sH in sim_H_values)
print(f"\nsim_H is BINARY (0 or 1): {binary}")

print(f"\n(H, sim_H) distribution:")
for key in sorted(H_sH_dist.keys()):
    H, sH = key
    print(f"  (H={H:2d}, sH={sH}): {H_sH_dist[key]:5d}")

# Score sequences with sim_H=1
print(f"\nScore sequences with sim_H=1:")
for scores in sorted(score_sH.keys()):
    if 1 in score_sH[scores]:
        total_with_score = sum(score_sH[scores].values())
        n1 = score_sH[scores][1]
        n0 = score_sH[scores].get(0, 0)
        print(f"  scores={scores}: {n1}/{total_with_score} have sim_H=1 ({n1/total_with_score*100:.0f}%)")

# Fraction with sim_H = 1 as function of n
print(f"\nFraction with sim_H=1:")
print(f"  n=4: 48/64 = {48/64:.4f}")
print(f"  n=5: 240/1024 = {240/1024:.4f}")
n1_n6 = sH_dist.get(1, 0)
print(f"  n=6: {n1_n6}/{total} = {n1_n6/total:.4f}")

# Conjecture: does P(sim_H=1) -> 0 as n -> infinity?
# Or does it stabilize?
ratios = [48/64, 240/1024, n1_n6/total]
print(f"\n  Ratios: {[f'{r:.4f}' for r in ratios]}")
print(f"  Declining? {ratios[0] > ratios[1] > ratios[2]}")
if ratios[1] > 0 and ratios[2] > 0:
    log_ratio = [ratios[i+1]/ratios[i] for i in range(len(ratios)-1)]
    print(f"  Successive ratios: {[f'{r:.4f}' for r in log_ratio]}")

print("\nDone!")
