#!/usr/bin/env python3
"""
digraph_metagraph_s20dm.py — The metagraph hierarchy: tournaments -> oriented -> digraphs
kind-pasteur-2026-03-23-S20dm

HIERARCHY:
  Level 0: TOURNAMENT — 2 states per tile (fwd/bwd), m = C(n,2) arcs
  Level 1: ORIENTED GRAPH — 3 states per tile (fwd/bwd/missing), m = C(n,2) positions
  Level 2: DIGRAPH — 4 states per pair (missing/fwd/bwd/both), m = C(n,2) positions

At each level, the metagraph is the quotient of the state space by S_n:
  Level 0: G_n = {0,1}^m / S_n (binary hypercube quotient)
  Level 1: O_n = {0,1,2}^m / S_n (ternary hypercube quotient)
  Level 2: D_n = {0,1,2,3}^m / S_n (quaternary hypercube quotient)

Compute the metagraph at each level for small n, compare structure.
"""

import sys
from math import comb, factorial
from itertools import permutations, product
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  METAGRAPH HIERARCHY: TOURNAMENT -> ORIENTED -> DIGRAPH")
print("  kind-pasteur-2026-03-23-S20dm")
print("=" * 80)

def analyze_level(n, num_states, level_name):
    """Analyze the metagraph at a given level."""
    t0 = time.time()

    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(pairs)
    total = num_states ** m

    if total > 200000:
        print(f"\n  {level_name} at n={n}: {total} states = {num_states}^{m} — TOO LARGE, skipping")
        return None

    perms = list(permutations(range(n)))

    def state_to_adj(state_tuple):
        """Convert state tuple to adjacency."""
        # State 0: no edge (missing)
        # State 1: i->j (forward)
        # State 2: j->i (backward)
        # State 3: both i->j and j->i (bidirectional) — only for digraphs
        A = [[0]*n for _ in range(n)]
        for k, (i,j) in enumerate(pairs):
            s = state_tuple[k]
            if num_states == 2:
                # Tournament: 0=fwd, 1=bwd
                if s == 0: A[i][j] = 1
                else: A[j][i] = 1
            elif num_states == 3:
                # Oriented: 0=missing, 1=fwd, 2=bwd
                if s == 1: A[i][j] = 1
                elif s == 2: A[j][i] = 1
            elif num_states == 4:
                # Digraph: 0=missing, 1=fwd, 2=bwd, 3=both
                if s == 1: A[i][j] = 1
                elif s == 2: A[j][i] = 1
                elif s == 3: A[i][j] = 1; A[j][i] = 1
        return A

    def canon(A):
        best = None
        for p in perms:
            s = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
            if best is None or s < best: best = s
        return best

    def Hdp(A):
        """Count directed Hamiltonian paths."""
        dp = {}
        for v in range(n): dp[(1<<v, v)] = 1
        for S in range(1, 1<<n):
            for v in range(n):
                if not (S & (1<<v)): continue
                val = dp.get((S,v), 0)
                if val == 0: continue
                for u in range(n):
                    if S & (1<<u): continue
                    if A[v][u]:
                        dp[(S|(1<<u), u)] = dp.get((S|(1<<u), u), 0) + val
        return sum(dp.get(((1<<n)-1, v), 0) for v in range(n))

    # Enumerate all states
    classes = {}
    for idx, state in enumerate(product(range(num_states), repeat=m)):
        A = state_to_adj(state)
        cn = canon(A)
        if cn not in classes:
            H = Hdp(A)
            edge_count = sum(sum(row) for row in A)
            classes[cn] = {'H': H, 'edges': edge_count, 'count': 0}
        classes[cn]['count'] += 1

    # Compute metagraph edges (single-position flip to any other state)
    # An edge exists between two iso classes if some single-position state change connects them
    V = len(classes)
    class_list = sorted(classes.keys())
    class_idx = {cn: i for i, cn in enumerate(class_list)}

    # For small state spaces, compute edges
    edges = set()
    if total <= 50000:
        state_to_class = {}
        for state in product(range(num_states), repeat=m):
            A = state_to_adj(state)
            cn = canon(A)
            state_to_class[state] = cn

        for state in product(range(num_states), repeat=m):
            cn1 = state_to_class[state]
            ci1 = class_idx[cn1]
            for k in range(m):
                for new_val in range(num_states):
                    if new_val == state[k]: continue
                    state2 = list(state)
                    state2[k] = new_val
                    state2 = tuple(state2)
                    cn2 = state_to_class[state2]
                    ci2 = class_idx[cn2]
                    if ci1 != ci2:
                        edges.add((min(ci1,ci2), max(ci1,ci2)))

    E = len(edges)

    # H distribution
    H_dist = Counter(classes[cn]['H'] for cn in class_list)
    edge_dist = Counter(classes[cn]['edges'] for cn in class_list)

    print(f"\n  {level_name} at n={n}:")
    print(f"    States per position: {num_states}")
    print(f"    Total states: {total} = {num_states}^{m}")
    print(f"    Iso classes (V): {V}")
    print(f"    Edges (E): {E}")
    print(f"    H values: min={min(H_dist)}, max={max(H_dist)}, distinct={len(H_dist)}")
    print(f"    H=0 classes (no Ham path): {H_dist.get(0, 0)}")
    print(f"    H odd classes: {sum(1 for cn in class_list if classes[cn]['H'] % 2 == 1)}")
    print(f"    H even classes: {sum(1 for cn in class_list if classes[cn]['H'] % 2 == 0)}")
    print(f"    Time: {time.time()-t0:.1f}s")

    # Edge density distribution
    print(f"    Edge density: {dict(sorted(edge_dist.items())[:10])}...")

    return {'V': V, 'E': E, 'H_dist': H_dist, 'level': level_name, 'n': n, 'states': num_states}


# ============================================================================
# MAIN
# ============================================================================

all_results = {}

for n in [3, 4, 5]:
    print(f"\n{'#'*60}")
    print(f"  n = {n}")
    print(f"{'#'*60}")

    # Tournament (2 states)
    r = analyze_level(n, 2, "TOURNAMENT")
    if r: all_results[(n, 2)] = r

    # Oriented graph (3 states)
    r = analyze_level(n, 3, "ORIENTED")
    if r: all_results[(n, 3)] = r

    # Digraph (4 states) — only for small n
    r = analyze_level(n, 4, "DIGRAPH")
    if r: all_results[(n, 4)] = r


# ============================================================================
# CROSS-LEVEL COMPARISON
# ============================================================================

print(f"\n\n{'='*60}")
print(f"  CROSS-LEVEL COMPARISON")
print(f"{'='*60}")

print(f"\n  {'n':>3} {'level':>12} {'states':>7} {'V':>8} {'E':>8} {'H=0':>5} {'H_odd':>6} {'H_even':>6}")
for key in sorted(all_results.keys()):
    r = all_results[key]
    n = r['n']
    hd = r['H_dist']
    h0 = hd.get(0, 0)
    h_odd = sum(v for k, v in hd.items() if k > 0 and k % 2 == 1)
    h_even = sum(v for k, v in hd.items() if k > 0 and k % 2 == 0)
    print(f"  {n:3d} {r['level']:>12} {r['states']:7d} {r['V']:8d} {r['E']:8d} {h0:5d} {h_odd:6d} {h_even:6d}")

# Key sequence: OEIS lookup candidates
print(f"\n  V(tournament) = A000568: {[all_results.get((n,2),{}).get('V','?') for n in [3,4,5]]}")
print(f"  V(oriented) = A001174: {[all_results.get((n,3),{}).get('V','?') for n in [3,4,5]]}")
print(f"  V(digraph) = A000273: {[all_results.get((n,4),{}).get('V','?') for n in [3,4,5]]}")

print(f"\n  DONE.")
print("=" * 80)
