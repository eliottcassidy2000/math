"""
alpha4_forcing_analysis.py -- kind-pasteur-2026-03-14-S66

Alpha_1=4 forcing theorem: at n=8, alpha_2 is in {0,4} ONLY.
alpha_2=1,2,3 never occur. This is a binary phase transition.

Check this at n=7, n=9, n=10 to see if universal.
Also: analyze the disjointness graph topology for alpha_2=4.
"""

import numpy as np
from itertools import combinations
from collections import Counter

def count_directed_hamcycles(A, vertices):
    k = len(vertices)
    if k < 3 or k % 2 == 0:
        return 0
    vlist = list(vertices)
    sub = [[0]*k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            if i != j:
                sub[i][j] = int(A[vlist[i]][vlist[j]])
    full = (1 << k) - 1
    dp = [[0]*k for _ in range(1 << k)]
    dp[1][0] = 1
    for mask in range(1, 1 << k):
        for v in range(k):
            if dp[mask][v] == 0:
                continue
            for u in range(1, k):
                if mask & (1 << u):
                    continue
                if sub[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    total = 0
    for v in range(1, k):
        if dp[full][v] and sub[v][0]:
            total += dp[full][v]
    return total

def find_cycles(A, n, max_size=None):
    """Find all directed odd cycles in tournament A on n vertices."""
    if max_size is None:
        max_size = n
    cycles = []
    for size in range(3, max_size+1, 2):
        for subset in combinations(range(n), size):
            cnt = count_directed_hamcycles(A, list(subset))
            if cnt > 0:
                cycles.append((frozenset(subset), cnt, size))
    return cycles

def disjointness_graph(cycles):
    """Return list of disjoint pairs (i,j) with i<j."""
    pairs = []
    for i in range(len(cycles)):
        for j in range(i+1, len(cycles)):
            if len(cycles[i][0] & cycles[j][0]) == 0:
                pairs.append((i,j))
    return pairs

def classify_graph_4v(edges):
    """Classify a graph on 4 vertices by its edge count and structure."""
    adj = {i: set() for i in range(4)}
    for i, j in edges:
        adj[i].add(j)
        adj[j].add(i)
    degrees = sorted([len(adj[i]) for i in range(4)])
    ne = len(edges)
    if ne == 0:
        return "empty"
    elif ne == 1:
        return "K2"
    elif ne == 2:
        if degrees == [0, 0, 1, 1]:
            return "2K2(matching)"
        elif degrees == [0, 1, 1, 2]:
            return "P3"
        else:
            return f"2e-{degrees}"
    elif ne == 3:
        if degrees == [0, 2, 2, 2]:
            return "K3+K1"
        elif degrees == [1, 1, 1, 3]:
            return "K1,3(star)"
        elif degrees == [1, 1, 2, 2]:
            return "P4(path)"
        else:
            return f"3e-{degrees}"
    elif ne == 4:
        if degrees == [1, 2, 2, 3]:
            return "paw"
        elif degrees == [2, 2, 2, 2]:
            return "C4"
        else:
            return f"4e-{degrees}"
    elif ne == 5:
        return "K4-e"
    elif ne == 6:
        return "K4"
    else:
        return f"{ne}e"

def main():
    print("=" * 70)
    print("ALPHA_1=4 FORCING ANALYSIS: IS alpha_2 IN {0,4} UNIVERSAL?")
    print("=" * 70)

    for n in [7, 8, 9]:
        print(f"\n{'='*50}")
        print(f"n = {n}")
        print(f"{'='*50}")

        rng = np.random.default_rng(2026_03_14_00 + n)
        a2_values = []
        topology_counts = Counter()

        num_samples = {7: 500000, 8: 500000, 9: 200000}[n]

        for trial in range(num_samples):
            A = np.zeros((n, n), dtype=np.int8)
            for i in range(n):
                for j in range(i+1, n):
                    if rng.random() < 0.5:
                        A[i][j] = 1
                    else:
                        A[j][i] = 1

            # Quick c3 filter
            c3 = 0
            for a, b, c in combinations(range(n), 3):
                if (A[a][b] and A[b][c] and A[c][a]) or \
                   (A[a][c] and A[c][b] and A[b][a]):
                    c3 += 1
            if c3 > 8:  # Generous filter for alpha_1=4
                continue

            cycles = find_cycles(A, n)
            alpha_1 = sum(cnt for _, cnt, _ in cycles)

            if alpha_1 == 4:
                # Compute alpha_2
                disj_pairs = disjointness_graph(cycles)
                alpha_2 = sum(cycles[i][1] * cycles[j][1] for i, j in disj_pairs)
                a2_values.append(alpha_2)

                # Classify topology if exactly 4 cycle vertex-sets
                # (could be fewer sets if some have multiple directed cycles)
                cycle_sets = [(vs, cnt, sz) for vs, cnt, sz in cycles]

                if len(cycle_sets) <= 4:
                    edges = [(i,j) for i in range(len(cycle_sets))
                             for j in range(i+1, len(cycle_sets))
                             if len(cycle_sets[i][0] & cycle_sets[j][0]) == 0]
                    if len(cycle_sets) == 4:
                        topo = classify_graph_4v(edges)
                    else:
                        topo = f"{len(cycle_sets)}sets_{len(edges)}disj"
                    topology_counts[topo] += 1

            if (trial+1) % 100000 == 0:
                print(f"  {trial+1}/{num_samples}, alpha_1=4 found: {len(a2_values)}")

        print(f"\nalpha_1=4 at n={n}: {len(a2_values)} found")
        if a2_values:
            a2_counts = Counter(a2_values)
            print("alpha_2 distribution:")
            for a2 in sorted(a2_counts.keys()):
                pct = 100 * a2_counts[a2] / len(a2_values)
                print(f"  alpha_2={a2}: {a2_counts[a2]} ({pct:.1f}%)")
            print(f"\nDisjointness graph topologies:")
            for topo in sorted(topology_counts.keys()):
                print(f"  {topo}: {topology_counts[topo]}")
        else:
            print("  No alpha_1=4 found!")

    # Deep dive at n=8: what EXACTLY are the alpha_2=4 tournaments?
    print(f"\n{'='*70}")
    print("DEEP DIVE: alpha_2=4 at n=8")
    print("=" * 70)

    n = 8
    rng3 = np.random.default_rng(2026_03_14_88)
    examples_a2_4 = []
    examples_a2_0 = []

    for trial in range(500000):
        A = np.zeros((n, n), dtype=np.int8)
        for i in range(n):
            for j in range(i+1, n):
                if rng3.random() < 0.5:
                    A[i][j] = 1
                else:
                    A[j][i] = 1

        c3 = 0
        for a, b, c in combinations(range(n), 3):
            if (A[a][b] and A[b][c] and A[c][a]) or \
               (A[a][c] and A[c][b] and A[b][a]):
                c3 += 1
        if c3 > 8:
            continue

        cycles = find_cycles(A, n)
        alpha_1 = sum(cnt for _, cnt, _ in cycles)

        if alpha_1 == 4:
            disj_pairs = disjointness_graph(cycles)
            alpha_2 = sum(cycles[i][1] * cycles[j][1] for i, j in disj_pairs)

            if alpha_2 == 4 and len(examples_a2_4) < 5:
                scores = tuple(sorted(int(sum(A[i])) for i in range(n)))
                cycle_info = [(sorted(vs), cnt, sz) for vs, cnt, sz in cycles]
                examples_a2_4.append((scores, cycle_info, len(disj_pairs)))
            elif alpha_2 == 0 and len(examples_a2_0) < 5:
                scores = tuple(sorted(int(sum(A[i])) for i in range(n)))
                cycle_info = [(sorted(vs), cnt, sz) for vs, cnt, sz in cycles]
                examples_a2_0.append((scores, cycle_info, 0))

    print(f"\nalpha_2=0 examples (first 5):")
    for scores, cycles, nd in examples_a2_0:
        print(f"  Score: {scores}")
        for vs, cnt, sz in cycles:
            print(f"    Cycle: {vs} (size {sz}, {cnt} directed)")
        print()

    print(f"alpha_2=4 examples (first 5):")
    for scores, cycles, nd in examples_a2_4:
        print(f"  Score: {scores}, disjoint_pairs={nd}")
        for vs, cnt, sz in cycles:
            print(f"    Cycle: {vs} (size {sz}, {cnt} directed)")
        # Check which pairs are disjoint
        sets = [frozenset(vs) for vs, cnt, sz in cycles]
        for i in range(len(sets)):
            for j in range(i+1, len(sets)):
                d = "DISJOINT" if len(sets[i] & sets[j]) == 0 else f"share {sets[i] & sets[j]}"
                print(f"    {i}-{j}: {d}")
        print()

    # Check: with alpha_1=4, can we THEORETICALLY have alpha_2=3?
    print("=" * 70)
    print("THEORETICAL ANALYSIS: WHY alpha_2 IN {0,4}?")
    print("=" * 70)
    print()
    print("alpha_1=4 means exactly 4 directed odd cycles (counting multiplicity).")
    print("These live on some set of vertex-SETS. Multiple directed cycles can")
    print("live on the same vertex set (e.g., a 5-vertex set with 3 Hamiltonian cycles).")
    print()
    print("Case analysis:")
    print("  If 4 cycles on 4 distinct vertex sets:")
    print("    alpha_2 = sum of products over disjoint pairs")
    print("    Each cycle has multiplicity 1. alpha_2 = #disjoint_pairs.")
    print("    So alpha_2 in {0,1,2,3,4,5,6} a priori.")
    print("    But data shows only {0,4}.")
    print()
    print("  If 4 cycles on fewer vertex sets (e.g., 2 cycles on one 5-vertex set):")
    print("    alpha_2 = sum c_i*c_j for disjoint set-pairs")
    print("    Products can be >1.")

if __name__ == "__main__":
    main()
