"""
ocf_independence_crossover.py — Why the OCF independence structure creates crossover

Key observation: the 2^k-weighted cycle count Σ C_k * 2^k favors interval
even at p=7 where Paley actually wins H. This means the OCF's INDEPENDENCE
POLYNOMIAL (not raw cycle count) is the crucial quantity.

H(T) = I(Ω(T), 2) = Σ_{k=0}^∞ α_k * 2^k

where α_k = number of independent sets of size k in the odd-cycle graph Ω(T).
Two cycles are "independent" iff they are vertex-disjoint.

This script computes:
1. The full odd-cycle complex for Paley and Interval at p=7, 11
2. The independence polynomial I(Ω, x)
3. The breakdown of α_k by cycle length composition
4. Which α_k terms give Paley vs Interval the advantage

Author: opus-2026-03-12-S60
"""
import sys
import time
import numpy as np
from collections import defaultdict, Counter
from itertools import combinations
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)


def circulant_adj(n, S):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for s in S:
            A[i][(i+s)%n] = 1
    return A


def paley_set(p):
    return frozenset(pow(x, 2, p) for x in range(1, p))


def interval_set(p):
    m = (p - 1) // 2
    return frozenset(range(1, m + 1))


def find_all_directed_cycles(A, n, max_len=None):
    """Find all directed cycles using DFS. Returns list of (length, vertex_set)."""
    if max_len is None:
        max_len = n
    cycles = []

    for start in range(n):
        # DFS from start, only visiting vertices > start (to avoid double-counting)
        stack = [(start, [start], {start})]
        while stack:
            v, path, visited = stack.pop()
            for w in range(n):
                if not A[v][w]:
                    continue
                if w == start and len(path) >= 3:
                    # Found a cycle
                    cycles.append((len(path), frozenset(path)))
                elif w > start and w not in visited and len(path) < max_len:
                    stack.append((w, path + [w], visited | {w}))

    return cycles


def find_odd_cycles(A, n):
    """Find all odd directed cycles."""
    all_cycles = find_all_directed_cycles(A, n)
    return [(length, vset) for length, vset in all_cycles if length % 2 == 1]


def independence_polynomial(cycles, n):
    """Compute the independence polynomial of the cycle graph.

    Two cycles are adjacent (dependent) if they share a vertex.
    α_k = number of independent sets of size k.

    Returns dict {k: α_k}.
    """
    m = len(cycles)
    if m == 0:
        return {0: 1}

    # Build adjacency: cycles[i] and cycles[j] are dependent iff they share a vertex
    adj = [set() for _ in range(m)]
    for i in range(m):
        for j in range(i + 1, m):
            if cycles[i][1] & cycles[j][1]:  # shared vertices
                adj[i].add(j)
                adj[j].add(i)

    # Compute α_k by enumeration for small m
    alpha = defaultdict(int)
    alpha[0] = 1

    if m <= 25:
        # Direct enumeration: try all subsets
        for size in range(1, m + 1):
            for subset in combinations(range(m), size):
                # Check independence
                ok = True
                for i in range(len(subset)):
                    for j in range(i + 1, len(subset)):
                        if subset[j] in adj[subset[i]]:
                            ok = False
                            break
                    if not ok:
                        break
                if ok:
                    alpha[size] += 1
    else:
        # For larger m, use greedy approximation for α_1 and α_2
        alpha[1] = m  # Every single cycle is an independent set of size 1

        # α_2 = number of pairs of vertex-disjoint cycles
        count2 = 0
        for i in range(m):
            for j in range(i + 1, m):
                if j not in adj[i]:
                    count2 += 1
        alpha[2] = count2

        # α_3: triple of mutually vertex-disjoint cycles
        count3 = 0
        for i in range(m):
            for j in range(i + 1, m):
                if j in adj[i]:
                    continue
                for k in range(j + 1, m):
                    if k not in adj[i] and k not in adj[j]:
                        count3 += 1
        alpha[3] = count3

    return dict(alpha)


def independence_polynomial_by_composition(cycles, n):
    """Like independence_polynomial but track which cycle lengths are in each IS."""
    m = len(cycles)
    if m == 0:
        return {}

    adj = [set() for _ in range(m)]
    for i in range(m):
        for j in range(i + 1, m):
            if cycles[i][1] & cycles[j][1]:
                adj[i].add(j)
                adj[j].add(i)

    # For each independent set, record the sorted tuple of cycle lengths
    composition_counts = defaultdict(int)

    if m <= 22:
        for size in range(1, min(m + 1, 5)):  # up to size 4
            for subset in combinations(range(m), size):
                ok = True
                for i in range(len(subset)):
                    for j in range(i + 1, len(subset)):
                        if subset[j] in adj[subset[i]]:
                            ok = False
                            break
                    if not ok:
                        break
                if ok:
                    lengths = tuple(sorted(cycles[s][0] for s in subset))
                    composition_counts[lengths] += 1

    return dict(composition_counts)


def hamiltonian_paths_dp(A, n):
    dp = defaultdict(lambda: defaultdict(int))
    for v in range(n):
        dp[1 << v][v] = 1
    full = (1 << n) - 1
    for mask in range(1, full + 1):
        if not dp[mask]:
            continue
        for v in dp[mask]:
            if dp[mask][v] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    dp[mask | (1 << w)][w] += dp[mask][v]
    return sum(dp[full][v] for v in range(n))


def main():
    print("OCF INDEPENDENCE STRUCTURE: PALEY vs INTERVAL")
    print("=" * 75)

    for p in [7, 11]:
        print(f"\n{'='*75}")
        print(f"p = {p}")
        print(f"{'='*75}")

        m = (p - 1) // 2
        S_paley = paley_set(p)
        S_interval = interval_set(p)

        for name, S in [("Paley", S_paley), ("Interval", S_interval)]:
            print(f"\n  {name} (S = {sorted(S)}):")
            A = circulant_adj(p, S)
            H = hamiltonian_paths_dp(A, p)

            # Find all odd cycles
            t0 = time.time()
            odd_cycles = find_odd_cycles(A, p)
            elapsed = time.time() - t0
            print(f"    {len(odd_cycles)} odd cycles found ({elapsed:.2f}s)")

            # Cycle length distribution
            len_dist = Counter(c[0] for c in odd_cycles)
            for length in sorted(len_dist.keys()):
                print(f"      C_{length} = {len_dist[length]}")

            # Independence polynomial
            t0 = time.time()
            alpha = independence_polynomial(odd_cycles, p)
            elapsed = time.time() - t0
            print(f"    Independence polynomial ({elapsed:.2f}s):")
            H_from_ocf = sum(v * (2 ** k) for k, v in alpha.items())
            for k in sorted(alpha.keys()):
                contribution = alpha[k] * (2 ** k)
                print(f"      α_{k} = {alpha[k]:>10}, 2^{k}·α_{k} = {contribution:>10}")
            print(f"    I(Ω, 2) = {H_from_ocf}")
            print(f"    H (DP)  = {H}")
            print(f"    Match: {H_from_ocf == H}")

            # Composition breakdown (which cycle lengths contribute to α_k)
            if len(odd_cycles) <= 22:
                print(f"    Composition breakdown (cycle lengths in each IS):")
                comp = independence_polynomial_by_composition(odd_cycles, p)
                for lengths in sorted(comp.keys()):
                    k = len(lengths)
                    count = comp[lengths]
                    contribution = count * (2 ** k)
                    print(f"      {lengths}: {count:>6} ISets, 2^{k} contribution = {contribution:>8}")

    # Part 2: The crucial difference
    print(f"\n{'='*75}")
    print("PART 2: WHERE DOES PALEY GAIN ITS ADVANTAGE?")
    print("-" * 75)

    for p in [7, 11]:
        print(f"\n  p = {p}:")
        m = (p - 1) // 2

        # Compute both
        results = {}
        for name, S in [("Paley", paley_set(p)), ("Interval", interval_set(p))]:
            A = circulant_adj(p, S)
            odd_cycles = find_odd_cycles(A, p)
            alpha = independence_polynomial(odd_cycles, p)
            results[name] = alpha

        # Compare α_k values
        all_k = sorted(set(results["Paley"].keys()) | set(results["Interval"].keys()))
        print(f"  {'k':>4} {'α_k(P)':>10} {'α_k(I)':>10} {'Δα_k':>10} {'2^k·Δα_k':>12}")
        total_diff = 0
        for k in all_k:
            aP = results["Paley"].get(k, 0)
            aI = results["Interval"].get(k, 0)
            delta = aP - aI
            weighted = delta * (2 ** k)
            total_diff += weighted
            print(f"  {k:>4} {aP:>10} {aI:>10} {delta:>+10} {weighted:>+12}")
        print(f"  Total Δ = {total_diff} = H(Paley) - H(Interval)")

    # Part 3: At p=11, which cycle length compositions drive the α difference?
    print(f"\n{'='*75}")
    print("PART 3: COMPOSITION-LEVEL ANALYSIS AT p=11")
    print("-" * 75)

    p = 11
    comp_results = {}
    for name, S in [("Paley", paley_set(p)), ("Interval", interval_set(p))]:
        A = circulant_adj(p, S)
        odd_cycles = find_odd_cycles(A, p)
        comp = independence_polynomial_by_composition(odd_cycles, p)
        comp_results[name] = comp

    all_comps = sorted(set(comp_results["Paley"].keys()) | set(comp_results["Interval"].keys()))
    print(f"  {'Composition':>20} {'Paley':>8} {'Interval':>8} {'Δ':>8} {'2^k·Δ':>10}")
    for lengths in all_comps:
        cP = comp_results["Paley"].get(lengths, 0)
        cI = comp_results["Interval"].get(lengths, 0)
        delta = cP - cI
        k = len(lengths)
        weighted = delta * (2 ** k)
        if abs(delta) > 0:
            print(f"  {str(lengths):>20} {cP:>8} {cI:>8} {delta:>+8} {weighted:>+10}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
