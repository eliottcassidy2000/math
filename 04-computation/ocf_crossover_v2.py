"""
ocf_crossover_v2.py — OCF independence analysis, optimized

For p=7 (80 cycles), we can do full independence polynomial.
For p=11 (21K+ cycles), we need the DP-verified H values and just
compare the α_k breakdown where feasible.

The key insight from p=7:
  Paley: α_1=80, α_2=7  → H = 1 + 160 + 28 = 189
  Interval: α_1=59, α_2=14 → H = 1 + 118 + 56 = 175

So Paley wins because it has MORE odd cycles (80 vs 59) at level α_1,
even though Interval has MORE independent pairs (14 vs 7) at level α_2.
The 2^1 weighting at α_1 dominates the 2^2 weighting at α_2.

This means: at p=7, Paley's TOTAL cycle count advantage overwhelms
Interval's better independence structure.

At larger p, the independence structure becomes richer (more vertex-disjoint
cycle configurations), so higher α_k terms dominate, eventually favoring
Interval's numerous long cycles.

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


def find_odd_cycles(A, n):
    """Find all odd directed cycles using DFS."""
    cycles = []
    for start in range(n):
        stack = [(start, [start], 1 << start)]
        while stack:
            v, path, mask = stack.pop()
            for w in range(n):
                if not A[v][w]:
                    continue
                if w == start and len(path) >= 3 and len(path) % 2 == 1:
                    cycles.append((len(path), frozenset(path)))
                elif w > start and not (mask & (1 << w)) and len(path) < n:
                    stack.append((w, path + [w], mask | (1 << w)))
    return cycles


def independence_polynomial_small(cycles, n):
    """Full independence polynomial for small cycle sets."""
    m = len(cycles)
    adj = [set() for _ in range(m)]
    for i in range(m):
        for j in range(i + 1, m):
            if cycles[i][1] & cycles[j][1]:
                adj[i].add(j)
                adj[j].add(i)

    alpha = defaultdict(int)
    alpha[0] = 1

    for size in range(1, m + 1):
        count = 0
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
                count += 1
        if count == 0:
            break
        alpha[size] = count
    return dict(alpha)


def alpha_low(cycles, n, max_k=3):
    """Compute α_1, α_2, α_3 efficiently for large cycle sets."""
    m = len(cycles)
    alpha = {0: 1, 1: m}

    if max_k >= 2:
        # α_2 = number of vertex-disjoint cycle pairs
        count2 = 0
        for i in range(m):
            for j in range(i + 1, m):
                if not (cycles[i][1] & cycles[j][1]):
                    count2 += 1
            if i % 1000 == 999:
                print(f"      α_2 progress: {i+1}/{m}...")
        alpha[2] = count2

    if max_k >= 3 and m < 5000:
        # α_3 using adjacency lists
        adj_list = [[] for _ in range(m)]
        non_adj = [set() for _ in range(m)]
        for i in range(m):
            for j in range(i + 1, m):
                if cycles[i][1] & cycles[j][1]:
                    adj_list[i].append(j)
                    adj_list[j].append(i)

        # Convert to non-adjacency for efficiency
        count3 = 0
        for i in range(m):
            adj_i = set(adj_list[i])
            indep_after_i = [j for j in range(i + 1, m) if j not in adj_i]
            for idx_j, j in enumerate(indep_after_i):
                adj_j = set(adj_list[j])
                for k_idx in range(idx_j + 1, len(indep_after_i)):
                    k = indep_after_i[k_idx]
                    if k not in adj_j:
                        count3 += 1
            if i % 500 == 499:
                print(f"      α_3 progress: {i+1}/{m}...")
        alpha[3] = count3

    return alpha


def main():
    print("OCF CROSSOVER v2: INDEPENDENCE STRUCTURE")
    print("=" * 75)

    # Part 1: p=7 full analysis
    print(f"\n{'='*75}")
    print(f"p = 7: FULL ANALYSIS")
    print(f"{'='*75}")

    p = 7
    for name, S in [("Paley", paley_set(p)), ("Interval", interval_set(p))]:
        A = circulant_adj(p, S)
        H = hamiltonian_paths_dp(A, p)
        odd_cycles = find_odd_cycles(A, p)
        alpha = independence_polynomial_small(odd_cycles, p)

        print(f"\n  {name}: S = {sorted(S)}, H = {H}")
        len_dist = Counter(c[0] for c in odd_cycles)
        print(f"    Cycle counts: {dict(sorted(len_dist.items()))}")
        for k in sorted(alpha.keys()):
            print(f"    α_{k} = {alpha[k]:>6}  (2^{k}·α_{k} = {alpha[k] * 2**k:>6})")
        print(f"    I(Ω, 2) = {sum(v * 2**k for k, v in alpha.items())}")

    # The key insight at p=7
    print(f"\n  KEY INSIGHT at p=7:")
    print(f"  Paley wins via α_1: 80*2 = 160 vs 59*2 = 118  (Δ = +42)")
    print(f"  Interval wins via α_2: 7*4 = 28 vs 14*4 = 56   (Δ = -28)")
    print(f"  Net: Paley +42-28 = +14 = 189-175 ✓")
    print(f"")
    print(f"  Paley has MORE total odd cycles (80 vs 59)")
    print(f"  but FEWER independent pairs (7 vs 14)")
    print(f"  Why? Paley cycles overlap MORE — the 'random' connectivity")
    print(f"  of QR creates dense overlap in the cycle graph.")
    print(f"  Interval cycles, being 'local' (small steps), overlap LESS.")

    # Part 2: p=11 approximate analysis
    print(f"\n{'='*75}")
    print(f"p = 11: APPROXIMATE ANALYSIS")
    print(f"{'='*75}")

    p = 11
    for name, S in [("Paley", paley_set(p)), ("Interval", interval_set(p))]:
        A = circulant_adj(p, S)
        H = hamiltonian_paths_dp(A, p)
        t0 = time.time()
        odd_cycles = find_odd_cycles(A, p)
        elapsed = time.time() - t0
        print(f"\n  {name}: S = {sorted(S)}, H = {H}")
        print(f"    {len(odd_cycles)} odd cycles ({elapsed:.2f}s)")
        len_dist = Counter(c[0] for c in odd_cycles)
        print(f"    Cycle counts: {dict(sorted(len_dist.items()))}")

        # Compute α_1, α_2 (α_3 may be slow)
        t0 = time.time()
        alpha = alpha_low(odd_cycles, p, max_k=2)
        elapsed = time.time() - t0
        for k in sorted(alpha.keys()):
            print(f"    α_{k} = {alpha[k]:>10}  (2^{k}·α_{k} = {alpha[k] * 2**k:>10})")
        H_approx = sum(v * 2**k for k, v in alpha.items())
        print(f"    Partial I(Ω, 2) from α_0,1,2 = {H_approx}")
        print(f"    Missing higher terms = {H - H_approx}")

    # Part 3: Interpretation
    print(f"\n{'='*75}")
    print("INTERPRETATION: THE CROSSOVER MECHANISM")
    print("=" * 75)

    print("""
  At p=7:
    Paley has 80 cycles, 7 disjoint pairs → dense overlap
    Interval has 59 cycles, 14 disjoint pairs → sparse overlap
    α_1 dominance: 2*Δα_1 > 4*Δα_2 → Paley wins

  At p=11:
    Both have ~21K cycles and many disjoint pairs
    The higher α_k terms (k≥3) involve MULTIPLE disjoint cycles
    Each disjoint cycle set gets 2^k weight
    Interval's LOCALITY creates more disjoint configurations:
      - Short-range cycles (local steps) are naturally disjoint
      - Paley's QR-based cycles span the whole vertex set, creating overlaps

  At p=19:
    The independence polynomial is dominated by terms with k ≥ 3
    (because p=19 allows many vertex-disjoint cycle packings)
    Interval's local structure creates exponentially more packings
    → Interval wins H

  GENERAL PRINCIPLE:
    Paley = "random" connectivity → many cycles but dense overlap
    Interval = "local" connectivity → fewer cycles but sparse overlap

    H = Σ α_k * 2^k is dominated by the LARGEST k terms
    Large k terms count vertex-disjoint cycle PACKINGS
    Local tournaments create better packings
    → At large p, local (interval) beats global (Paley)
    """)

    # Part 4: Verify the overlap density claim
    print(f"{'='*75}")
    print("PART 4: OVERLAP DENSITY (average neighbors per cycle in Ω)")
    print("=" * 75)

    for p in [7, 11]:
        print(f"\n  p = {p}:")
        for name, S in [("Paley", paley_set(p)), ("Interval", interval_set(p))]:
            A = circulant_adj(p, S)
            cycles = find_odd_cycles(A, p)
            m = len(cycles)

            # Average overlap: how many other cycles does each cycle intersect?
            total_adj = 0
            for i in range(m):
                for j in range(i + 1, m):
                    if cycles[i][1] & cycles[j][1]:
                        total_adj += 2  # both directions
            avg_deg = total_adj / m if m > 0 else 0
            density = total_adj / (m * (m - 1)) if m > 1 else 0

            print(f"    {name}: {m} cycles, avg overlap degree = {avg_deg:.2f}, "
                  f"overlap density = {density:.4f}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
