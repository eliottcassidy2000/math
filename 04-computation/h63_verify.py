"""
Search for tournaments on n=8 vertices with exactly H(T)=63 Hamiltonian paths.
Uses DP (bitmask) approach for counting Hamiltonian paths.
Samples random tournaments and checks structured ones.
"""

import numpy as np
import random
import itertools
import time

def count_hamiltonian_paths(adj, n):
    """Count Hamiltonian paths using DP with bitmask.
    dp[mask][v] = number of paths visiting exactly the vertices in mask, ending at v.
    """
    dp = [[0] * n for _ in range(1 << n)]
    # Base case: single vertex paths
    for v in range(n):
        dp[1 << v][v] = 1

    full = (1 << n) - 1
    for mask in range(1, full + 1):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]

    return sum(dp[full][v] for v in range(n))


def random_tournament(n):
    """Generate a random tournament on n vertices."""
    adj = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if random.random() < 0.5:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
    return adj


def transitive_tournament(n):
    """Transitive tournament: i beats j iff i < j."""
    adj = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            adj[i][j] = 1
    return adj


def circulant_tournament(n, S):
    """Circulant tournament C_n^S: vertex i beats j iff (j-i) mod n in S."""
    adj = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and (j - i) % n in S:
                adj[i][j] = 1
    return adj


def near_regular_tournament(n):
    """Start from a regular-ish tournament and perturb slightly."""
    # For even n, no regular tournament exists. Use circulant with half the residues.
    S = set(range(1, n // 2 + 1))
    return circulant_tournament(n, S)


def print_adj(adj, n):
    """Print adjacency matrix."""
    for i in range(n):
        print("  " + " ".join(str(adj[i][j]) for j in range(n)))


def adj_to_tuple(adj, n):
    """Convert adjacency to a hashable tuple (upper triangle)."""
    bits = []
    for i in range(n):
        for j in range(i + 1, n):
            bits.append(adj[i][j])
    return tuple(bits)


def main():
    n = 8
    target = 63

    print(f"Searching for tournaments on n={n} with H(T)={target}")
    print(f"=" * 60)

    found = []
    seen = set()

    # Check structured tournaments first
    print("\n--- Structured tournaments ---")

    # Transitive
    adj = transitive_tournament(n)
    h = count_hamiltonian_paths(adj, n)
    print(f"Transitive: H = {h}")
    if h == target:
        found.append(("Transitive", adj))

    # Circulant with various generating sets
    for S_tuple in itertools.combinations(range(1, n), n // 2):
        S = set(S_tuple)
        # Check it's a valid tournament (S and complement partition {1,...,n-1})
        complement = set(range(1, n)) - S
        # For tournament: if d in S then n-d must be in complement
        valid = all((n - d) % n not in S for d in S if (n - d) % n != 0)
        if not valid:
            continue
        adj = circulant_tournament(n, S)
        key = adj_to_tuple(adj, n)
        if key in seen:
            continue
        seen.add(key)
        h = count_hamiltonian_paths(adj, n)
        if h == target:
            print(f"Circulant S={S_tuple}: H = {h} *** FOUND ***")
            found.append((f"Circulant S={S_tuple}", adj))
        elif abs(h - target) <= 5:
            print(f"Circulant S={S_tuple}: H = {h} (close)")

    # Near-regular
    adj = near_regular_tournament(n)
    key = adj_to_tuple(adj, n)
    if key not in seen:
        seen.add(key)
        h = count_hamiltonian_paths(adj, n)
        print(f"Near-regular (circulant {{1..4}}): H = {h}")
        if h == target:
            found.append(("Near-regular", adj))

    # Random search
    print(f"\n--- Random search (up to 100,000 tournaments) ---")
    t0 = time.time()
    h_counts = {}
    num_random = 100_000

    for trial in range(num_random):
        adj = random_tournament(n)
        h = count_hamiltonian_paths(adj, n)
        h_counts[h] = h_counts.get(h, 0) + 1

        if h == target:
            key = adj_to_tuple(adj, n)
            if key not in seen:
                seen.add(key)
                found.append((f"Random trial {trial}", adj))
                if len(found) <= 5:
                    print(f"  Trial {trial}: FOUND H={target}!")

        if (trial + 1) % 20000 == 0:
            elapsed = time.time() - t0
            print(f"  ...{trial+1} checked ({elapsed:.1f}s)")

    elapsed = time.time() - t0
    print(f"  Completed {num_random} random trials in {elapsed:.1f}s")

    # Distribution summary
    print(f"\n--- H distribution from random sample ---")
    for h_val in sorted(h_counts.keys()):
        count = h_counts[h_val]
        pct = 100.0 * count / num_random
        marker = " <-- TARGET" if h_val == target else ""
        if count >= 10 or h_val == target or abs(h_val - target) <= 5:
            print(f"  H={h_val:4d}: {count:6d} ({pct:5.2f}%){marker}")

    # Results
    print(f"\n{'=' * 60}")
    print(f"Total tournaments found with H={target}: {len(found)}")

    if found:
        # Print first few
        for i, (label, adj) in enumerate(found[:3]):
            print(f"\n--- Example {i+1}: {label} ---")
            print(f"Adjacency matrix (adj[i][j]=1 means i beats j):")
            print_adj(adj, n)
            h = count_hamiltonian_paths(adj, n)
            print(f"Verification: H = {h}")
    else:
        print("\nNo tournament with H=63 found.")
        # Show closest values
        closest = sorted(h_counts.keys(), key=lambda x: abs(x - target))[:5]
        print(f"Closest H values observed: {closest}")
        for h_val in closest:
            print(f"  H={h_val}: {h_counts[h_val]} occurrences")


if __name__ == "__main__":
    random.seed(42)
    main()
