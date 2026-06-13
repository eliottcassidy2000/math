#!/usr/bin/env python3
"""
Hamiltonian cycle threshold at n=6.
opus-2026-03-14-S84

At n=5: H≥9 iff tournament has a directed Hamiltonian cycle (HC).
H∈{1,3,5}: 0% HC. H∈{9,11,13,15}: 100% HC.
Threshold is exactly H=9 = KEY2² = 3².

Does a similar threshold exist at n=6?
By Camion's theorem: strongly connected tournament iff has HC.
So the question becomes: at n=6, which H values have ALL tournaments
strongly connected?

Also: explore the connection between #HCs and H.
"""

from itertools import permutations, combinations
from collections import Counter, defaultdict
import sys

def compute_all_data_n6():
    """Compute H, HC existence, #SCCs for all n=6 tournaments."""
    n = 6
    m = 15
    N = 1 << m
    arcs = [(i, j) for i in range(n) for j in range(i+1, n)]
    all_perms = list(permutations(range(n)))

    data = []
    for bits in range(N):
        if bits % 5000 == 0:
            print(f"  {bits}/{N} ({100*bits/N:.1f}%)", file=sys.stderr)

        adj = [[0]*n for _ in range(n)]
        for k, (i, j) in enumerate(arcs):
            if (bits >> k) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1

        # Count Hamiltonian paths
        H = 0
        for p in all_perms:
            valid = True
            for i in range(n-1):
                if adj[p[i]][p[i+1]] != 1:
                    valid = False
                    break
            if valid:
                H += 1

        # Check for Hamiltonian cycle
        has_hc = False
        for p in all_perms:
            if all(adj[p[i]][p[(i+1)%n]] == 1 for i in range(n)):
                has_hc = True
                break

        # Count SCCs using simple reachability
        def reachable(start, adj, n):
            visited = {start}
            queue = [start]
            while queue:
                v = queue.pop(0)
                for w in range(n):
                    if adj[v][w] and w not in visited:
                        visited.add(w)
                        queue.append(w)
            return visited

        # Strongly connected iff every vertex reachable from every other
        sc = len(reachable(0, adj, n)) == n

        data.append((H, has_hc, sc))

    return data

print("Computing all n=6 tournaments...")
data = compute_all_data_n6()

# ============================================================
# Part 1: HC existence by H value
# ============================================================
print("\n" + "=" * 70)
print("PART 1: HAMILTONIAN CYCLE EXISTENCE BY H VALUE")
print("=" * 70)

hc_by_H = defaultdict(lambda: [0, 0])  # [with_hc, without_hc]
sc_by_H = defaultdict(lambda: [0, 0])  # [sc, not_sc]

for H, has_hc, sc in data:
    hc_by_H[H][0 if has_hc else 1] += 1
    sc_by_H[H][0 if sc else 1] += 1

print(f"\n{'H':>4} {'with_HC':>8} {'no_HC':>8} {'%HC':>8}  {'SC':>8} {'not_SC':>8} {'%SC':>8}")
for h in sorted(hc_by_H.keys()):
    whc, nhc = hc_by_H[h]
    wsc, nsc = sc_by_H[h]
    total = whc + nhc
    pct_hc = 100 * whc / total
    pct_sc = 100 * wsc / total
    print(f"{h:4d} {whc:8d} {nhc:8d} {pct_hc:7.1f}%  {wsc:8d} {nsc:8d} {pct_sc:7.1f}%")

# ============================================================
# Part 2: Threshold analysis
# ============================================================
print("\n" + "=" * 70)
print("PART 2: THRESHOLD ANALYSIS")
print("=" * 70)

# Find the H value where ALL tournaments have HC
threshold_hc = None
for h in sorted(hc_by_H.keys()):
    whc, nhc = hc_by_H[h]
    if nhc == 0:
        if threshold_hc is None:
            threshold_hc = h
    else:
        threshold_hc = None  # reset if not all have HC

print(f"Lowest H where ALL tournaments have HC: {threshold_hc}")

# Find where HC first appears
first_hc = None
for h in sorted(hc_by_H.keys()):
    if hc_by_H[h][0] > 0:
        first_hc = h
        break
print(f"Lowest H where HC exists: {first_hc}")

# Find where ALL are strongly connected
threshold_sc = None
for h in sorted(sc_by_H.keys()):
    if sc_by_H[h][1] == 0:
        if threshold_sc is None:
            threshold_sc = h
    else:
        threshold_sc = None

print(f"Lowest H where ALL are strongly connected: {threshold_sc}")

# ============================================================
# Part 3: Number of Hamiltonian cycles vs H
# ============================================================
print("\n" + "=" * 70)
print("PART 3: NUMBER OF HAMILTONIAN CYCLES BY H")
print("=" * 70)

n = 6
arcs = [(i, j) for i in range(n) for j in range(i+1, n)]
all_perms = list(permutations(range(n)))

# Count HCs for a sample of tournaments at each H value
# Full HC count for all 32768 is expensive. Let's do H count for
# tournaments grouped by H value.

# Actually let me count HCs directly during the main loop.
# Re-do with HC counts (lighter weight than full recompute)
print("\nRe-computing with HC counts...")
hc_count_by_H = defaultdict(list)

for bits in range(1 << 15):
    if bits % 5000 == 0:
        print(f"  {bits}/32768", file=sys.stderr)

    adj = [[0]*n for _ in range(n)]
    for k, (i, j) in enumerate(arcs):
        if (bits >> k) & 1:
            adj[i][j] = 1
        else:
            adj[j][i] = 1

    H = sum(1 for p in all_perms if all(adj[p[i]][p[i+1]] == 1 for i in range(n-1)))

    # Count directed Hamiltonian cycles (each cycle counted n times for rotations)
    hc_raw = sum(1 for p in all_perms if all(adj[p[i]][p[(i+1)%n]] == 1 for i in range(n)))
    # Divide by n for rotation equivalence
    hc_count = hc_raw // n

    hc_count_by_H[H].append(hc_count)

print(f"\n{'H':>4} {'mean_HC':>8} {'min_HC':>7} {'max_HC':>7} {'HC=0':>6} {'total':>6}")
for h in sorted(hc_count_by_H.keys()):
    vals = hc_count_by_H[h]
    mean_hc = sum(vals) / len(vals)
    min_hc = min(vals)
    max_hc = max(vals)
    zero_count = vals.count(0)
    print(f"{h:4d} {mean_hc:8.2f} {min_hc:7d} {max_hc:7d} {zero_count:6d} {len(vals):6d}")

# ============================================================
# Part 4: HC/HP ratio
# ============================================================
print("\n" + "=" * 70)
print("PART 4: HC/HP RATIO = PROBABILITY OF CLOSING THE PATH")
print("=" * 70)

# For a tournament with H Hamiltonian paths and C Hamiltonian cycles,
# the ratio C*n/H is the "cycle closure probability":
# given a random HP, what fraction of them close to a cycle?

for h in sorted(hc_count_by_H.keys()):
    vals = hc_count_by_H[h]
    if h > 0:
        ratios = [v * n / h if h > 0 else 0 for v in vals]
        mean_ratio = sum(ratios) / len(ratios)
        print(f"  H={h:2d}: mean(C*n/H) = {mean_ratio:.4f} ({mean_ratio*100:.1f}%)")

# ============================================================
# SYNTHESIS
# ============================================================
print("\n" + "=" * 70)
print("SYNTHESIS — HC THRESHOLD")
print("=" * 70)
print(f"""
At n=5: HC threshold = H=9 = 3² (SHARP: all H≥9 have HC, all H<9 don't)
At n=6: Threshold behavior analyzed above.

Key questions answered:
1. Is there a sharp threshold at n=6?
2. How does the number of HCs scale with H?
3. What fraction of HPs can be closed to cycles?
""")
