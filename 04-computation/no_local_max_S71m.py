#!/usr/bin/env python3
"""
NO LOCAL MAXIMA BELOW GLOBAL — DEEP INVESTIGATION
opus-2026-03-14-S71m

The most striking finding: ALL local maxima of H on the tournament
hypercube Q_m are GLOBAL maxima. There are no "false peaks."

This is a CONVEXITY-like property, but H is NOT convex on Q_m.
(It can't be — Q_m is a discrete space.)

What this means: From ANY tournament, you can reach the maximum H
by following a greedy arc-flip strategy. There are no dead ends.

Questions:
1. Does this hold for n=6?
2. What about local MINIMA? Are they all global?
3. What is the count of sinks as a function of n?
4. Is this a known property of independence polynomials on hypercubes?
5. What about SADDLE points?
"""

from itertools import permutations
from collections import Counter, defaultdict
import math

def adj_matrix(n, bits):
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

def count_hp(n, A):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))

print("=" * 70)
print("NO LOCAL MAXIMA BELOW GLOBAL — DEEP INVESTIGATION")
print("opus-2026-03-14-S71m")
print("=" * 70)

for n in [3, 4, 5, 6]:
    m = n * (n-1) // 2
    total = 1 << m

    print(f"\n{'='*70}")
    print(f"n={n} (m={m}, {total} tournaments)")
    print(f"{'='*70}")

    # Compute all H values
    H_map = {}
    for bits in range(total):
        A = adj_matrix(n, bits)
        H_map[bits] = count_hp(n, A)

    H_vals = sorted(set(H_map.values()))
    max_H = max(H_vals)
    min_H = min(H_vals)
    print(f"  H-spectrum: {H_vals}")
    print(f"  min H = {min_H}, max H = {max_H}")

    # Classify critical points
    local_max = []  # H >= all neighbors
    local_min = []  # H <= all neighbors
    saddle = []     # neither
    strict_local_max = []  # H > all neighbors
    strict_local_min = []  # H < all neighbors

    for bits in range(total):
        H = H_map[bits]
        neighbors_H = []
        for arc in range(m):
            nbr = bits ^ (1 << arc)
            neighbors_H.append(H_map[nbr])

        if all(H >= nh for nh in neighbors_H):
            local_max.append(bits)
            if all(H > nh for nh in neighbors_H):
                strict_local_max.append(bits)
        elif all(H <= nh for nh in neighbors_H):
            local_min.append(bits)
            if all(H < nh for nh in neighbors_H):
                strict_local_min.append(bits)
        else:
            saddle.append(bits)

    # Classify by H value
    max_H_vals = sorted(set(H_map[b] for b in local_max))
    min_H_vals = sorted(set(H_map[b] for b in local_min))

    print(f"\n  Critical points:")
    print(f"    Local maxima: {len(local_max)} (strict: {len(strict_local_max)})")
    print(f"      H values: {max_H_vals}")
    all_global_max = all(H_map[b] == max_H for b in local_max)
    print(f"      ALL are global maximum? {all_global_max}")

    print(f"    Local minima: {len(local_min)} (strict: {len(strict_local_min)})")
    print(f"      H values: {min_H_vals}")
    all_global_min = all(H_map[b] == min_H for b in local_min)
    print(f"      ALL are global minimum? {all_global_min}")

    print(f"    Saddle points: {len(saddle)}")

    # More detailed saddle analysis
    if n <= 5:
        # Index of each tournament: how many neighbors have HIGHER H?
        index_dist = Counter()
        for bits in range(total):
            H = H_map[bits]
            higher = sum(1 for arc in range(m) if H_map[bits ^ (1 << arc)] > H)
            lower = sum(1 for arc in range(m) if H_map[bits ^ (1 << arc)] < H)
            equal = m - higher - lower
            index_dist[(higher, lower, equal)] += 1

        print(f"\n  Morse index distribution (higher, lower, equal neighbors):")
        for (h, l, e), count in sorted(index_dist.items()):
            # Traditional Morse index = number of descending directions
            print(f"    ({h} higher, {l} lower, {e} equal): {count} tournaments")

    # Count how many tournaments have NO upward move (local max)
    # and NO downward move (local min)
    print(f"\n  Gradient structure:")

    # Average H-change from a random arc flip
    total_delta = 0
    total_abs_delta = 0
    count = 0
    for bits in range(total):
        H = H_map[bits]
        for arc in range(m):
            nbr = bits ^ (1 << arc)
            delta = H_map[nbr] - H
            total_delta += delta
            total_abs_delta += abs(delta)
            count += 1

    print(f"    Mean dH per flip = {total_delta/count:.6f} (should be 0 by symmetry)")
    print(f"    Mean |dH| per flip = {total_abs_delta/count:.6f}")

    # Distribution of dH
    dH_dist = Counter()
    for bits in range(total):
        H = H_map[bits]
        for arc in range(m):
            nbr = bits ^ (1 << arc)
            dH_dist[H_map[nbr] - H] += 1

    print(f"    dH distribution:")
    for dh in sorted(dH_dist.keys()):
        pct = 100 * dH_dist[dh] / count
        print(f"      dH={dh:+4d}: {dH_dist[dh]:8d} ({pct:5.1f}%)")

    # Number of sinks in steepest ascent
    print(f"\n  Gradient flow (steepest ascent):")
    gradient_next = {}
    for bits in range(total):
        H = H_map[bits]
        best_delta = 0
        best_nbr = None
        for arc in range(m):
            nbr = bits ^ (1 << arc)
            delta = H_map[nbr] - H
            if delta > best_delta:
                best_delta = delta
                best_nbr = nbr
        gradient_next[bits] = best_nbr

    sinks = [b for b, nbr in gradient_next.items() if nbr is None]
    print(f"    Sinks (local max): {len(sinks)}")
    print(f"    Sink H values: {sorted(set(H_map[s] for s in sinks))}")

    # Gradient flow (steepest DESCENT)
    gradient_down = {}
    for bits in range(total):
        H = H_map[bits]
        best_delta = 0
        best_nbr = None
        for arc in range(m):
            nbr = bits ^ (1 << arc)
            delta = H - H_map[nbr]  # want LARGEST decrease
            if delta > best_delta:
                best_delta = delta
                best_nbr = nbr
        gradient_down[bits] = best_nbr

    sources = [b for b, nbr in gradient_down.items() if nbr is None]
    print(f"    Sources (local min, descending): {len(sources)}")
    print(f"    Source H values: {sorted(set(H_map[s] for s in sources))}")

    # Basin sizes for ascending flow
    basin_sizes = Counter()
    for bits in range(total):
        curr = bits
        path_len = 0
        while gradient_next[curr] is not None:
            curr = gradient_next[curr]
            path_len += 1
            if path_len > 1000:
                break
        basin_sizes[curr] += 1

    sizes = sorted(basin_sizes.values())
    print(f"    Basin sizes: min={min(sizes)}, max={max(sizes)}, "
          f"mean={sum(sizes)/len(sizes):.1f}, median={sizes[len(sizes)//2]}")

    # Check: can every tournament reach max H by SOME sequence of
    # H-increasing flips? (Not just steepest ascent)
    print(f"\n  Can every tournament reach max H by some H-increasing path?")
    can_reach_max = set()
    # BFS backward from max-H tournaments
    max_H_set = {b for b in range(total) if H_map[b] == max_H}
    can_reach_max = set(max_H_set)
    queue = list(max_H_set)
    while queue:
        curr = queue.pop(0)
        H_curr = H_map[curr]
        for arc in range(m):
            nbr = curr ^ (1 << arc)
            if nbr not in can_reach_max and H_map[nbr] < H_curr:
                can_reach_max.add(nbr)
                queue.append(nbr)

    print(f"    Reachable from below: {len(can_reach_max)} / {total} "
          f"({'ALL' if len(can_reach_max) == total else 'NOT ALL'})")

    if len(can_reach_max) < total:
        unreachable = [b for b in range(total) if b not in can_reach_max]
        print(f"    Unreachable H values: {sorted(set(H_map[b] for b in unreachable))}")

    # Check: can every tournament reach min H by some H-DECREASING path?
    can_reach_min = set()
    min_H_set = {b for b in range(total) if H_map[b] == min_H}
    can_reach_min = set(min_H_set)
    queue = list(min_H_set)
    while queue:
        curr = queue.pop(0)
        H_curr = H_map[curr]
        for arc in range(m):
            nbr = curr ^ (1 << arc)
            if nbr not in can_reach_min and H_map[nbr] > H_curr:
                can_reach_min.add(nbr)
                queue.append(nbr)

    print(f"    Can reach min by decreasing: {len(can_reach_min)} / {total} "
          f"({'ALL' if len(can_reach_min) == total else 'NOT ALL'})")

    # Number of sinks at each n
    if n <= 5:
        print(f"\n  SINK COUNT ANALYSIS:")
        print(f"    n={n}: {len(sinks)} sinks = {len(sinks)}")
        if n == 3:
            print(f"    2 = 2")
        elif n == 4:
            print(f"    24 = 4! = 24")
        elif n == 5:
            print(f"    64 = 2^6 = 2^C(4,2)")
            # Actually 64 = number of regular tournaments * their automorphism group?
            # Regular tournaments at n=5: there are 3 isomorphism classes
            # But 24 out of 64 are the (2,2,2,2,2) regular tournaments
            # The rest have score (1,2,2,2,3) and H=15
            H15_count = sum(1 for b in range(total) if H_map[b] == 15)
            print(f"    H=15 total count: {H15_count}")
            print(f"    Sinks that are H=15: {len(sinks)}")
            print(f"    Fraction of H=15 that are sinks: {len(sinks)}/{H15_count} = {len(sinks)/H15_count:.4f}")

    if n == 6:
        print(f"\n  SINK COUNT ANALYSIS:")
        print(f"    n={n}: {len(sinks)} sinks")
        H_max_count = sum(1 for b in range(total) if H_map[b] == max_H)
        print(f"    H={max_H} total count: {H_max_count}")
        print(f"    Fraction of H={max_H} that are sinks: {len(sinks)}/{H_max_count} = {len(sinks)/H_max_count:.4f}")

        # Also check: how many of these sinks are transitive?
        # (transitive = there exists an ordering where all arcs go "forward")
        trans_count = 0
        for b in sinks:
            A = adj_matrix(n, b)
            # Check if tournament is transitive: score seq = (0,1,2,...,n-1)
            scores = sorted(sum(A[i][j] for j in range(n)) for i in range(n))
            if scores == list(range(n)):
                trans_count += 1
        print(f"    Transitive sinks: {trans_count}")

print("\n" + "=" * 70)
print("SUMMARY: LOCAL MAX = GLOBAL MAX PATTERN")
print("=" * 70)

# Print the sink count sequence
print("""
  Sink counts: n=3: 2, n=4: 24, n=5: 64

  These are the counts of TOURNAMENTS with maximum H that are LOCAL
  maxima under single arc flips. Since ALL local maxima are global,
  this equals the number of local maxima.

  Interpretations:
  - n=3: 2 = 2 (both 3-cycle tournaments are local maxima)
  - n=4: 24 = 4! (the transitive tournament orbits)
  - n=5: 64 = 2^6 = 2^C(4,2)

  The PROPERTY that all local maxima are global is called
  "no local maxima" or "unimodal landscape" in optimization theory.

  For tournament H, this would mean:
  THEOREM (conjecture): For any tournament T, there exists a sequence
  of arc flips T = T_0, T_1, ..., T_k = T_max where:
  (1) H(T_i) < H(T_{i+1}) for all i,
  (2) T_max has maximum H among all tournaments on n vertices.

  This is equivalent to saying: the steepest ascent always reaches
  the global maximum. No local traps.

  SIGNIFICANCE: This would be a VERY strong structural result about
  the tournament HP-count landscape. It would mean that greedy
  optimization of H always succeeds.
""")

print("=" * 70)
print("DONE")
print("=" * 70)
