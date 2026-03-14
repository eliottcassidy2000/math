#!/usr/bin/env python3
"""
alpha_ratio_n7_fast.py — opus-2026-03-14-S75

Fast version: only count 3-cycles for CG(T) at n=7.
5-cycles are rare and expensive to enumerate.
At n=7, most cycles are 3-cycles anyway.

Verify α₁ ≥ α₂ using only 3-cycle conflict graph.
Since adding 5-cycles would only INCREASE α₁ and might increase or
decrease α₂, this is a partial verification.

Actually: since 5-cycles are additional vertices in CG(T),
ignoring them gives FEWER vertices → LOWER α₁.
But the disjointness structure is complex.
So this gives a LOWER BOUND on α₁ and might miss some α₂ pairs.

Better approach: use the OCF formula directly.
H = 1 + 2·dc3 + 4·(dc3_pairs) for n≤8 (where dc5 pairs don't exist at n=7
since 5+5=10>7, so all disjoint pairs are (3,3) or (3,5) or (5,5)).
Wait: at n=7, we can have (3-cycle, 3-cycle) pairs using 6 vertices,
or (3-cycle, no room for 5-cycle since 3+5=8>7).
So at n=7: α₂ = number of disjoint 3-cycle pairs only!
No (3,5) or (5,5) disjoint pairs are possible.

Similarly: α₃ = number of disjoint 3-cycle triples.
But 3+3+3=9>7, so α₃ = 0 at n=7!

CONCLUSION at n=7:
- α₁ = dc3 + dc5 + dc7 (3-cycles + 5-cycles + 7-cycles)
- α₂ = number of disjoint (3,3) pairs only (since 3+5=8>7)
- α₃ = 0

So we only need to count 3-cycles and their disjoint pairs!
"""

from itertools import combinations
import sys

n = 7
edges = [(i, j) for i in range(n) for j in range(i+1, n)]
num_edges = len(edges)  # 21
total = 2 ** num_edges  # 2097152

print(f"n={n}: {total} tournaments, {num_edges} edges")
print(f"At n=7: α₃=0 (9>7), disjoint pairs only from 3-cycles (3+5=8>7)")
print()

max_ratio = 0
max_ratio_info = None
violations = 0
alpha_dist = {}
im1_dist = {}
h_alpha_pairs = {}

chunk = total // 10

for bits in range(total):
    if bits % chunk == 0:
        pct = 100 * bits // total
        print(f"  Progress: {bits}/{total} ({pct}%), violations: {violations}, max α₂/α₁: {max_ratio:.4f}", flush=True)

    # Build adjacency
    adj = [0] * n  # bitmask adjacency
    for idx, (i, j) in enumerate(edges):
        if bits & (1 << idx):
            adj[i] |= (1 << j)
        else:
            adj[j] |= (1 << i)

    # Count 3-cycles
    cycles_3 = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                # Check if i→j→k→i or i→k→j→i (one of the two cyclic orders)
                if (adj[i] & (1<<j)) and (adj[j] & (1<<k)) and (adj[k] & (1<<i)):
                    cycles_3.append((i, j, k))
                elif (adj[i] & (1<<k)) and (adj[k] & (1<<j)) and (adj[j] & (1<<i)):
                    cycles_3.append((i, j, k))

    dc3 = len(cycles_3)

    # Count disjoint 3-cycle pairs
    # Two 3-cycles are disjoint iff they share no vertex
    # At n=7, a disjoint pair uses 6 vertices, leaving 1 unused
    dc3_pairs = 0
    for a in range(dc3):
        for b in range(a+1, dc3):
            v1 = set(cycles_3[a])
            v2 = set(cycles_3[b])
            if len(v1 & v2) == 0:
                dc3_pairs += 1

    # At n=7: α₁ = dc3 + (5-cycles) + (7-cycles)
    # But we're only counting 3-cycles, so α₁_lower = dc3
    # α₂ = dc3_pairs (exact, since no other disjoint pairs possible)
    # For the inequality α₁ ≥ α₂:
    # We check dc3 ≥ dc3_pairs (lower bound on α₁ vs exact α₂)
    # If this holds, then α₁ ≥ dc3 ≥ dc3_pairs = α₂ ✓

    a1_lower = dc3  # actual α₁ ≥ dc3
    a2_exact = dc3_pairs

    if a2_exact > a1_lower:
        # This means dc3 < dc3_pairs, but actual α₁ might still be ≥ α₂
        # Count 5-cycles to get actual α₁
        violations += 1
        if violations <= 10:
            # Count 5-cycles for this tournament
            dc5 = 0
            for combo in combinations(range(n), 5):
                verts = list(combo)
                # Check all 12 directed 5-cycles on these vertices
                from itertools import permutations as perms
                found = False
                for perm in perms(verts):
                    if found:
                        break
                    is_cycle = True
                    for idx in range(5):
                        if not (adj[perm[idx]] & (1 << perm[(idx+1) % 5])):
                            is_cycle = False
                            break
                    if is_cycle:
                        # Check chordless
                        chordless = True
                        for idx in range(5):
                            v1 = perm[idx]
                            v2 = perm[(idx+2) % 5]
                            if adj[v1] & (1 << v2):
                                chordless = False
                                break
                        if chordless:
                            dc5 += 1
                            found = True

            actual_a1 = dc3 + dc5  # ignoring 7-cycles
            print(f"    dc3={dc3} < dc3_pairs={dc3_pairs}, but dc5={dc5}, actual α₁≥{actual_a1}, α₂={a2_exact}, ok? {actual_a1 >= a2_exact}")

    if dc3 > 0:
        ratio = a2_exact / a1_lower
        if ratio > max_ratio:
            max_ratio = ratio
            max_ratio_info = (bits, dc3, a2_exact)

    key = (dc3, a2_exact)
    alpha_dist[key] = alpha_dist.get(key, 0) + 1

    im1 = 1 - dc3 + a2_exact  # approximate I(-1) ignoring 5/7-cycles
    im1_dist[im1] = im1_dist.get(im1, 0) + 1

print()
print(f"RESULTS for n={n} (3-cycle only analysis):")
print(f"  Total tournaments: {total}")
print(f"  dc3 < dc3_pairs ('violations'): {violations}")
print(f"  Max ratio dc3_pairs/dc3: {max_ratio:.6f}")
if max_ratio_info:
    b, a1, a2 = max_ratio_info
    print(f"    Achieved at bits={b}, dc3={a1}, dc3_pairs={a2}")

print()
print("  (dc3, dc3_pairs) distribution (top 30 by frequency):")
sorted_dist = sorted(alpha_dist.items(), key=lambda x: -x[1])
for (a1, a2), count in sorted_dist[:30]:
    ratio = a2/a1 if a1 > 0 else 0
    im1 = 1 - a1 + a2
    print(f"    dc3={a1:3d}, pairs={a2:3d}: {count:8d} tours, ratio={ratio:.4f}, I(-1)_approx={im1}")

print()
print("  I(-1) approximate distribution (3-cycle only):")
for im1 in sorted(im1_dist.keys()):
    if im1_dist[im1] > 1000:  # only show significant
        print(f"    I(-1) ≈ {im1:4d}: {im1_dist[im1]:8d} tournaments")

print()
print("  Summary:")
print(f"    Max dc3: {max(k[0] for k in alpha_dist.keys())}")
print(f"    Max dc3_pairs: {max(k[1] for k in alpha_dist.keys())}")
print(f"    dc3 < pairs violations: {violations}")
if violations == 0:
    print("    CONFIRMED: dc3 ≥ dc3_pairs for ALL n=7 tournaments")
    print("    Since α₁ ≥ dc3 and α₂ = dc3_pairs, this implies α₁ ≥ α₂ ✓")
