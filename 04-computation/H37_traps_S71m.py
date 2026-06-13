#!/usr/bin/env python3
"""
H=37 TRAPS AT n=6: WHAT MAKES THEM SPECIAL?
opus-2026-03-14-S71m

At n=6, there are 720 tournaments with H=37 that are local maxima
(no single arc flip increases H), even though the global max is H=45.

Questions:
1. What are the score sequences of H=37 traps?
2. How many arc flips needed to escape the trap?
3. What is the structure of the "valley" around these traps?
4. Are the traps related to specific tournament structures (e.g., near-regular)?
5. What is the relationship between H=37 traps and the H=45 global maxima?
6. Connection to the Morse theory on the hypercube
"""

from itertools import permutations, combinations
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

def score_seq(n, A):
    return tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))

def count_3cycles(n, A):
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j] + A[j][k] + A[k][i] == 3:
                    count += 1
                if A[i][k] + A[k][j] + A[j][i] == 3:
                    count += 1
    return count

print("=" * 70)
print("H=37 TRAPS AT n=6: WHAT MAKES THEM SPECIAL?")
print("opus-2026-03-14-S71m")
print("=" * 70)

n = 6
m = n * (n-1) // 2
total = 1 << m

# Compute all H values
print(f"\n  Computing all tournaments for n={n}...")
H_map = {}
for bits in range(total):
    A = adj_matrix(n, bits)
    H_map[bits] = count_hp(n, A)
print(f"  Done. {total} tournaments computed.")

H_vals = sorted(set(H_map.values()))
max_H = max(H_vals)

# Find local maxima
local_max = []
for bits in range(total):
    H = H_map[bits]
    is_max = True
    for arc in range(m):
        nbr = bits ^ (1 << arc)
        if H_map[nbr] > H:
            is_max = False
            break
    if is_max:
        local_max.append(bits)

# Separate into H=45 (global) and H=37 (traps)
global_max = [b for b in local_max if H_map[b] == max_H]
traps = [b for b in local_max if H_map[b] != max_H]
trap_H = sorted(set(H_map[b] for b in traps))

print(f"\n  Local maxima: {len(local_max)}")
print(f"    Global max (H={max_H}): {len(global_max)}")
print(f"    Traps: {len(traps)} at H values {trap_H}")

# ======================================================================
# ANALYSIS 1: Score sequences and 3-cycle counts
# ======================================================================
print(f"\n{'='*70}")
print("ANALYSIS 1: SCORE SEQUENCES AND 3-CYCLES")
print(f"{'='*70}")

# Score sequences of traps vs global maxima
trap_scores = Counter()
global_scores = Counter()
trap_t3 = Counter()
global_t3 = Counter()

for b in traps:
    A = adj_matrix(n, b)
    ss = score_seq(n, A)
    t3 = count_3cycles(n, A)
    trap_scores[ss] += 1
    trap_t3[t3] += 1

for b in global_max:
    A = adj_matrix(n, b)
    ss = score_seq(n, A)
    t3 = count_3cycles(n, A)
    global_scores[ss] += 1
    global_t3[t3] += 1

print(f"\n  Trap score sequences (H=37):")
for ss, count in sorted(trap_scores.items()):
    print(f"    {ss}: {count} tournaments")

print(f"\n  Global max score sequences (H=45):")
for ss, count in sorted(global_scores.items()):
    print(f"    {ss}: {count} tournaments")

print(f"\n  Trap 3-cycle counts:")
for t3, count in sorted(trap_t3.items()):
    print(f"    t3={t3}: {count} tournaments")

print(f"\n  Global max 3-cycle counts:")
for t3, count in sorted(global_t3.items()):
    print(f"    t3={t3}: {count} tournaments")

# ======================================================================
# ANALYSIS 2: Minimum flips to escape
# ======================================================================
print(f"\n{'='*70}")
print("ANALYSIS 2: MINIMUM FLIPS TO ESCAPE TRAP")
print(f"{'='*70}")

# For each trap, what is the minimum Hamming distance to a tournament
# with H > 37?
higher_H = {b for b in range(total) if H_map[b] > 37}

# Sample traps for distance calculation
sample_traps = traps[:100]

min_dists = []
for trap in sample_traps:
    min_d = m + 1
    for b in higher_H:
        d = bin(trap ^ b).count('1')
        if d < min_d:
            min_d = d
            if d == 2:
                break  # 1 is impossible (we know all single flips decrease)
    min_dists.append(min_d)

dist_counter = Counter(min_dists)
print(f"\n  Minimum flips to reach H > 37 from trap (sample of {len(sample_traps)}):")
for d, c in sorted(dist_counter.items()):
    print(f"    {d} flips: {c} traps")

# For the traps needing 2 flips, what does H look like along the way?
print(f"\n  2-flip escape paths (if any):")
escape_count = 0
for trap in sample_traps[:20]:
    if trap not in traps:
        continue
    # Try all 2-flip paths
    best_intermediate = None
    best_final = 0
    for a1 in range(m):
        mid = trap ^ (1 << a1)
        H_mid = H_map[mid]
        for a2 in range(a1+1, m):
            final = mid ^ (1 << a2)
            H_final = H_map[final]
            if H_final > 37:
                if best_final < H_final or best_intermediate is None:
                    best_intermediate = H_mid
                    best_final = H_final
    if best_intermediate is not None:
        escape_count += 1
        print(f"    Trap H=37 -> intermediate H={best_intermediate} -> H={best_final}")
        if escape_count >= 5:
            break

# ======================================================================
# ANALYSIS 3: Neighbor H distribution for traps
# ======================================================================
print(f"\n{'='*70}")
print("ANALYSIS 3: NEIGHBOR H DISTRIBUTION")
print(f"{'='*70}")

# For H=37 traps: what H values do their neighbors have?
trap_neighbor_H = Counter()
for trap in traps:
    for arc in range(m):
        nbr = trap ^ (1 << arc)
        trap_neighbor_H[H_map[nbr]] += 1

print(f"\n  Neighbors of H=37 traps (total: {sum(trap_neighbor_H.values())} edges):")
for h, count in sorted(trap_neighbor_H.items()):
    print(f"    H={h:3d}: {count:5d} ({100*count/sum(trap_neighbor_H.values()):5.1f}%)")

# Same for H=45 global maxima
global_neighbor_H = Counter()
for gm in global_max:
    for arc in range(m):
        nbr = gm ^ (1 << arc)
        global_neighbor_H[H_map[nbr]] += 1

print(f"\n  Neighbors of H=45 global max (total: {sum(global_neighbor_H.values())} edges):")
for h, count in sorted(global_neighbor_H.items()):
    print(f"    H={h:3d}: {count:5d} ({100*count/sum(global_neighbor_H.values()):5.1f}%)")

# ======================================================================
# ANALYSIS 4: Are traps self-complementary?
# ======================================================================
print(f"\n{'='*70}")
print("ANALYSIS 4: COMPLEMENT AND REVERSAL PROPERTIES")
print(f"{'='*70}")

# T^op (complement): flip all arcs
complement_mask = (1 << m) - 1

trap_set = set(traps)
global_set = set(global_max)

# Check if complements of traps are also traps
comp_also_trap = sum(1 for b in traps if (b ^ complement_mask) in trap_set)
comp_also_global = sum(1 for b in traps if (b ^ complement_mask) in global_set)

print(f"\n  Complement of H=37 trap is:")
print(f"    Also H=37 trap: {comp_also_trap}/{len(traps)}")
print(f"    H=45 global max: {comp_also_global}/{len(traps)}")

# What H does complement have?
comp_H = Counter(H_map[b ^ complement_mask] for b in traps)
print(f"    H distribution of complements: {dict(sorted(comp_H.items()))}")

# ======================================================================
# ANALYSIS 5: Which H=37 are NOT traps?
# ======================================================================
print(f"\n{'='*70}")
print("ANALYSIS 5: H=37 TRAPS vs NON-TRAPS")
print(f"{'='*70}")

all_H37 = [b for b in range(total) if H_map[b] == 37]
H37_not_trap = [b for b in all_H37 if b not in trap_set]

print(f"\n  Total H=37 tournaments: {len(all_H37)}")
print(f"  H=37 traps (local max): {len(traps)}")
print(f"  H=37 non-traps: {len(H37_not_trap)}")
print(f"  Fraction that are traps: {len(traps)/len(all_H37):.4f}")

# Score sequences of non-traps
nontrap_scores = Counter()
for b in H37_not_trap:
    A = adj_matrix(n, b)
    ss = score_seq(n, A)
    nontrap_scores[ss] += 1

print(f"\n  Non-trap H=37 score sequences:")
for ss, count in sorted(nontrap_scores.items()):
    print(f"    {ss}: {count} tournaments")

# For non-traps: what higher-H neighbor do they have?
if H37_not_trap:
    for b in H37_not_trap[:5]:
        higher_nbrs = []
        for arc in range(m):
            nbr = b ^ (1 << arc)
            if H_map[nbr] > 37:
                higher_nbrs.append((arc, H_map[nbr]))
        print(f"    Non-trap: has {len(higher_nbrs)} upward neighbors -> {[h for _, h in higher_nbrs]}")

# ======================================================================
# ANALYSIS 6: The H=37 "plateau" structure
# ======================================================================
print(f"\n{'='*70}")
print("ANALYSIS 6: IS THERE A PLATEAU?")
print(f"{'='*70}")

# How many H=37 neighbors does a typical H=37 trap have?
trap_to_trap = Counter()
trap_to_H37 = Counter()
for trap in traps[:200]:
    same_count = 0
    for arc in range(m):
        nbr = trap ^ (1 << arc)
        if H_map[nbr] == 37:
            same_count += 1
    trap_to_H37[same_count] += 1

print(f"\n  Number of H=37 neighbors of H=37 traps:")
for count, freq in sorted(trap_to_H37.items()):
    print(f"    {count} neighbors at H=37: {freq} traps")

# How connected is the H=37 trap subgraph?
print(f"\n  Connectivity of H=37 trap subgraph:")
trap_set_full = set(traps)
visited = set()
components = 0
for start in traps:
    if start in visited:
        continue
    components += 1
    queue = [start]
    visited.add(start)
    while queue:
        curr = queue.pop(0)
        for arc in range(m):
            nbr = curr ^ (1 << arc)
            if nbr in trap_set_full and nbr not in visited:
                visited.add(nbr)
                queue.append(nbr)

print(f"  {components} connected components among {len(traps)} traps")

# ======================================================================
# ANALYSIS 7: dH spectrum from traps
# ======================================================================
print(f"\n{'='*70}")
print("ANALYSIS 7: dH VALUES FROM TRAPS")
print(f"{'='*70}")

# For each trap, what are the possible dH values?
dH_from_trap = Counter()
for trap in traps:
    for arc in range(m):
        nbr = trap ^ (1 << arc)
        dH = H_map[nbr] - 37
        dH_from_trap[dH] += 1

print(f"\n  dH distribution from H=37 traps:")
for dh, count in sorted(dH_from_trap.items()):
    pct = 100 * count / sum(dH_from_trap.values())
    print(f"    dH={dh:+4d}: {count:6d} ({pct:5.1f}%)")

print(f"\n  Key: the MAXIMUM dH from any trap = {max(dH_from_trap.keys())}")
print(f"  This means no single flip can reach H > {37 + max(dH_from_trap.keys())}")
print(f"  from a trap. Traps are surrounded by H ≤ {37 + max(dH_from_trap.keys())}.")

# ======================================================================
# ANALYSIS 8: What does the H landscape look like?
# ======================================================================
print(f"\n{'='*70}")
print("ANALYSIS 8: LANDSCAPE STRUCTURE")
print(f"{'='*70}")

# Count tournaments at each H level
H_counter = Counter(H_map.values())
print(f"\n  H distribution:")
for h in sorted(H_counter.keys()):
    count = H_counter[h]
    bar = '#' * (count // 20)
    print(f"    H={h:3d}: {count:5d} {bar}")

# For each H level, count local max, local min, saddle
print(f"\n  Critical point analysis by H level:")
for h_target in sorted(H_counter.keys()):
    level_set = [b for b in range(total) if H_map[b] == h_target]
    loc_max = 0
    loc_min = 0
    for b in level_set:
        higher = any(H_map[b ^ (1<<a)] > H_map[b] for a in range(m))
        lower = any(H_map[b ^ (1<<a)] < H_map[b] for a in range(m))
        if not higher:
            loc_max += 1
        if not lower:
            loc_min += 1
    if loc_max > 0 or loc_min > 0:
        print(f"    H={h_target:3d}: {H_counter[h_target]:5d} total, "
              f"{loc_max:4d} local max, {loc_min:4d} local min")

print("\n" + "=" * 70)
print("DONE — H=37 TRAP ANALYSIS")
print("=" * 70)
