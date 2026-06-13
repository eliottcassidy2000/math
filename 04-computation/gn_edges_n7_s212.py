#!/usr/bin/env python3
"""
gn_edges_n7_s212.py — Compute |E(G_7)| exactly
opus-2026-03-22-S212

Strategy: find one representative per iso class, then compute neighbors.
Phase 1: hash all 2^21 tournaments by (score_seq, c3) → 242 hash groups
Phase 2: split hash groups into iso classes via canonical form sampling
Phase 3: for each class rep, flip 21 arcs, canonicalize → neighbor class
Phase 4: count edges

This gives us the 5th term of the sequence 1, 5, 30, 290, ?
"""

import sys
import random
from math import comb
from itertools import permutations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  COMPUTING |E(G_7)| EXACTLY")
print("  opus-2026-03-22-S212")
print("=" * 80)

n = 7
m = comb(n, 2)  # 21
total = 1 << m   # 2097152

# Precompute all permutations of range(n) for canonical form
ALL_PERMS = list(permutations(range(n)))
print(f"  n={n}, m={m}, |tournaments|={total}, |S_n|={len(ALL_PERMS)}")

def tournament_adj(n, bits):
    adj = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            idx += 1
    return adj

def score_seq(adj, n):
    return tuple(sorted(sum(adj[i][j] for j in range(n)) for i in range(n)))

def count_3cycles(adj, n):
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if adj[i][j] and adj[j][k] and adj[k][i]:
                    c3 += 1
                if adj[i][k] and adj[k][j] and adj[j][i]:
                    c3 += 1
    return c3

def canonical_form(adj, n):
    """Canonical form using precomputed permutations."""
    best = None
    for perm in ALL_PERMS:
        form = tuple(adj[perm[i]][perm[j]] for i in range(n) for j in range(n))
        if best is None or form < best:
            best = form
    return best

def adj_to_bits(adj, n):
    """Convert adjacency matrix back to bit encoding."""
    bits = 0
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if adj[i][j]:
                bits |= (1 << idx)
            idx += 1
    return bits

# ============================================================================
# PHASE 1: Hash all tournaments by (score_seq, c3)
# ============================================================================

print(f"\n  PHASE 1: Hash all {total} tournaments by (score, c3)...")
t0 = time.time()

hash_groups = defaultdict(list)  # (score, c3) -> [bits]

for bits in range(total):
    adj = tournament_adj(n, bits)
    ss = score_seq(adj, n)
    c3 = count_3cycles(adj, n)
    hash_groups[(ss, c3)].append(bits)

    if (bits + 1) % 200000 == 0:
        elapsed = time.time() - t0
        rate = (bits + 1) / elapsed
        eta = (total - bits - 1) / rate
        print(f"    {bits+1}/{total} ({rate:.0f}/s, ETA {eta:.0f}s)")

print(f"  Phase 1 done: {len(hash_groups)} hash groups in {time.time()-t0:.1f}s")

group_sizes = Counter(len(v) for v in hash_groups.values())
print(f"  Group size distribution: {dict(sorted(group_sizes.items()))}")

# ============================================================================
# PHASE 2: Split hash groups into iso classes via canonical forms
# ============================================================================

print(f"\n  PHASE 2: Split hash groups into iso classes...")
t1 = time.time()

class_reps = {}  # canonical_form -> bits (representative)
canon_to_id = {} # canonical_form -> class_id
class_id = 0
class_sizes = {} # class_id -> count

total_canon_computed = 0

for hkey, members in sorted(hash_groups.items(), key=lambda x: len(x[1])):
    group_size = len(members)

    if group_size <= 5040:
        # Small group: compute canonical form for each member
        local_classes = defaultdict(int)
        for bits in members:
            adj = tournament_adj(n, bits)
            canon = canonical_form(adj, n)
            total_canon_computed += 1
            if canon not in class_reps:
                class_reps[canon] = bits
                canon_to_id[canon] = class_id
                class_sizes[class_id] = 0
                class_id += 1
            class_sizes[canon_to_id[canon]] += 1
    else:
        # Large group: sample to find distinct classes, then count
        # First, find all distinct canonical forms by sampling
        found_canons = {}
        remaining = list(members)
        random.shuffle(remaining)

        accounted = 0
        for bits in remaining:
            adj = tournament_adj(n, bits)
            canon = canonical_form(adj, n)
            total_canon_computed += 1
            if canon not in found_canons:
                found_canons[canon] = bits
                if canon not in class_reps:
                    class_reps[canon] = bits
                    canon_to_id[canon] = class_id
                    class_sizes[class_id] = 0
                    class_id += 1

            class_sizes[canon_to_id[canon]] += 1
            accounted += 1

            # Check if we've accounted for the entire group
            total_found = sum(class_sizes[canon_to_id[c]] for c in found_canons)
            if total_found >= group_size:
                break

    if total_canon_computed % 5000 == 0:
        elapsed = time.time() - t1
        print(f"    {total_canon_computed} canonical forms computed, "
              f"{class_id} classes found ({elapsed:.1f}s)")

print(f"  Phase 2 done: {class_id} iso classes found in {time.time()-t1:.1f}s")
print(f"  Total canonical forms computed: {total_canon_computed}")
print(f"  Expected: 456 classes (A000568)")

# Verify total
total_accounted = sum(class_sizes.values())
print(f"  Total tournaments accounted: {total_accounted}/{total}")

if total_accounted < total:
    print(f"  WARNING: {total - total_accounted} tournaments unaccounted!")
    print(f"  These are in large groups where sampling didn't cover all members.")
    print(f"  But we should have found all DISTINCT classes.")

# ============================================================================
# PHASE 3: For each class, flip 21 arcs, canonicalize → neighbor
# ============================================================================

print(f"\n  PHASE 3: Compute neighbors for each class...")
t2 = time.time()

# Build lookup: canonical_form -> class_id
# Already have canon_to_id

edges = set()
self_loops = set()
degrees = Counter()

for ci, (canon, bits) in enumerate(class_reps.items()):
    cid = canon_to_id[canon]
    adj = tournament_adj(n, bits)

    neighbors = set()

    for arc_pos in range(m):
        # Flip this arc
        flipped_bits = bits ^ (1 << arc_pos)
        flipped_adj = tournament_adj(n, flipped_bits)
        flipped_canon = canonical_form(flipped_adj, n)

        if flipped_canon in canon_to_id:
            nbr_id = canon_to_id[flipped_canon]
            if nbr_id == cid:
                self_loops.add(cid)
            else:
                neighbors.add(nbr_id)
        else:
            print(f"    WARNING: flipped canon not found for class {cid}, arc {arc_pos}")

    for nbr_id in neighbors:
        edge = (min(cid, nbr_id), max(cid, nbr_id))
        edges.add(edge)

    degrees[cid] = len(neighbors)

    if (ci + 1) % 50 == 0:
        elapsed = time.time() - t2
        print(f"    {ci+1}/{class_id} classes processed ({elapsed:.1f}s)")

print(f"  Phase 3 done in {time.time()-t2:.1f}s")

# ============================================================================
# RESULTS
# ============================================================================

print(f"\n{'='*80}")
print(f"  RESULTS: G_7")
print(f"{'='*80}")

print(f"  Vertices (iso classes): {class_id}")
print(f"  Edges: {len(edges)}")
print(f"  Self-loop classes: {len(self_loops)}")
print(f"  Degree sequence: min={min(degrees.values())}, max={max(degrees.values())}, "
      f"avg={sum(degrees.values())/len(degrees):.2f}")
print(f"  Sum of degrees: {sum(degrees.values())} (should = 2 × {len(edges)} = {2*len(edges)})")

# Degree distribution
deg_dist = Counter(degrees.values())
print(f"  Degree distribution: {dict(sorted(deg_dist.items()))}")

# The sequence
print(f"\n  THE SEQUENCE |E(G_n)| for n=3..7:")
print(f"    n=3: 1")
print(f"    n=4: 5")
print(f"    n=5: 30")
print(f"    n=6: 290")
print(f"    n=7: {len(edges)}")

total_time = time.time() - t0
print(f"\n  Total time: {total_time:.1f}s")

print(f"\n{'='*80}")
print(f"  DONE")
print(f"{'='*80}")
