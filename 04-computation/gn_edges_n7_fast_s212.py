#!/usr/bin/env python3
"""
gn_edges_n7_fast_s212.py — Compute |E(G_7)| efficiently
opus-2026-03-22-S212

Optimized approach:
1. Hash all 2^21 tournaments by (score, c3, sorted deletion-H fingerprint)
   - deletion-H fingerprint: sorted tuple of H(T-v) for all v
   - This is O(7 × 2^6 × 6) = O(2700) per tournament (MUCH cheaper than canonical)
2. Each hash group is almost certainly a single iso class
3. For any ambiguous groups, use canonical form to split
4. For each class representative, flip 21 arcs, hash the result → neighbor class
5. Count edges

Expected: ~5 minutes total.
"""

import sys
from math import comb
from itertools import permutations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  COMPUTING |E(G_7)| EFFICIENTLY")
print("  opus-2026-03-22-S212")
print("=" * 80)

n = 7
m = comb(n, 2)  # 21
total = 1 << m   # 2097152

# Precompute pair indices
PAIRS = []
idx = 0
for i in range(n):
    for j in range(i+1, n):
        PAIRS.append((i, j))
        idx += 1

print(f"  n={n}, m={m}, |tournaments|={total}")

def bits_to_adj(bits):
    adj = [[0]*n for _ in range(n)]
    for k, (i, j) in enumerate(PAIRS):
        if bits & (1 << k):
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return adj

def H_dp_n(adj, nn):
    """Count Hamiltonian paths for a tournament of size nn."""
    dp = {}
    for v in range(nn):
        dp[(1 << v, v)] = 1
    for S in range(1, 1 << nn):
        for v in range(nn):
            if not (S & (1 << v)):
                continue
            val = dp.get((S, v), 0)
            if val == 0:
                continue
            for u in range(nn):
                if S & (1 << u):
                    continue
                if adj[v][u]:
                    key = (S | (1 << u), u)
                    dp[key] = dp.get(key, 0) + val
    full = (1 << nn) - 1
    return sum(dp.get((full, v), 0) for v in range(nn))

def score_seq(adj):
    return tuple(sorted(sum(adj[i][j] for j in range(n)) for i in range(n)))

def c3_count(adj):
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if adj[i][j] and adj[j][k] and adj[k][i]:
                    c3 += 1
                if adj[i][k] and adj[k][j] and adj[j][i]:
                    c3 += 1
    return c3

def deletion_fingerprint(adj):
    """Compute sorted tuple of H(T-v) for all v in [n].
    Invariant under isomorphism. Very discriminating."""
    h_vals = []
    for v in range(n):
        # Build T-v: tournament on n-1 vertices (excluding v)
        verts = [i for i in range(n) if i != v]
        nn = len(verts)
        sub_adj = [[adj[verts[i]][verts[j]] for j in range(nn)] for i in range(nn)]
        h = H_dp_n(sub_adj, nn)
        h_vals.append(h)
    return tuple(sorted(h_vals))

def canonical_form(adj):
    """Canonical form (expensive: n! permutations)."""
    best = None
    for perm in permutations(range(n)):
        form = tuple(adj[perm[i]][perm[j]] for i in range(n) for j in range(n))
        if best is None or form < best:
            best = form
    return best

# ============================================================================
# PHASE 1: Hash all tournaments by (score, c3, deletion_fingerprint)
# ============================================================================

print(f"\n  PHASE 1: Computing fingerprints for all {total} tournaments...")
t0 = time.time()

# Store: hash -> [bits]
hash_groups = defaultdict(list)

for bits in range(total):
    adj = bits_to_adj(bits)
    ss = score_seq(adj)
    c3 = c3_count(adj)
    dfp = deletion_fingerprint(adj)

    hkey = (ss, c3, dfp)
    hash_groups[hkey].append(bits)

    if (bits + 1) % 50000 == 0:
        elapsed = time.time() - t0
        rate = (bits + 1) / elapsed
        eta = (total - bits - 1) / rate
        print(f"    {bits+1}/{total} ({rate:.0f}/s, ETA {eta:.0f}s, "
              f"{len(hash_groups)} groups so far)")

print(f"  Phase 1 done: {len(hash_groups)} hash groups in {time.time()-t0:.1f}s")

# Verify
total_members = sum(len(v) for v in hash_groups.values())
print(f"  Total tournaments: {total_members} (should be {total})")

group_sizes = Counter(len(v) for v in hash_groups.values())
print(f"  Group size distribution: {dict(sorted(group_sizes.items()))}")

# Check: if all groups have sizes that are multiples of n!/gcd,
# they're likely single iso classes
nfact = 5040  # 7!
ambiguous = 0
for hkey, members in hash_groups.items():
    sz = len(members)
    if nfact % (nfact // 1) != 0:  # always true
        pass
    # Check if size is n!/k for some integer k (= |Aut|)
    if nfact % sz != 0 and sz % nfact != 0:
        ambiguous += 1

# Better check: expected class sizes divide n!
# A class with |Aut|=a has size n!/a
# Valid sizes: 5040, 2520, 1680, 1260, 1008, 840, 720, 630, 560, 504, ...
# Groups whose size isn't a valid class size must contain multiple classes

valid_sizes = set()
for a in range(1, nfact + 1):
    if nfact % a == 0:
        valid_sizes.add(nfact // a)

needs_splitting = []
single_class_groups = 0
for hkey, members in hash_groups.items():
    sz = len(members)
    if sz in valid_sizes:
        single_class_groups += 1
    else:
        needs_splitting.append((hkey, sz))

print(f"  Single-class groups (size divides n!): {single_class_groups}")
print(f"  Groups needing further splitting: {len(needs_splitting)}")
if needs_splitting:
    for hkey, sz in needs_splitting[:10]:
        # Find possible decompositions
        print(f"    size={sz}: possible splits = ", end="")
        # sz = sum of valid_sizes
        decomps = []
        for vs in sorted(valid_sizes, reverse=True):
            if sz % vs == 0:
                decomps.append(f"{sz//vs}×{vs}")
        print(", ".join(decomps[:5]))

# ============================================================================
# PHASE 2: Assign class IDs
# ============================================================================

print(f"\n  PHASE 2: Assigning class IDs...")
t1 = time.time()

class_id_map = {}  # hash_key -> class_id (for single-class groups)
                    # or: canonical_form -> class_id (for split groups)
class_reps = {}    # class_id -> bits
class_sizes_map = {}  # class_id -> size
cid = 0

hash_to_cid = {}  # hash_key -> class_id (for fast lookup later)

for hkey, members in hash_groups.items():
    sz = len(members)
    if sz in valid_sizes:
        # Single class
        hash_to_cid[hkey] = cid
        class_reps[cid] = members[0]
        class_sizes_map[cid] = sz
        cid += 1
    else:
        # Need to split using canonical form
        sub_canons = {}
        for bits in members:
            adj = bits_to_adj(bits)
            canon = canonical_form(adj)
            if canon not in sub_canons:
                sub_canons[canon] = cid
                class_reps[cid] = bits
                class_sizes_map[cid] = 0
                cid += 1
            class_sizes_map[sub_canons[canon]] += 1

        # Store mapping for these members
        # (We'll need to look up each member individually for edge computation)
        # Store the sub_canons for this hash group
        hash_to_cid[hkey] = sub_canons  # dict, not int

print(f"  Phase 2 done: {cid} iso classes in {time.time()-t1:.1f}s")
print(f"  (Expected 456 from A000568)")

# ============================================================================
# PHASE 3: Compute edges
# ============================================================================

print(f"\n  PHASE 3: Computing edges for {cid} classes...")
t2 = time.time()

# For each class representative, flip 21 arcs, find neighbor class
edges = set()
self_loops = set()
degree_counter = Counter()

def lookup_class(bits):
    """Look up the class ID for a given tournament (bits)."""
    adj = bits_to_adj(bits)
    ss = score_seq(adj)
    c3 = c3_count(adj)
    dfp = deletion_fingerprint(adj)
    hkey = (ss, c3, dfp)

    entry = hash_to_cid.get(hkey)
    if entry is None:
        return -1  # not found
    if isinstance(entry, int):
        return entry  # single class
    else:
        # Need canonical form to disambiguate
        canon = canonical_form(adj)
        return entry.get(canon, -1)

for ci in range(cid):
    bits = class_reps[ci]
    neighbors = set()

    for arc_pos in range(m):
        flipped = bits ^ (1 << arc_pos)
        nbr_cid = lookup_class(flipped)

        if nbr_cid == -1:
            print(f"    WARNING: class {ci} arc {arc_pos} → unknown class!")
            continue

        if nbr_cid == ci:
            self_loops.add(ci)
        else:
            neighbors.add(nbr_cid)

    for nbr in neighbors:
        edge = (min(ci, nbr), max(ci, nbr))
        edges.add(edge)

    degree_counter[ci] = len(neighbors)

    if (ci + 1) % 50 == 0:
        elapsed = time.time() - t2
        rate = (ci + 1) / elapsed
        eta = (cid - ci - 1) / rate
        print(f"    {ci+1}/{cid} classes ({rate:.1f}/s, ETA {eta:.0f}s, "
              f"{len(edges)} edges so far)")

print(f"  Phase 3 done in {time.time()-t2:.1f}s")

# ============================================================================
# RESULTS
# ============================================================================

print(f"\n{'='*80}")
print(f"  RESULTS: G_7")
print(f"{'='*80}")

print(f"  Vertices: {cid}")
print(f"  EDGES: {len(edges)}")
print(f"  Self-loop classes: {len(self_loops)}/{cid}")
print(f"  Sum of degrees: {sum(degree_counter.values())} (should = 2 × {len(edges)})")
print(f"  Degree: min={min(degree_counter.values())}, max={max(degree_counter.values())}, "
      f"avg={sum(degree_counter.values())/cid:.2f}")

# Degree distribution
deg_dist = Counter(degree_counter.values())
print(f"  Degree distribution: {dict(sorted(deg_dist.items()))}")

# The complete sequence
print(f"\n  *** THE SEQUENCE |E(G_n)| ***")
print(f"  n=3: |V|=2,   |E|=1")
print(f"  n=4: |V|=4,   |E|=5")
print(f"  n=5: |V|=12,  |E|=30")
print(f"  n=6: |V|=56,  |E|=290")
print(f"  n=7: |V|={cid}, |E|={len(edges)}")

# Key ratios
avg_deg = sum(degree_counter.values()) / cid
print(f"\n  avg_deg = {avg_deg:.4f}")
print(f"  avg_deg / m = {avg_deg/m:.4f}")
print(f"  |E| / |V| = {len(edges)/cid:.4f}")
print(f"  density = {2*len(edges)/(cid*(cid-1)):.6f}")

total_time = time.time() - t0
print(f"\n  Total time: {total_time:.1f}s")
print(f"{'='*80}")
