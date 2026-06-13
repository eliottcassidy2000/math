#!/usr/bin/env python3
"""
gn_n8_nauty_s216.py — Compute |E(G_8)| using nauty for class enumeration
opus-2026-03-22-S216

nauty's gentourng generates all 6880 non-isomorphic tournaments on 8 vertices
in 0.004 seconds. This bypasses the need to enumerate 268M tournaments!

Strategy:
1. Read 6880 class representatives from nauty
2. Compute fast hash (score, c3, H) for each → build lookup table
3. For each class, flip 28 arcs, look up neighbor → build edges
4. For ambiguous hash groups, use canonical form

Also compute: SC status, |Aut|, complement pairing, merged graph.
"""

import sys
import subprocess
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  G_8 VIA NAUTY: |E(G_8)|, SC_8, MERGED GRAPH")
print("  opus-2026-03-22-S216")
print("=" * 80)

n = 8
m = comb(n, 2)  # 28
nfact = factorial(n)  # 40320

PAIRS = [(i, j) for i in range(n) for j in range(i+1, n)]

# ============================================================================
# PHASE 1: Read tournaments from nauty
# ============================================================================

print(f"\n  PHASE 1: Reading tournaments from gentourng n={n}...")
t0 = time.time()

result = subprocess.run(['gentourng', str(n)], capture_output=True, text=True)
lines = result.stdout.strip().split('\n') if result.stdout else result.stderr.strip().split('\n')

# Filter: tournament lines are binary strings of length m
tournaments = []
for line in lines:
    line = line.strip()
    if len(line) == m and all(c in '01' for c in line):
        bits = int(line, 2)
        tournaments.append(bits)

print(f"  Read {len(tournaments)} tournaments in {time.time()-t0:.3f}s")

# ============================================================================
# HELPERS
# ============================================================================

def bits_to_adj(bits):
    adj = [[0]*n for _ in range(n)]
    for k, (i, j) in enumerate(PAIRS):
        if bits & (1 << (m - 1 - k)):  # nauty: MSB first
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return adj

def adj_to_bits_nauty(adj):
    """Convert adj to nauty-style bits (MSB first)."""
    bits = 0
    for k, (i, j) in enumerate(PAIRS):
        if adj[i][j]:
            bits |= (1 << (m - 1 - k))
    return bits

def score_seq(adj):
    return tuple(sorted(sum(adj[i][j] for j in range(n)) for i in range(n)))

def c3_count(adj):
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if adj[i][j] and adj[j][k] and adj[k][i]: c3 += 1
                if adj[i][k] and adj[k][j] and adj[j][i]: c3 += 1
    return c3

def H_dp(adj):
    dp = {}
    for v in range(n): dp[(1 << v, v)] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            if not (S & (1 << v)): continue
            val = dp.get((S, v), 0)
            if val == 0: continue
            for u in range(n):
                if S & (1 << u): continue
                if adj[v][u]:
                    key = (S | (1 << u), u)
                    dp[key] = dp.get(key, 0) + val
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def canonical_form(adj):
    """Canonical form via degree-based pruning."""
    # Sort vertices by out-degree first
    scores = [sum(adj[i][j] for j in range(n)) for i in range(n)]
    # Group vertices by score
    score_groups = defaultdict(list)
    for v in range(n):
        score_groups[scores[v]].append(v)

    # Generate only permutations compatible with score ordering
    # For full generality, fall back to all perms
    # (optimization: if all scores distinct, only 1 perm to check)
    sorted_scores = sorted(set(scores))
    groups = [score_groups[s] for s in sorted_scores]

    if all(len(g) == 1 for g in groups):
        # All scores distinct: unique canonical ordering
        perm = [g[0] for g in groups]
        return tuple(adj[perm[i]][perm[j]] for i in range(n) for j in range(n))

    # Generate compatible permutations
    from itertools import product
    best = None
    def gen_perms_from_groups(groups):
        if not groups:
            yield []
            return
        for p in permutations(groups[0]):
            for rest in gen_perms_from_groups(groups[1:]):
                yield list(p) + rest

    for perm in gen_perms_from_groups(groups):
        form = tuple(adj[perm[i]][perm[j]] for i in range(n) for j in range(n))
        if best is None or form < best:
            best = form
    return best

def complement_adj(adj):
    return [[1 - adj[i][j] if i != j else 0 for j in range(n)] for i in range(n)]

# ============================================================================
# PHASE 2: Build class database with fast hashes
# ============================================================================

print(f"\n  PHASE 2: Computing invariants for {len(tournaments)} classes...")
t1 = time.time()

class_data = []  # list of dicts
hash_to_cids = defaultdict(list)  # (score, c3, H) -> [class_ids]
canon_to_cid = {}

for cid, bits in enumerate(tournaments):
    adj = bits_to_adj(bits)
    ss = score_seq(adj)
    c3 = c3_count(adj)
    h = H_dp(adj)
    canon = canonical_form(adj)

    class_data.append({
        'bits': bits, 'adj': adj, 'score': ss, 'c3': c3, 'H': h,
        'canon': canon, 'cid': cid
    })
    hash_to_cids[(ss, c3, h)].append(cid)
    canon_to_cid[canon] = cid

    if (cid + 1) % 500 == 0:
        elapsed = time.time() - t1
        print(f"    {cid+1}/{len(tournaments)} ({elapsed:.1f}s)")

print(f"  Phase 2 done in {time.time()-t1:.1f}s")
print(f"  Hash groups: {len(hash_to_cids)}")
ambiguous = sum(1 for v in hash_to_cids.values() if len(v) > 1)
print(f"  Ambiguous groups (>1 class): {ambiguous}")

# ============================================================================
# PHASE 3: Compute SC status, complement pairing
# ============================================================================

print(f"\n  PHASE 3: SC status and complement pairing...")
t2 = time.time()

sc_status = {}
complement_of = {}

for cid, data in enumerate(class_data):
    comp = complement_adj(data['adj'])
    comp_canon = canonical_form(comp)
    sc_status[cid] = (comp_canon == data['canon'])
    complement_of[cid] = canon_to_cid.get(comp_canon, -1)

sc_count = sum(1 for v in sc_status.values() if v)
print(f"  SC classes: {sc_count} (Burnside prediction: 176)")
print(f"  MATCH: {'YES ✓' if sc_count == 176 else 'NO ✗'}")

v_merged = (len(tournaments) + sc_count) // 2
print(f"  V_merged = ({len(tournaments)} + {sc_count}) / 2 = {v_merged} (prediction: 3528)")
print(f"  MATCH: {'YES ✓' if v_merged == 3528 else 'NO ✗'}")

# ============================================================================
# PHASE 4: Compute edges
# ============================================================================

print(f"\n  PHASE 4: Computing edges...")
t3 = time.time()

def lookup_class(adj):
    """Look up class ID for a tournament."""
    ss = score_seq(adj)
    c3 = c3_count(adj)
    h = H_dp(adj)
    hkey = (ss, c3, h)
    cids = hash_to_cids.get(hkey)
    if cids is None:
        return -1
    if len(cids) == 1:
        return cids[0]
    # Ambiguous: use canonical form
    canon = canonical_form(adj)
    return canon_to_cid.get(canon, -1)

edges = set()
self_loops = set()
degrees = Counter()

for cid in range(len(tournaments)):
    bits = class_data[cid]['bits']
    adj = class_data[cid]['adj']
    nbrs = set()

    for arc_pos in range(m):
        # Flip arc at position arc_pos
        flipped_bits = bits ^ (1 << (m - 1 - arc_pos))
        flipped_adj = bits_to_adj(flipped_bits)

        nb = lookup_class(flipped_adj)
        if nb == cid:
            self_loops.add(cid)
        elif nb >= 0:
            nbrs.add(nb)

    for nb in nbrs:
        edges.add((min(cid, nb), max(cid, nb)))
    degrees[cid] = len(nbrs)

    if (cid + 1) % 500 == 0:
        elapsed = time.time() - t3
        rate = (cid + 1) / elapsed
        eta = (len(tournaments) - cid - 1) / rate
        print(f"    {cid+1}/{len(tournaments)} ({elapsed:.1f}s, ETA {eta:.0f}s, "
              f"{len(edges)} edges)")

print(f"  Phase 4 done in {time.time()-t3:.1f}s")

# ============================================================================
# PHASE 5: Build merged graph
# ============================================================================

print(f"\n  PHASE 5: Building merged graph...")

merged_node = {}
mid = 0
for cid in range(len(tournaments)):
    if cid in merged_node:
        continue
    if sc_status[cid]:
        merged_node[cid] = mid
        mid += 1
    else:
        comp = complement_of[cid]
        merged_node[cid] = mid
        if comp != cid and comp >= 0:
            merged_node[comp] = mid
        mid += 1

merged_edges = set()
collapsed = 0
for (a, b) in edges:
    ma = merged_node[a]
    mb = merged_node[b]
    if ma == mb:
        collapsed += 1
    else:
        merged_edges.add((min(ma, mb), max(ma, mb)))

e_merged = len(merged_edges)
twin = len(edges) - e_merged - collapsed

# ============================================================================
# RESULTS
# ============================================================================

print(f"\n{'='*80}")
print(f"  RESULTS: G_8")
print(f"{'='*80}")
print(f"  Vertices: {len(tournaments)}")
print(f"  EDGES: {len(edges)}")
print(f"  Self-loops: {len(self_loops)}")
print(f"  avg_deg: {sum(degrees.values())/len(tournaments):.2f}")
print(f"  SC classes: {sc_count}")
print(f"  V_merged: {mid}")
print(f"  E_merged: {e_merged}")
print(f"  Collapsed: {collapsed}")
print(f"  Twin: {twin}")

# Burnside predictions
T8 = 188288
print(f"\n  Burnside predictions:")
print(f"    T_8 = {T8}")
print(f"    T_8 / E = {T8/len(edges):.4f}")
print(f"    T_8 / (2*E) = {T8/(2*len(edges)):.4f}")

# The complete sequence
print(f"\n  *** THE SEQUENCE |E(G_n)| ***")
print(f"  n=3: |V|=    2,  |E|=      1")
print(f"  n=4: |V|=    4,  |E|=      5")
print(f"  n=5: |V|=   12,  |E|=     30")
print(f"  n=6: |V|=   56,  |E|=    290")
print(f"  n=7: |V|=  456,  |E|=   4086")
print(f"  n=8: |V|= 6880,  |E|= {len(edges)}")

# Edge decomposition
print(f"\n  EDGE DECOMPOSITION:")
print(f"  {'n':>3} {'E_orig':>8} {'E_merged':>8} {'collapsed':>9} {'twin':>6}")
for nn, eo, em, co, tw in [
    (3, 1, 1, 0, 0), (4, 5, 3, 0, 2), (5, 30, 21, 0, 9),
    (6, 290, 143, 5, 142), (7, 4086, 2123, 0, 1963),
    (8, len(edges), e_merged, collapsed, twin)]:
    print(f"  {nn:3d} {eo:8d} {em:8d} {co:9d} {tw:6d}")

total_time = time.time() - t0
print(f"\n  Total time: {total_time:.1f}s")
print(f"{'='*80}")
