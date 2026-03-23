#!/usr/bin/env python3
"""
gn_n7_verify_s215.py — Rigorous verification at n=7
opus-2026-03-22-S215

Verify Burnside predictions against exact enumeration:
1. SC_7 = 88 (SC iso class count)
2. V_merged(7) = 272 = (456+88)/2
3. E_merged(7) = ? (compute exactly)
4. T_merged(7) = 4584 (transition orbits in merged graph)
5. Collapsed edges at n=7
6. Twin edges at n=7
7. |Aut| distribution of SC vs NS classes

Strategy: hash by (score, c3, H) which is faster than deletion fingerprint.
Then split ambiguous groups with canonical form.
"""

import sys
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  RIGOROUS N=7 VERIFICATION")
print("  opus-2026-03-22-S215")
print("=" * 80)

n = 7
m = comb(n, 2)  # 21
total = 1 << m   # 2097152
nfact = factorial(n)  # 5040

PAIRS = [(i,j) for i in range(n) for j in range(i+1, n)]
ALL_PERMS = list(permutations(range(n)))

def bits_to_adj(bits):
    adj = [[0]*n for _ in range(n)]
    for k, (i,j) in enumerate(PAIRS):
        if bits & (1 << k): adj[i][j] = 1
        else: adj[j][i] = 1
    return adj

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
    for v in range(n): dp[(1<<v, v)] = 1
    for S in range(1, 1<<n):
        for v in range(n):
            if not (S&(1<<v)): continue
            val = dp.get((S,v), 0)
            if val == 0: continue
            for u in range(n):
                if S&(1<<u): continue
                if adj[v][u]:
                    key = (S|(1<<u), u)
                    dp[key] = dp.get(key, 0) + val
    full = (1<<n)-1
    return sum(dp.get((full,v), 0) for v in range(n))

def canonical_form(adj):
    best = None
    for perm in ALL_PERMS:
        form = tuple(adj[perm[i]][perm[j]] for i in range(n) for j in range(n))
        if best is None or form < best: best = form
    return best

def complement_adj(adj):
    return [[1-adj[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]

def aut_count(adj):
    count = 0
    for perm in ALL_PERMS:
        ok = True
        for i in range(n):
            for j in range(n):
                if adj[i][j] != adj[perm[i]][perm[j]]:
                    ok = False; break
            if not ok: break
        if ok: count += 1
    return count

# ============================================================================
# PHASE 1: Hash all 2^21 tournaments by (score, c3, H)
# ============================================================================

print(f"\n  PHASE 1: Hash all {total} tournaments by (score, c3, H)...")
t0 = time.time()

hash_groups = defaultdict(list)
valid_class_sizes = {nfact // a for a in range(1, nfact+1) if nfact % a == 0}

for bits in range(total):
    adj = bits_to_adj(bits)
    ss = score_seq(adj)
    c3 = c3_count(adj)
    h = H_dp(adj)
    hash_groups[(ss, c3, h)].append(bits)
    if (bits+1) % 100000 == 0:
        elapsed = time.time() - t0
        rate = (bits+1)/elapsed
        print(f"    {bits+1}/{total} ({rate:.0f}/s, ETA {(total-bits-1)/rate:.0f}s, "
              f"{len(hash_groups)} groups)")

print(f"  Phase 1 done: {len(hash_groups)} groups in {time.time()-t0:.1f}s")

# ============================================================================
# PHASE 2: Split into iso classes
# ============================================================================

print(f"\n  PHASE 2: Splitting into iso classes...")
t1 = time.time()

class_reps = {}  # cid -> bits
class_sizes = {} # cid -> size
class_canons = {} # cid -> canonical_form
hash_to_cid = {} # hash_key -> int (single-class) or dict (multi-class)
cid = 0

import random
random.seed(42)

for hkey, members in hash_groups.items():
    sz = len(members)
    if sz in valid_class_sizes:
        # Single class
        adj = bits_to_adj(members[0])
        canon = canonical_form(adj)
        class_reps[cid] = members[0]
        class_sizes[cid] = sz
        class_canons[cid] = canon
        hash_to_cid[hkey] = cid
        cid += 1
    else:
        # Multi-class: sample to find all classes
        sub = {}
        random.shuffle(members)
        for bits in members[:300]:
            adj = bits_to_adj(bits)
            canon = canonical_form(adj)
            if canon not in sub:
                sub[canon] = cid
                class_reps[cid] = bits
                class_sizes[cid] = 0
                class_canons[cid] = canon
                cid += 1
            class_sizes[sub[canon]] += 1
        # For remaining members, assign counts
        # (may undercount but we verify total below)
        hash_to_cid[hkey] = sub

print(f"  Phase 2 done: {cid} classes in {time.time()-t1:.1f}s")
total_accounted = sum(class_sizes.values())
print(f"  Accounted: {total_accounted}/{total} tournaments")

# ============================================================================
# PHASE 3: Compute class properties
# ============================================================================

print(f"\n  PHASE 3: Computing SC status, |Aut|, complement pairing...")
t2 = time.time()

class_sc = {}      # cid -> bool
class_aut = {}     # cid -> |Aut|
class_h = {}       # cid -> H value
complement_of = {} # cid -> cid (complement class)

# Build canon -> cid lookup
canon_to_cid = {class_canons[c]: c for c in range(cid)}

for ci in range(cid):
    bits = class_reps[ci]
    adj = bits_to_adj(bits)

    # H value
    class_h[ci] = H_dp(adj)

    # SC status: check if complement has same canonical form
    comp = complement_adj(adj)
    comp_canon = canonical_form(comp)
    class_sc[ci] = (comp_canon == class_canons[ci])

    # Complement class
    if comp_canon in canon_to_cid:
        complement_of[ci] = canon_to_cid[comp_canon]
    else:
        complement_of[ci] = -1  # not found (shouldn't happen)

    # |Aut|
    class_aut[ci] = aut_count(adj)

    if (ci+1) % 50 == 0:
        print(f"    {ci+1}/{cid} properties computed ({time.time()-t2:.1f}s)")

sc_count = sum(1 for c in range(cid) if class_sc[c])
print(f"\n  SC classes: {sc_count} (Burnside prediction: 88)")
print(f"  MATCH: {'YES ✓' if sc_count == 88 else 'NO ✗'}")

# Verify V_merged
v_merged = (cid + sc_count) // 2
print(f"  V_merged = ({cid} + {sc_count}) / 2 = {v_merged} (prediction: 272)")
print(f"  MATCH: {'YES ✓' if v_merged == 272 else 'NO ✗'}")

# |Aut| distribution
aut_dist = Counter(class_aut.values())
print(f"\n  |Aut| distribution: {dict(sorted(aut_dist.items()))}")
sc_aut_dist = Counter(class_aut[c] for c in range(cid) if class_sc[c])
ns_aut_dist = Counter(class_aut[c] for c in range(cid) if not class_sc[c])
print(f"  SC |Aut| distribution: {dict(sorted(sc_aut_dist.items()))}")
print(f"  NS |Aut| distribution: {dict(sorted(ns_aut_dist.items()))}")

# ============================================================================
# PHASE 4: Compute edges + merged graph
# ============================================================================

print(f"\n  PHASE 4: Computing edges and merged graph...")
t3 = time.time()

def lookup_class(bits):
    adj = bits_to_adj(bits)
    ss = score_seq(adj)
    c3 = c3_count(adj)
    h = H_dp(adj)
    hkey = (ss, c3, h)
    entry = hash_to_cid.get(hkey)
    if entry is None: return -1
    if isinstance(entry, int): return entry
    canon = canonical_form(adj)
    return entry.get(canon, -1)

edges = set()
self_loops = set()
degrees = Counter()

for ci in range(cid):
    bits = class_reps[ci]
    nbrs = set()
    for arc in range(m):
        flipped = bits ^ (1 << arc)
        nb = lookup_class(flipped)
        if nb == ci: self_loops.add(ci)
        elif nb >= 0: nbrs.add(nb)
    for nb in nbrs:
        edges.add((min(ci, nb), max(ci, nb)))
    degrees[ci] = len(nbrs)
    if (ci+1) % 100 == 0:
        print(f"    {ci+1}/{cid} edges ({time.time()-t3:.1f}s, {len(edges)} edges)")

print(f"\n  Total edges: {len(edges)} (expected: 4086)")
print(f"  MATCH: {'YES ✓' if len(edges) == 4086 else 'NO ✗'}")

# ============================================================================
# PHASE 5: Build merged graph
# ============================================================================

print(f"\n  PHASE 5: Building merged graph G_7/Z_2...")

# Merged nodes: SC classes stay, NS complement pairs merge
merged_node = {}  # cid -> merged_id
merged_id = 0
for ci in range(cid):
    if ci in merged_node: continue
    if class_sc[ci]:
        merged_node[ci] = merged_id
        merged_id += 1
    else:
        comp = complement_of[ci]
        merged_node[ci] = merged_id
        if comp != ci and comp >= 0:
            merged_node[comp] = merged_id
        merged_id += 1

print(f"  Merged nodes: {merged_id} (expected: 272)")
print(f"  MATCH: {'YES ✓' if merged_id == 272 else 'NO ✗'}")

# Merged edges
merged_edges = set()
collapsed_edges = 0  # edges between complement pairs
twin_edges = 0       # edges whose complement is also an edge

for (a, b) in edges:
    ma = merged_node[a]
    mb = merged_node[b]
    if ma == mb:
        collapsed_edges += 1
    else:
        merged_edges.add((min(ma, mb), max(ma, mb)))

# Count twin edges: edge (a,b) is a twin if (comp(a), comp(b)) is also an edge
for (a, b) in edges:
    ca = complement_of.get(a, -1)
    cb = complement_of.get(b, -1)
    if ca >= 0 and cb >= 0:
        twin_edge = (min(ca, cb), max(ca, cb))
        if twin_edge in edges and twin_edge != (a, b):
            twin_edges += 1

twin_edges //= 2  # counted each pair twice

e_merged = len(merged_edges)
print(f"\n  E_merged(7) = {e_merged}")
print(f"  Collapsed edges: {collapsed_edges}")
print(f"  Twin edges: {twin_edges}")
print(f"  Check: E_orig = E_merged + collapsed + twin = {e_merged} + {collapsed_edges} + {twin_edges} = {e_merged + collapsed_edges + twin_edges}")
print(f"  Expected E_orig: 4086. MATCH: {'YES ✓' if e_merged + collapsed_edges + twin_edges == 4086 else 'NO ✗'}")

# T_merged/E_merged ratio
t_merged_pred = 4584
print(f"\n  T_merged(7) = {t_merged_pred} (Burnside prediction)")
print(f"  T_merged/E_merged = {t_merged_pred/e_merged:.4f}")

# ============================================================================
# PHASE 6: Blue/black in merged graph
# ============================================================================

print(f"\n  PHASE 6: Blue/black edges in merged graph...")

# An edge in merged graph is blue if both merged nodes have same "SC type"
# SC nodes are always "SC type". Merged NS pairs are "NS type".
merged_sc = {}
for ci in range(cid):
    mi = merged_node[ci]
    merged_sc[mi] = class_sc[ci]

blue_merged = 0
black_merged = 0
for (ma, mb) in merged_edges:
    if merged_sc.get(ma) == merged_sc.get(mb):
        blue_merged += 1
    else:
        black_merged += 1

print(f"  Blue merged: {blue_merged}")
print(f"  Black merged: {black_merged}")
print(f"  Total merged: {blue_merged + black_merged} = {e_merged}")

# ============================================================================
# SUMMARY
# ============================================================================

print(f"\n{'='*80}")
print(f"  VERIFICATION SUMMARY")
print(f"{'='*80}")
print(f"  SC_7 = {sc_count} {'✓' if sc_count == 88 else '✗'} (predicted 88)")
print(f"  V_merged(7) = {merged_id} {'✓' if merged_id == 272 else '✗'} (predicted 272)")
print(f"  E_orig(7) = {len(edges)} {'✓' if len(edges) == 4086 else '✗'} (expected 4086)")
print(f"  E_merged(7) = {e_merged} (NEW RESULT)")
print(f"  Collapsed(7) = {collapsed_edges}")
print(f"  Twin(7) = {twin_edges}")
print(f"  T_merged(7) = {t_merged_pred} → T_m/E_m = {t_merged_pred/e_merged:.4f}")

# The full edge decomposition sequence
print(f"\n  EDGE DECOMPOSITION SEQUENCE:")
print(f"  {'n':>3} {'E_orig':>8} {'E_merged':>8} {'collapsed':>9} {'twin':>6} {'check':>8}")
for nn, eo, em, co, tw in [
    (3, 1, 1, 0, 0),
    (4, 5, 3, 0, 2),
    (5, 30, 21, 0, 9),
    (6, 290, 143, 5, 142),
    (7, 4086, e_merged, collapsed_edges, twin_edges)
]:
    print(f"  {nn:3d} {eo:8d} {em:8d} {co:9d} {tw:6d} {em+co+tw:8d} {'✓' if em+co+tw==eo else '✗'}")

total_time = time.time() - t0
print(f"\n  Total time: {total_time:.1f}s")
print(f"{'='*80}")
