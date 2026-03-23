#!/usr/bin/env python3
"""
overhead_n7_s20cn.py -- kind-pasteur-2026-03-22-S20cn

Compute EXACT overhead decomposition T_7 - 2|E_7| = SL + MW at n=7.

Uses the fast fingerprint approach from gn_edges_n7_fast_s212.py,
then additionally computes:
  - Aut(T) for each iso class representative
  - Arc orbits under Aut
  - Self-loop vs cross orbit classification
  - Multi-weight edge detection

Expected: ~30 minutes total (Phase 1 dominates).

Author: kind-pasteur-2026-03-22-S20cn
"""
import sys
import time
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter

sys.stdout.reconfigure(line_buffering=True)

n = 7
m = comb(n, 2)  # 21
total = 1 << m   # 2097152
nfact = factorial(n)  # 5040

PAIRS = []
for i in range(n):
    for j in range(i+1, n):
        PAIRS.append((i, j))

# Precompute pair lookup
PAIR_INDEX = {}
for k, (i, j) in enumerate(PAIRS):
    PAIR_INDEX[(i, j)] = k
    PAIR_INDEX[(j, i)] = k

print("=" * 70)
print("  OVERHEAD DECOMPOSITION AT n=7: T-2E = SL + MW")
print("  kind-pasteur-2026-03-22-S20cn")
print("=" * 70)

def bits_to_adj(bits):
    adj = [[0]*n for _ in range(n)]
    for k, (i, j) in enumerate(PAIRS):
        if bits & (1 << k):
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return adj

def H_dp(adj, nn):
    """Count Hamiltonian paths via DP."""
    dp = [0] * ((1 << nn) * nn)
    for v in range(nn):
        dp[(1 << v) * nn + v] = 1
    for S in range(1, 1 << nn):
        for v in range(nn):
            if not (S & (1 << v)):
                continue
            val = dp[S * nn + v]
            if val == 0:
                continue
            for u in range(nn):
                if S & (1 << u):
                    continue
                if adj[v][u]:
                    dp[(S | (1 << u)) * nn + u] += val
    full = (1 << nn) - 1
    return sum(dp[full * nn + v] for v in range(nn))

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
    h_vals = []
    for v in range(n):
        verts = [i for i in range(n) if i != v]
        nn = len(verts)
        sub_adj = [[adj[verts[i]][verts[j]] for j in range(nn)] for i in range(nn)]
        h = H_dp(sub_adj, nn)
        h_vals.append(h)
    return tuple(sorted(h_vals))

def canonical_form(adj):
    best = None
    for perm in permutations(range(n)):
        form = tuple(adj[perm[i]][perm[j]] for i in range(n) for j in range(n))
        if best is None or form < best:
            best = form
    return best

# ============================================================================
# PHASE 1: Hash all tournaments
# ============================================================================
print(f"\n  PHASE 1: Computing fingerprints for {total} tournaments...")
t0 = time.time()

hash_groups = defaultdict(list)

for bits in range(total):
    adj = bits_to_adj(bits)
    ss = score_seq(adj)
    c3 = c3_count(adj)
    dfp = deletion_fingerprint(adj)
    hkey = (ss, c3, dfp)
    hash_groups[hkey].append(bits)

    if (bits + 1) % 100000 == 0:
        elapsed = time.time() - t0
        rate = (bits + 1) / elapsed
        eta = (total - bits - 1) / rate
        print(f"    {bits+1}/{total} ({rate:.0f}/s, ETA {eta:.0f}s, {len(hash_groups)} groups)")

print(f"  Phase 1 done: {len(hash_groups)} hash groups in {time.time()-t0:.1f}s")

# ============================================================================
# PHASE 2: Assign class IDs, build bits_to_class
# ============================================================================
print(f"\n  PHASE 2: Identifying classes...")
t1 = time.time()

valid_sizes = set()
for a in range(1, nfact + 1):
    if nfact % a == 0:
        valid_sizes.add(nfact // a)

class_reps = {}
class_sizes = {}
class_auts = {}
bits_to_class = {}
cid = 0

for hkey, members in hash_groups.items():
    sz = len(members)
    if sz in valid_sizes:
        class_reps[cid] = members[0]
        class_sizes[cid] = sz
        class_auts[cid] = nfact // sz
        for b in members:
            bits_to_class[b] = cid
        cid += 1
    else:
        # Split using canonical form
        sub_canons = {}
        for b in members:
            adj = bits_to_adj(b)
            canon = canonical_form(adj)
            if canon not in sub_canons:
                sub_canons[canon] = cid
                class_reps[cid] = b
                class_sizes[cid] = 0
                cid += 1
            class_sizes[sub_canons[canon]] += 1
            bits_to_class[b] = sub_canons[canon]
        for c in sub_canons.values():
            class_auts[c] = nfact // class_sizes[c]

N_classes = cid
print(f"  Phase 2 done: {N_classes} iso classes in {time.time()-t1:.1f}s")
print(f"  (Expected 456)")

# ============================================================================
# PHASE 3: Compute Aut, arc orbits, SL/MW for each class
# ============================================================================
print(f"\n  PHASE 3: Computing arc orbit decomposition for {N_classes} classes...")
t2 = time.time()

total_SL = 0
total_MW = 0
total_degree = 0
total_arc_orbits = 0

SL_by_aut = Counter()
MW_by_aut = Counter()
degree_list = []

for ci in range(N_classes):
    bits = class_reps[ci]
    adj = bits_to_adj(bits)
    aut_size = class_auts[ci]

    # Find Aut(T) by checking all n! permutations
    aut_perms = []
    for perm in permutations(range(n)):
        is_aut = True
        for i in range(n):
            for j in range(n):
                if adj[perm[i]][perm[j]] != adj[i][j]:
                    is_aut = False
                    break
            if not is_aut:
                break
        if is_aut:
            aut_perms.append(perm)

    # Compute arc orbits under Aut
    arc_orbit_id = [-1] * m
    orbit_count = 0
    for k in range(m):
        if arc_orbit_id[k] >= 0:
            continue
        orbit = set()
        for perm in aut_perms:
            i, j = PAIRS[k]
            pi, pj = perm[i], perm[j]
            if pi < pj:
                new_k = PAIR_INDEX[(pi, pj)]
            else:
                new_k = PAIR_INDEX[(pj, pi)]
            orbit.add(new_k)
        for a in orbit:
            arc_orbit_id[a] = orbit_count
        orbit_count += 1

    n_orbits = orbit_count

    # For each orbit: pick representative arc, flip, check neighbor class
    orbit_rep = {}
    for k in range(m):
        oid = arc_orbit_id[k]
        if oid not in orbit_rep:
            orbit_rep[oid] = k

    orbit_targets = {}
    for oid, rep_arc in orbit_rep.items():
        flipped_bits = bits ^ (1 << rep_arc)
        nb_class = bits_to_class[flipped_bits]
        orbit_targets[oid] = nb_class

    # Count SL, degree, MW
    self_loop_orbits = sum(1 for oid in range(orbit_count) if orbit_targets[oid] == ci)
    cross_targets = [orbit_targets[oid] for oid in range(orbit_count) if orbit_targets[oid] != ci]
    neighbor_classes = set(cross_targets)
    degree = len(neighbor_classes)
    cross_orbits = len(cross_targets)
    multi_weight = cross_orbits - degree

    total_SL += self_loop_orbits
    total_MW += multi_weight
    total_degree += degree
    total_arc_orbits += n_orbits

    SL_by_aut[aut_size] += self_loop_orbits
    MW_by_aut[aut_size] += multi_weight
    degree_list.append(degree)

    if (ci + 1) % 100 == 0:
        elapsed = time.time() - t2
        rate = (ci + 1) / elapsed
        eta = (N_classes - ci - 1) / rate
        print(f"    {ci+1}/{N_classes} ({rate:.1f}/s, ETA {eta:.0f}s)")

print(f"  Phase 3 done in {time.time()-t2:.1f}s")

# ============================================================================
# RESULTS
# ============================================================================
edges = total_degree // 2
overhead = total_SL + total_MW
T_n = total_arc_orbits

print(f"\n{'='*70}")
print(f"  RESULTS: n = 7")
print(f"{'='*70}")
print(f"  N = {N_classes} iso classes")
print(f"  T_n = {T_n} transition orbits (expected 8912)")
print(f"  |E| = {edges} edges (expected 4086)")
print(f"  2|E| = {total_degree}")
print(f"  overhead = T - 2|E| = {overhead} (expected 740)")
print(f"    SL (self-loop orbits) = {total_SL}")
print(f"    MW (multi-weight)     = {total_MW}")
print(f"    SL + MW = {total_SL + total_MW}")
print(f"  Check: T = SL + 2E + MW = {total_SL + total_degree + total_MW} (should be {T_n})")

print(f"\n  Breakdown by |Aut|:")
all_auts = sorted(set(class_auts.values()))
for a in all_auts:
    count = sum(1 for ci in range(N_classes) if class_auts[ci] == a)
    sl = SL_by_aut.get(a, 0)
    mw = MW_by_aut.get(a, 0)
    print(f"    |Aut|={a}: {count} classes, SL={sl}, MW={mw}")

print(f"\n  Self-loop analysis:")
print(f"    SL_7 / V_7 = {total_SL / N_classes:.4f}")
print(f"    MW_7 / V_7 = {total_MW / N_classes:.4f}")
print(f"    SL_7 / T_7 = {total_SL / T_n:.4f}")
print(f"    MW_7 / T_7 = {total_MW / T_n:.4f}")

# Degree distribution
deg_dist = Counter(degree_list)
print(f"\n  Degree distribution:")
for d in sorted(deg_dist.keys()):
    print(f"    deg={d:>3d}: {deg_dist[d]} classes")

# ============================================================================
# COMPLETE SUMMARY TABLE (including n=3..6 from prior computation)
# ============================================================================
print(f"\n{'='*70}")
print(f"  COMPLETE OVERHEAD DECOMPOSITION TABLE")
print(f"{'='*70}")
print(f"  {'n':>3s} {'V':>6s} {'m':>3s} {'T':>7s} {'|E|':>6s} {'OH':>5s} {'SL':>5s} {'MW':>5s}")
print(f"  {3:>3d} {2:>6d} {3:>3d} {4:>7d} {1:>6d} {2:>5d} {2:>5d} {0:>5d}")
print(f"  {4:>3d} {4:>6d} {6:>3d} {16:>7d} {5:>6d} {6:>5d} {6:>5d} {0:>5d}")
print(f"  {5:>3d} {12:>6d} {10:>3d} {88:>7d} {30:>6d} {28:>5d} {16:>5d} {12:>5d}")
print(f"  {6:>3d} {56:>6d} {15:>3d} {704:>7d} {290:>6d} {124:>5d} {58:>5d} {66:>5d}")
print(f"  {7:>3d} {N_classes:>6d} {m:>3d} {T_n:>7d} {edges:>6d} {overhead:>5d} {total_SL:>5d} {total_MW:>5d}")

total_time = time.time() - t0
print(f"\n  Total time: {total_time:.1f}s ({total_time/60:.1f} minutes)")
print(f"{'='*70}")
