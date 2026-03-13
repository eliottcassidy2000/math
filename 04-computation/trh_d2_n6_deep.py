#!/usr/bin/env python3
"""
trh_d2_n6_deep.py — Find the combinatorial condition for d²=0 at n=6.

At n=5: d²≠0 ⟺ t3=2 (PERFECT).
At n=6: t3 alone is insufficient. Need to find additional invariant.

Key observation: all failures are at degree m=4.
Score sequence + t3 is also insufficient (6 mixed classes).

Strategy: look at the STRUCTURE of the 3-cycle graph, not just the count.
"""

import numpy as np
from itertools import combinations, permutations
from collections import Counter

def get_regular_paths(A, m):
    n = A.shape[0]
    paths = []
    def dfs(path, depth, prev):
        if depth == m:
            paths.append(tuple(path))
            return
        last = path[-1]
        for v in range(n):
            if v in path: continue
            if not A[last][v]: continue
            if depth >= 1 and not A[prev][v]: continue
            path.append(v)
            dfs(path, depth+1, last)
            path.pop()
    for start in range(n):
        dfs([start], 0, -1)
    return paths

def interior_boundary_matrix(paths_m, paths_m1):
    if not paths_m or not paths_m1:
        return np.zeros((len(paths_m1) if paths_m1 else 0,
                         len(paths_m) if paths_m else 0), dtype=int)
    idx = {p: i for i, p in enumerate(paths_m1)}
    m = len(paths_m[0]) - 1
    B = np.zeros((len(paths_m1), len(paths_m)), dtype=int)
    for j, path in enumerate(paths_m):
        for i in range(1, m):
            face = path[:i] + path[i+1:]
            if face in idx:
                B[idx[face], j] += (-1)**i
    return B

def check_d2(A):
    n = A.shape[0]
    all_paths = {m: get_regular_paths(A, m) for m in range(n)}
    for m in range(2, n):
        p_m = all_paths[m]
        p_m1 = all_paths[m-1]
        p_m2 = all_paths[m-2]
        if p_m and p_m1 and p_m2:
            B_m = interior_boundary_matrix(p_m, p_m1)
            B_m1 = interior_boundary_matrix(p_m1, p_m2)
            d2 = B_m1 @ B_m
            if np.max(np.abs(d2)) > 0:
                return False
    return True

def score_sequence(A):
    n = A.shape[0]
    return tuple(sorted([int(sum(A[v])) for v in range(n)]))

def count_3cycles(A):
    n = A.shape[0]
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (A[i][j] and A[j][k] and A[k][i]) or \
                   (A[j][i] and A[i][k] and A[k][j]):
                    count += 1
    return count

def get_3cycles(A):
    """Return list of 3-cycle vertex sets."""
    n = A.shape[0]
    cycles = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (A[i][j] and A[j][k] and A[k][i]) or \
                   (A[j][i] and A[i][k] and A[k][j]):
                    cycles.append(frozenset([i, j, k]))
    return cycles

def all_tournaments(n):
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    for bits in range(2**len(pairs)):
        A = np.zeros((n, n), dtype=int)
        for k, (i, j) in enumerate(pairs):
            if (bits >> k) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

def subtournament_iso_type(A):
    n = A.shape[0]
    min_repr = None
    for perm in permutations(range(n)):
        repr_bits = []
        for i in range(n):
            for j in range(i+1, n):
                repr_bits.append(int(A[perm[i]][perm[j]]))
        repr_tuple = tuple(repr_bits)
        if min_repr is None or repr_tuple < min_repr:
            min_repr = repr_tuple
    return min_repr

# ============================================================
print("=" * 70)
print("n=6: STRUCTURE OF 3-CYCLE GRAPH")
print("=" * 70)

# For mixed t3 values, check what STRUCTURAL property of the 3-cycle
# arrangement distinguishes pass from fail

n = 6

# 3-cycle graph properties
def cycle_graph_properties(A):
    """Compute properties of the 3-cycle intersection graph."""
    cycles = get_3cycles(A)
    t3 = len(cycles)
    if t3 == 0:
        return {'t3': 0, 'edges': 0, 'max_deg': 0, 'components': 0,
                'max_vertex_coverage': 0, 'edge_in_cycles': 0}

    # Edge count in cycle graph (cycles sharing 2 vertices)
    edges = 0
    for i in range(t3):
        for j in range(i+1, t3):
            if len(cycles[i] & cycles[j]) == 2:
                edges += 1

    # Max degree in cycle graph
    deg = [0] * t3
    for i in range(t3):
        for j in range(i+1, t3):
            if len(cycles[i] & cycles[j]) == 2:
                deg[i] += 1
                deg[j] += 1
    max_deg = max(deg) if deg else 0

    # How many vertices are in at least one 3-cycle?
    covered = set()
    for c in cycles:
        covered |= c
    coverage = len(covered)

    # How many edges are in at least one 3-cycle?
    n_vtx = A.shape[0]
    edge_in_cycles = 0
    for i in range(n_vtx):
        for j in range(i+1, n_vtx):
            for c in cycles:
                if i in c and j in c:
                    edge_in_cycles += 1
                    break

    return {'t3': t3, 'edges': edges, 'max_deg': max_deg,
            'coverage': coverage, 'edge_in_cycles': edge_in_cycles}

# Collect data
print("\nCollecting cycle graph properties for all n=6 tournaments...")
data = []
for idx, A in enumerate(all_tournaments(n)):
    d2_ok = check_d2(A)
    t3 = count_3cycles(A)
    props = cycle_graph_properties(A)
    ss = score_sequence(A)
    data.append({
        'd2_ok': d2_ok,
        't3': t3,
        'ss': ss,
        **props,
    })
    if (idx + 1) % 5000 == 0:
        print(f"  {idx+1}/32768...")

print(f"Done. Total: {len(data)}")

# For each mixed t3 value, find what distinguishes pass from fail
mixed_t3 = [2, 4, 5, 6, 7]

for t3_val in mixed_t3:
    print(f"\n{'='*50}")
    print(f"t3={t3_val}: PASS vs FAIL analysis")
    print(f"{'='*50}")

    subset = [d for d in data if d['t3'] == t3_val]
    pass_sub = [d for d in subset if d['d2_ok']]
    fail_sub = [d for d in subset if not d['d2_ok']]

    print(f"  {len(pass_sub)} pass, {len(fail_sub)} fail")

    # Compare properties
    for prop in ['edges', 'max_deg', 'coverage', 'edge_in_cycles']:
        pass_vals = Counter(d[prop] for d in pass_sub)
        fail_vals = Counter(d[prop] for d in fail_sub)
        all_vals = sorted(set(list(pass_vals.keys()) + list(fail_vals.keys())))
        print(f"\n  {prop}:")
        for v in all_vals:
            p = pass_vals.get(v, 0)
            f = fail_vals.get(v, 0)
            print(f"    {v}: {p} pass, {f} fail")

# ============================================================
print(f"\n{'=' * 70}")
print("TESTING: vertex-level cycle participation")
print("=" * 70)

# For each vertex, count how many 3-cycles it's in.
# Sort this to get a "cycle participation sequence".

def cycle_participation_seq(A):
    n_vtx = A.shape[0]
    cycles = get_3cycles(A)
    counts = [0] * n_vtx
    for c in cycles:
        for v in c:
            counts[v] += 1
    return tuple(sorted(counts))

# Check if (t3, cycle_participation_seq) characterizes d²=0
print("\nTesting: (t3, cycle_participation_seq) as characterization...")
key_d2 = {}
for d in data:
    A_idx = data.index(d)  # slow but ok for 32k
    # Reconstruct for cycle participation
    # Actually we need the tournament... let's just recompute
pass

# Redo more efficiently
print("Re-scanning with cycle participation...")
key_counts = {}  # (t3, cps) -> [pass, fail]

for A in all_tournaments(n):
    d2_ok = check_d2(A)
    t3 = count_3cycles(A)
    cps = cycle_participation_seq(A)
    key = (t3, cps)
    if key not in key_counts:
        key_counts[key] = [0, 0]
    key_counts[key][0 if d2_ok else 1] += 1

mixed_keys = [(k, v) for k, v in key_counts.items() if v[0] > 0 and v[1] > 0]
print(f"\nTotal (t3, cps) classes: {len(key_counts)}")
print(f"Mixed classes: {len(mixed_keys)}")

if mixed_keys:
    print("\n  Mixed (t3, cps) classes:")
    for (t3, cps), (p, f) in sorted(mixed_keys):
        print(f"    t3={t3}, cps={cps}: {p} pass, {f} fail")

    # Try adding more invariants
    print("\n\nTrying (t3, cps, cycle_edge_count)...")
    key_counts2 = {}
    for A in all_tournaments(n):
        d2_ok = check_d2(A)
        t3 = count_3cycles(A)
        cps = cycle_participation_seq(A)
        props = cycle_graph_properties(A)
        key = (t3, cps, props['edges'])
        if key not in key_counts2:
            key_counts2[key] = [0, 0]
        key_counts2[key][0 if d2_ok else 1] += 1

    mixed_keys2 = [(k, v) for k, v in key_counts2.items() if v[0] > 0 and v[1] > 0]
    print(f"  Mixed classes with edges added: {len(mixed_keys2)}")
    for (t3, cps, e), (p, f) in sorted(mixed_keys2):
        print(f"    t3={t3}, cps={cps}, edges={e}: {p} pass, {f} fail")
else:
    print("  *** (t3, cycle_participation_seq) PERFECTLY characterizes d²=0! ***")

# ============================================================
print(f"\n{'=' * 70}")
print("TESTING: isomorphism class as characterization")
print("=" * 70)

print("\nGrouping by tournament isomorphism class...")
iso_counts = {}
for A in all_tournaments(n):
    d2_ok = check_d2(A)
    canon = subtournament_iso_type(A)
    if canon not in iso_counts:
        iso_counts[canon] = [0, 0]
    iso_counts[canon][0 if d2_ok else 1] += 1

mixed_isos = [(k, v) for k, v in iso_counts.items() if v[0] > 0 and v[1] > 0]
print(f"Total iso classes: {len(iso_counts)}")
print(f"Mixed iso classes: {len(mixed_isos)}")

if mixed_isos:
    print("\n  *** d²=0 is NOT an isomorphism invariant at n=6! ***")
    for canon, (p, f) in mixed_isos:
        A_tmp = np.zeros((n, n), dtype=int)
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                A_tmp[i][j] = canon[idx]
                A_tmp[j][i] = 1 - canon[idx]
                idx += 1
        ss = score_sequence(A_tmp)
        t3 = count_3cycles(A_tmp)
        print(f"    ss={ss}, t3={t3}: {p} pass, {f} fail")
else:
    print("\n  d²=0 IS an isomorphism invariant at n=6 (as expected)")
    pure_fail_isos = sum(1 for _, (p, f) in iso_counts.items() if p == 0 and f > 0)
    pure_pass_isos = sum(1 for _, (p, f) in iso_counts.items() if f == 0 and p > 0)
    print(f"  Pass-only iso classes: {pure_pass_isos}")
    print(f"  Fail-only iso classes: {pure_fail_isos}")

    # For each failing iso class, show score seq and t3
    print("\n  Failing iso classes:")
    for canon, (p, f) in sorted(iso_counts.items(), key=lambda x: -x[1][1]):
        if f == 0:
            continue
        A_tmp = np.zeros((n, n), dtype=int)
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                A_tmp[i][j] = canon[idx]
                A_tmp[j][i] = 1 - canon[idx]
                idx += 1
        ss = score_sequence(A_tmp)
        t3 = count_3cycles(A_tmp)
        cps = cycle_participation_seq(A_tmp)
        print(f"    ss={ss}, t3={t3}, cps={cps}: {f} labeled")

# ============================================================
print(f"\n{'=' * 70}")
print("THE KEY QUESTION: What local pattern causes d²≠0?")
print("=" * 70)

# Since d²≠0 always at m=4, let's look at degree 4 paths specifically.
# A regular 4-path (a,b,c,d,e) requires:
#   a→b, b→c, c→d, d→e (edges)
#   a→c, b→d, c→e (skip-one)
# Interior boundary deletes positions 1,2,3:
#   Delete b: (a,c,d,e) — needs a→c (yes), c→d (yes), d→e (yes), a→d(?), c→e (yes)
#     New skip-one needed: a→d
#   Delete c: (a,b,d,e) — needs a→b (yes), b→d (yes), d→e (yes), a→d(?), b→e(?)
#     New skip-ones: a→d, b→e
#   Delete d: (a,b,c,e) — needs a→b (yes), b→c (yes), c→e (yes), a→c (yes), b→e(?)
#     New skip-one: b→e

# For d²=0, in standard simplicial:
#   ∂(delete b) produces faces deleting c,d from (a,c,d,e)
#   ∂(delete c) produces faces deleting b,d from (a,b,d,e)
#   ∂(delete d) produces faces deleting b,c from (a,b,c,e)
# The cancellation requires specific face pairings.

# The critical pattern: when does a face FAIL to be regular?
# For delete b → (a,c,d,e): fails if ¬(a→d)
# For delete c → (a,b,d,e): fails if ¬(a→d) or ¬(b→e)
# For delete d → (a,b,c,e): fails if ¬(b→e)

# So the OBSTACLE edges are a→d and b→e (the "skip-TWO" edges)
# d²=0 requires: when skip-two edges are missing, the signed
# sum of remaining face-of-face contributions still cancels.

print("""
For a regular 4-path (a,b,c,d,e), interior faces are:
  Delete pos 1 (b): (a,c,d,e) — regular iff a→d ("skip-two")
  Delete pos 2 (c): (a,b,d,e) — regular iff a→d AND b→e ("skip-two" x2)
  Delete pos 3 (d): (a,b,c,e) — regular iff b→e ("skip-two")

The NEW regularity conditions are exactly the SKIP-TWO edges: a→d and b→e.

In standard simplicial d²=0: deleting (i,j) pairs with deleting (j,i-1).
Here, if a face is MISSING, the pairing is broken.

The failure pattern: a regular 4-path where exactly one of {a→d, b→e}
holds but not the other, creating an UNPAIRED face-of-face.
""")

# Count skip-two patterns
print("Skip-two edge analysis for all n=5 regular 4-paths:")
for A in all_tournaments(5):
    d2_ok = check_d2(A)
    paths4 = get_regular_paths(A, 4)
    if not paths4:
        continue
    t3 = count_3cycles(A)
    if t3 != 2 and t3 != 3:  # compare failing (t3=2) with passing (t3=3)
        continue

    skip2_pattern = []
    for path in paths4:
        a, b, c, d, e = path
        ad = bool(A[a][d])
        be = bool(A[b][e])
        skip2_pattern.append((ad, be))

    pattern_counts = Counter(skip2_pattern)
    print(f"  t3={t3}, d²={'OK' if d2_ok else 'FAIL'}, "
          f"paths4={len(paths4)}, skip2 patterns: {dict(pattern_counts)}")
    break  # one example per class

# Do full analysis
print("\n\nFull skip-two pattern analysis:")
skip2_d2_counts = {}
for A in all_tournaments(5):
    d2_ok = check_d2(A)
    t3 = count_3cycles(A)
    paths4 = get_regular_paths(A, 4)

    skip2_summary = Counter()
    for path in paths4:
        a, b, c, d, e = path
        ad = bool(A[a][d])
        be = bool(A[b][e])
        skip2_summary[(ad, be)] += 1

    key = (t3, tuple(sorted(skip2_summary.items())))
    if key not in skip2_d2_counts:
        skip2_d2_counts[key] = [0, 0]
    skip2_d2_counts[key][0 if d2_ok else 1] += 1

print("\n  (t3, skip2_patterns) -> pass/fail:")
for (t3, pats), (p, f) in sorted(skip2_d2_counts.items()):
    status = "PASS" if f == 0 else ("FAIL" if p == 0 else "MIXED")
    print(f"    t3={t3}, {dict(pats)}: {p} pass, {f} fail [{status}]")

print("\nDONE.")
