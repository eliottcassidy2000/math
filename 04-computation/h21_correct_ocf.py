#!/usr/bin/env python3
"""
h21_correct_ocf.py — Correct the OCF understanding.

CRITICAL: Ω(T) has vertices = DIRECTED odd cycles, not cycle VERTEX SETS.
Multiple directed cycles on same vertex set are SEPARATE vertices.
Two directed cycles are adjacent iff they share ≥1 vertex.

So for a 5-element subtournament with 2 directed Hamiltonian cycles,
those are 2 separate vertices in Ω, and they're adjacent (same vertices!).

This changes everything about the independence polynomial.

opus-2026-03-14-S71e
"""

import sys
from itertools import combinations, permutations
from collections import defaultdict

sys.stdout.reconfigure(line_buffering=True)

def fast_hp(A, n):
    """DP Hamiltonian path count."""
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    full = (1 << n) - 1
    for mask in range(1, 1 << n):
        for v in range(n):
            c = dp.get((mask, v), 0)
            if not c or not (mask & (1 << v)):
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + c
    return sum(dp.get((full, v), 0) for v in range(n))

def get_directed_odd_cycles(A, n):
    """Get all DIRECTED odd cycles as (vertex_tuple, vertex_set) pairs."""
    cycles = []
    for length in range(3, n+1, 2):
        for verts in combinations(range(n), length):
            v0 = verts[0]
            for perm in permutations(verts[1:]):
                cycle = (v0,) + perm
                ok = True
                for i in range(length):
                    if A[cycle[i]][cycle[(i+1) % length]] != 1:
                        ok = False
                        break
                if ok:
                    # Canonical form: start from v0, so different orderings
                    # of same directed cycle are the same.
                    # But v0 is fixed, and rotation gives length equivalent cycles.
                    # For a directed cycle a→b→c→a, the rotations are:
                    # a→b→c→a, b→c→a→b, c→a→b→c — all same directed cycle.
                    # We fixed v0 = min(verts), so each directed cycle appears once.
                    cycles.append((cycle, frozenset(verts)))
    return cycles

def compute_ocf_ip(dir_cycles):
    """Compute I(Ω,2) where Ω has vertices = directed cycles,
    edges = shared vertex (adjacent in conflict graph)."""
    n_dc = len(dir_cycles)
    if n_dc == 0:
        return 1

    # Build adjacency: two directed cycles adjacent iff share ≥1 vertex
    # Two cycles on SAME vertex set are ALWAYS adjacent
    adj = [[False]*n_dc for _ in range(n_dc)]
    for i in range(n_dc):
        for j in range(i+1, n_dc):
            if dir_cycles[i][1] & dir_cycles[j][1]:  # share vertex
                adj[i][j] = adj[j][i] = True

    # I(Ω,2) = sum over independent sets S of 2^|S|
    total = 0
    for mask in range(2**n_dc):
        verts = [i for i in range(n_dc) if mask & (1<<i)]
        independent = True
        for i in range(len(verts)):
            for j in range(i+1, len(verts)):
                if adj[verts[i]][verts[j]]:
                    independent = False
                    break
            if not independent:
                break
        if independent:
            total += 2**len(verts)

    return total

print("=" * 70)
print("CORRECT OCF: DIRECTED CYCLES AS VERTICES")
print("=" * 70)

# ═══════════════════════════════════════════════════════════════════
# Part 1: Verify OCF at small n
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 1: Verify H = I(Ω,2) at n=5 ---")

n = 5
edges = [(i,j) for i in range(n) for j in range(i+1,n)]
ne = len(edges)

match = 0
mismatch = 0

for bits in range(2**ne):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if bits & (1 << idx):
            A[i][j] = 1
        else:
            A[j][i] = 1

    H = fast_hp(A, n)
    dc = get_directed_odd_cycles(A, n)
    I2 = compute_ocf_ip(dc)

    if H == I2:
        match += 1
    else:
        mismatch += 1
        if mismatch <= 3:
            print(f"  MISMATCH: bits={bits}, H={H}, I(Ω,2)={I2}")
            print(f"    Directed cycles: {len(dc)}")
            vsets = set(c[1] for c in dc)
            print(f"    Vertex sets: {len(vsets)}")

print(f"  n=5: {match} match, {mismatch} mismatch out of {2**ne}")

# ═══════════════════════════════════════════════════════════════════
# Part 2: Check where directed cycles differ from vertex sets
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 2: Directed cycles vs vertex sets ---")

n = 5
for bits in range(2**ne):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if bits & (1 << idx):
            A[i][j] = 1
        else:
            A[j][i] = 1

    dc = get_directed_odd_cycles(A, n)
    vsets = set(c[1] for c in dc)

    if len(dc) != len(vsets):
        # Multiple directed cycles on same vertex set
        from collections import Counter
        vset_count = Counter(c[1] for c in dc)
        multi = {str(sorted(k)): v for k,v in vset_count.items() if v > 1}
        print(f"  bits={bits}: {len(dc)} dir cycles, {len(vsets)} vertex sets. Multi: {multi}")
        break

print("  (Shows first tournament where dir cycle count ≠ vertex set count)")

# ═══════════════════════════════════════════════════════════════════
# Part 3: Check the (10,0) tournament at n=7
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 3: Revisiting bits=153 at n=7 ---")

n = 7
edges7 = [(i,j) for i in range(n) for j in range(i+1,n)]
ne7 = len(edges7)

bits = 153
A = [[0]*n for _ in range(n)]
for idx, (i,j) in enumerate(edges7):
    if bits & (1 << idx):
        A[i][j] = 1
    else:
        A[j][i] = 1

H = fast_hp(A, n)
dc = get_directed_odd_cycles(A, n)

print(f"  H = {H}")
print(f"  Directed odd cycles: {len(dc)}")

# Count by vertex set
from collections import Counter
vset_count = Counter(c[1] for c in dc)
vsets = set(c[1] for c in dc)
print(f"  Unique vertex sets: {len(vsets)}")

for vs, cnt in sorted(vset_count.items(), key=lambda x: (len(x[0]), x[0])):
    print(f"    {sorted(vs)} (len {len(vs)}): {cnt} directed cycle(s)")

# Now compute I(Ω,2) with correct directed-cycle definition
# This is expensive for many cycles, but let's try
if len(dc) <= 20:
    I2 = compute_ocf_ip(dc)
    print(f"\n  I(Ω,2) with directed cycles = {I2}")
    print(f"  H = {H}")
    print(f"  Match: {'YES' if H == I2 else 'NO'}")
else:
    print(f"\n  Too many cycles ({len(dc)}) for brute-force I.P. computation")
    # But we can compute since all cycles on same vertex set are adjacent
    # and cycles on intersecting vertex sets are also adjacent
    # An independent set can have at most one cycle per vertex set,
    # AND vertex sets must be pairwise disjoint.
    # So: independent sets in Ω ↔ collections of vertex-disjoint VERTEX SETS
    # with ONE directed cycle chosen per vertex set.
    # If vertex set S has d(S) directed cycles, then it contributes
    # d(S) choices (multiplied into the independent set count).

    print(f"\n  Using vertex-set decomposition:")
    print(f"  An independent set in Ω = collection of vertex-disjoint vertex sets,")
    print(f"  with one directed cycle chosen per vertex set.")
    print(f"  So i_k = Σ (over k vertex-disjoint sets) Π d(S_i)")

    # Group by vertex set
    vs_groups = defaultdict(int)
    for cycle_tuple, vs in dc:
        vs_groups[vs] += 1

    vs_list = list(vs_groups.items())
    n_vs = len(vs_list)

    # Build adjacency on vertex sets (intersecting = adjacent)
    vs_adj = [[False]*n_vs for _ in range(n_vs)]
    for i in range(n_vs):
        for j in range(i+1, n_vs):
            if vs_list[i][0] & vs_list[j][0]:
                vs_adj[i][j] = vs_adj[j][i] = True

    # I(Ω,2) = sum over independent sets of vertex sets:
    #   for each such set, multiply directed-cycle multiplicities × 2^k
    I2_correct = 0
    for mask in range(2**n_vs):
        verts_in = [i for i in range(n_vs) if mask & (1<<i)]
        independent = True
        for i in range(len(verts_in)):
            for j in range(i+1, len(verts_in)):
                if vs_adj[verts_in[i]][verts_in[j]]:
                    independent = False
                    break
            if not independent:
                break
        if independent:
            # Contribution: 2^k × product of multiplicities
            k = len(verts_in)
            mult = 1
            for idx in verts_in:
                mult *= vs_groups[vs_list[idx][0]]
            # Wait, this is wrong. Let me think...
            # An independent set of SIZE k in Ω means k directed cycles,
            # pairwise non-adjacent (pairwise vertex-disjoint).
            # For a set of k vertex-disjoint vertex sets,
            # the number of independent sets of size k in Ω
            # = product of d(S_i) over the k sets.
            # Contribution to I(Ω,2): product(d(S_i)) × 2^k
            I2_correct += mult * (2**k)

    print(f"  I(Ω,2) = {I2_correct}")
    print(f"  H = {H}")
    print(f"  Match: {'YES' if H == I2_correct else 'NO'}")
