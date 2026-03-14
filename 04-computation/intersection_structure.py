#!/usr/bin/env python3
"""
intersection_structure.py — opus-2026-03-14-S71d

KEY INSIGHT: (dc3, dc5, dc7) does NOT determine H because α₂ depends on
the intersection structure of directed cycles, not just their counts.

Two tournaments with the SAME (dc3, dc5, dc7) can have different α₂
(different numbers of vertex-disjoint pairs).

This script explores:
1. At n=5: Why is H=15 the ONLY ambiguous case?
2. At n=7: What structural feature distinguishes same-(dc3,dc5,dc7) tournaments?
3. The "intersection graph" of directed cycles and its invariants
4. Can we characterize which (dc3,dc5,dc7) tuples are α₂-ambiguous?
"""

import sys, time
import numpy as np
from itertools import combinations, permutations
from collections import defaultdict, Counter
sys.stdout.reconfigure(line_buffering=True)

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for ms in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != ms: continue
            for v in range(n):
                if not (mask & (1 << v)): continue
                pm = mask ^ (1 << v)
                t = 0
                for u in range(n):
                    if (pm & (1 << u)) and A[u][v]:
                        t += dp.get((pm, u), 0)
                if t:
                    dp[(mask, v)] = t
    return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))

def count_ham_cycles(A, n):
    if n < 3: return 0
    full_mask = (1 << n) - 1
    dp = {}
    dp[(1 << 0, 0)] = 1
    for ms in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != ms: continue
            if not (mask & 1): continue
            for v in range(n):
                if not (mask & (1 << v)): continue
                if v == 0 and ms < n: continue
                pm = mask ^ (1 << v)
                if not (pm & 1): continue
                t = 0
                for u in range(n):
                    if (pm & (1 << u)) and A[u][v]:
                        t += dp.get((pm, u), 0)
                if t:
                    dp[(mask, v)] = t
    total = 0
    for v in range(1, n):
        if A[v][0] and (full_mask, v) in dp:
            total += dp[(full_mask, v)]
    return total

def count_directed_k_cycles(A, n, k):
    if k > n: return 0
    total = 0
    for combo in combinations(range(n), k):
        verts = list(combo)
        sub = np.zeros((k, k), dtype=int)
        for i in range(k):
            for j in range(k):
                sub[i][j] = A[verts[i]][verts[j]]
        total += count_ham_cycles(sub, k)
    return total

def find_all_directed_3cycles(A, n):
    """Return list of vertex-sets that support a directed 3-cycle."""
    cycles = []
    for combo in combinations(range(n), 3):
        a, b, c = combo
        # Check both orientations
        if A[a][b] and A[b][c] and A[c][a]:
            cycles.append(frozenset(combo))
        elif A[a][c] and A[c][b] and A[b][a]:
            cycles.append(frozenset(combo))
    return cycles

def find_all_directed_5cycles(A, n):
    """Return list of (vertex-set, count) for directed 5-cycles."""
    cycles = []
    for combo in combinations(range(n), 5):
        verts = list(combo)
        sub = np.zeros((5, 5), dtype=int)
        for i in range(5):
            for j in range(5):
                sub[i][j] = A[verts[i]][verts[j]]
        hc = count_ham_cycles(sub, 5)
        if hc > 0:
            cycles.append((frozenset(combo), hc))
    return cycles

def count_disjoint_pairs(cycles_3, cycles_5=None):
    """Count vertex-disjoint pairs among all cycles.
    Each directed cycle on the same vertex-set counts separately."""
    all_cycles = []
    for vs in cycles_3:
        all_cycles.append(vs)  # 3-cycles: 1 directed per vertex-set
    if cycles_5:
        for vs, mult in cycles_5:
            for _ in range(mult):
                all_cycles.append(vs)

    pairs = 0
    for i in range(len(all_cycles)):
        for j in range(i+1, len(all_cycles)):
            if all_cycles[i].isdisjoint(all_cycles[j]):
                pairs += 1
    return pairs

# ======================================================================
# PART 1: n=5, why only H=15 is ambiguous
# ======================================================================
print("=" * 70)
print("PART 1: WHY H=15 IS THE ONLY AMBIGUOUS CASE AT n=5")
print("=" * 70)

n = 5
tb = n*(n-1)//2

# Collect all tournaments with H=15
h15_data = []
for bits in range(1 << tb):
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)
    if H != 15:
        continue
    dc3 = count_directed_k_cycles(A, n, 3)
    dc5 = count_ham_cycles(A, n)
    c3_sets = find_all_directed_3cycles(A, n)
    c5_data = find_all_directed_5cycles(A, n)

    # Count vertex-disjoint pairs (only 3-3 possible at n=5,6)
    dp33 = 0
    for i in range(len(c3_sets)):
        for j in range(i+1, len(c3_sets)):
            if c3_sets[i].isdisjoint(c3_sets[j]):
                dp33 += 1

    # Also count how many 3-cycles per vertex
    vertex_in_3cycle = Counter()
    for vs in c3_sets:
        for v in vs:
            vertex_in_3cycle[v] += 1

    a1 = dc3 + dc5
    a2 = (H - 1 - 2*a1) // 4

    h15_data.append({
        'bits': bits, 'dc3': dc3, 'dc5': dc5, 'a1': a1, 'a2': a2,
        'c3_sets': c3_sets, 'c5_data': c5_data, 'dp33': dp33,
        'vertex_coverage': dict(vertex_in_3cycle)
    })

print(f"\n  {len(h15_data)} tournaments with H=15")
groups = defaultdict(list)
for d in h15_data:
    groups[(d['dc3'], d['dc5'])].append(d)

for key in sorted(groups.keys()):
    g = groups[key]
    print(f"\n  (dc3={key[0]}, dc5={key[1]}): {len(g)} tournaments")
    for d in g[:3]:  # show first 3
        print(f"    bits={d['bits']}: α₂={d['a2']}, dp33={d['dp33']}")
        print(f"      3-cycle vertex-sets: {[set(vs) for vs in d['c3_sets']]}")
        print(f"      5-cycle data: {[(set(vs), m) for vs, m in d['c5_data']]}")

# ======================================================================
# PART 2: The intersection structure at n=7
# ======================================================================
print("\n" + "=" * 70)
print("PART 2: INTERSECTION STRUCTURE AT n=7 (ambiguous cases)")
print("=" * 70)

n = 7
tb = n*(n-1)//2
np.random.seed(42)

data7 = []
t0 = time.time()
for trial in range(500):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)
    dc3 = count_directed_k_cycles(A, n, 3)
    dc5 = count_directed_k_cycles(A, n, 5)
    dc7 = count_ham_cycles(A, n)
    a1 = dc3 + dc5 + dc7
    a2 = (H - 1 - 2*a1) // 4

    c3_sets = find_all_directed_3cycles(A, n)
    dp33 = 0
    for i in range(len(c3_sets)):
        for j in range(i+1, len(c3_sets)):
            if c3_sets[i].isdisjoint(c3_sets[j]):
                dp33 += 1

    data7.append({
        'H': H, 'dc3': dc3, 'dc5': dc5, 'dc7': dc7,
        'a1': a1, 'a2': a2, 'dp33': dp33, 'bits': bits,
        'n_3cycle_sets': len(c3_sets)
    })

    if trial % 100 == 0:
        print(f"  trial {trial}: {time.time()-t0:.1f}s")

print(f"  Done: {time.time()-t0:.1f}s")

# Ambiguous groups
groups_full = defaultdict(list)
for d in data7:
    groups_full[(d['dc3'], d['dc5'], d['dc7'])].append(d)

print(f"\n  AMBIGUOUS (dc3, dc5, dc7) → H groups:")
amb_count = 0
for key in sorted(groups_full.keys()):
    g = groups_full[key]
    h_vals = sorted(set(d['H'] for d in g))
    if len(h_vals) > 1:
        amb_count += 1
        a2_vals = sorted(set(d['a2'] for d in g))
        dp33_vals = sorted(set(d['dp33'] for d in g))
        print(f"    (dc3={key[0]}, dc5={key[1]}, dc7={key[2]}): H ∈ {h_vals}, α₂ ∈ {a2_vals}, dp33 ∈ {dp33_vals}")

total_groups = len(groups_full)
print(f"\n  {amb_count}/{total_groups} groups are ambiguous")

# ======================================================================
# PART 3: Does dp33 (disjoint 3-cycle pairs) determine α₂?
# ======================================================================
print("\n" + "=" * 70)
print("PART 3: DOES dp33 DETERMINE α₂?")
print("=" * 70)

# At n=7, α₂ = dp33 (only (3,3) disjoint pairs possible, since 3+5=8>7)
print(f"\n  At n=7: 3+5=8>7, so only (3,3) disjoint pairs contribute to α₂.")
print(f"  Therefore α₂ = dp33 (exact equality).")

mismatches = sum(1 for d in data7 if d['a2'] != d['dp33'])
print(f"  α₂ == dp33 check: {500-mismatches}/{500} match")

if mismatches > 0:
    for d in data7:
        if d['a2'] != d['dp33']:
            print(f"    MISMATCH: α₂={d['a2']}, dp33={d['dp33']}, dc3={d['dc3']}, dc5={d['dc5']}, dc7={d['dc7']}")
            break

# ======================================================================
# PART 4: What determines dp33 given dc3?
# ======================================================================
print("\n" + "=" * 70)
print("PART 4: WHAT DETERMINES dp33 GIVEN dc3?")
print("=" * 70)

groups_dc3 = defaultdict(list)
for d in data7:
    groups_dc3[d['dc3']].append(d)

print(f"\n  dc3 → possible dp33 values:")
for dc3 in sorted(groups_dc3.keys()):
    g = groups_dc3[dc3]
    dp33_vals = sorted(set(d['dp33'] for d in g))
    if len(dp33_vals) > 1:
        print(f"    dc3={dc3:2d}: dp33 ∈ {dp33_vals} ({len(dp33_vals)} values)")
    else:
        print(f"    dc3={dc3:2d}: dp33 = {dp33_vals[0]} (unique)")

# ======================================================================
# PART 5: The "spread" of 3-cycles determines dp33
# ======================================================================
print("\n" + "=" * 70)
print("PART 5: 3-CYCLE INTERSECTION GRAPH")
print("=" * 70)

# For tournaments with same dc3, what distinguishes high dp33 from low dp33?
# The key is the INTERSECTION GRAPH of 3-cycles:
# - Vertex = 3-cycle vertex-set
# - Edge = they share a vertex (non-disjoint)
# dp33 = #(non-edges) = C(dc3, 2) - #(edges in intersection graph)

print(f"\n  For each dc3 value, correlation between intersection density and dp33:")
for dc3 in sorted(groups_dc3.keys()):
    g = groups_dc3[dc3]
    if dc3 < 2 or len(g) < 5:
        continue
    # dp33 = C(dc3,2) - edges_in_int_graph
    # So dp33 range tells us about graph structure variety
    dp33_vals = [d['dp33'] for d in g]
    max_possible = dc3 * (dc3-1) // 2
    print(f"    dc3={dc3:2d}: dp33 range [{min(dp33_vals)}, {max(dp33_vals)}] out of max C({dc3},2)={max_possible}")
    print(f"             edges in int-graph range [{max_possible-max(dp33_vals)}, {max_possible-min(dp33_vals)}]")

# ======================================================================
# PART 6: Score sequence and dp33
# ======================================================================
print("\n" + "=" * 70)
print("PART 6: DOES SCORE SEQUENCE DETERMINE dp33?")
print("=" * 70)

def score_sequence(A, n):
    return tuple(sorted([sum(A[i]) for i in range(n)]))

groups_score = defaultdict(list)
for d in data7:
    A = bits_to_adj(d['bits'], n)
    ss = score_sequence(A, n)
    groups_score[ss].append(d)

print(f"\n  {len(groups_score)} distinct score sequences in 500 samples")

amb_score = 0
for ss in sorted(groups_score.keys()):
    g = groups_score[ss]
    dp33_vals = sorted(set(d['dp33'] for d in g))
    dc3_vals = sorted(set(d['dc3'] for d in g))
    if len(dp33_vals) > 1:
        amb_score += 1
        if amb_score <= 5:
            print(f"    score={ss}: dc3∈{dc3_vals}, dp33∈{dp33_vals}")

print(f"  {amb_score}/{len(groups_score)} score sequences have ambiguous dp33")

# The KEY question: does the score sequence determine dc3?
amb_dc3 = sum(1 for ss, g in groups_score.items()
              if len(set(d['dc3'] for d in g)) > 1)
print(f"  {amb_dc3}/{len(groups_score)} score sequences have ambiguous dc3")

# ======================================================================
# PART 7: The structural insight
# ======================================================================
print("\n" + "=" * 70)
print("INSIGHT: THE HIERARCHY OF DETERMINATION")
print("=" * 70)

print("""
  HIERARCHY (n=7):

  Full adjacency matrix A
       ↓ (determines everything)
  Score sequence (loses arc directions)
       ↓ (usually determines dc3 but not always)
  dc3 = #(directed 3-cycles)
       ↓ (DOES NOT determine dp33)
  dp33 = α₂ = #(vertex-disjoint 3-cycle pairs)
       ↓ (together with α₁)
  H = 1 + 2α₁ + 4α₂

  KEY INSIGHT: The OCF captures the INTERSECTION STRUCTURE
  of directed cycles, not just their counts.

  H encodes:
  - How many directed cycles exist (α₁)
  - How they are arranged in space (α₂, α₃, ...)

  This is why H is a richer invariant than any single cycle count.
  It's an independence polynomial, which captures the full
  clique/independence structure of the conflict graph.

  THE 2 CONNECTION: I(Ω, 2) = H means that evaluation at x=2
  picks up EXACTLY the right combination of intersection data
  to count Hamiltonian paths. The "2" is the unique evaluation
  point where cycle intersection structure = path structure.
""")

print("Done.")
