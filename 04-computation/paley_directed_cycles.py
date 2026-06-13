#!/usr/bin/env python3
"""
Fast directed cycle counting for regular n=7 tournaments using DP.
Counts directed Hamiltonian cycles on each subset via bitmask DP.

kind-pasteur-2026-03-06-S18h
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits, hamiltonian_path_count
from collections import defaultdict
from itertools import combinations

def build_paley(p):
    qr = set()
    for x in range(1, p):
        qr.add((x * x) % p)
    T = [[0]*p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in qr:
                T[i][j] = 1
    return T

def count_directed_cycles_all(T):
    """Count directed Hamiltonian cycles for ALL odd-sized subsets using DP.
    Returns dict: frozenset(verts) -> number of directed Ham cycles starting from min vertex."""
    n = len(T)
    result = {}

    for k in range(3, n+1, 2):
        for combo in combinations(range(n), k):
            verts = list(combo)
            v0 = verts[0]  # min vertex = start
            # DP: dp[mask][v] = #paths from v0 through exactly the vertices in mask, ending at v
            # mask is over LOCAL indices (0..k-1)
            dp = {}
            dp[(1, 0)] = 1  # start at local vertex 0 = v0
            for mask in range(1, 1 << k):
                if not (mask & 1):
                    continue
                for v in range(k):
                    if not (mask & (1 << v)):
                        continue
                    c = dp.get((mask, v), 0)
                    if c == 0:
                        continue
                    for u in range(k):
                        if mask & (1 << u):
                            continue
                        if T[verts[v]][verts[u]]:
                            key = (mask | (1 << u), u)
                            dp[key] = dp.get(key, 0) + c
            # Count cycles: paths visiting all vertices, ending at v with v->v0
            full = (1 << k) - 1
            num_cycles = 0
            for v in range(1, k):
                c = dp.get((full, v), 0)
                if c > 0 and T[verts[v]][verts[0]]:
                    num_cycles += c
            if num_cycles > 0:
                result[frozenset(combo)] = num_cycles
    return result

def find_3cycles(T):
    n = len(T)
    cycles = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if T[i][j] and T[j][k] and T[k][i]:
                    cycles.append((i,j,k))
                elif T[i][k] and T[k][j] and T[j][i]:
                    cycles.append((i,j,k))
    return cycles

def count_disjoint_pairs_c3(cycles):
    m = len(cycles)
    count = 0
    for i in range(m):
        si = set(cycles[i])
        for j in range(i+1, m):
            if not (si & set(cycles[j])):
                count += 1
    return count

# ============================================================
# Paley T_7 analysis
# ============================================================
print("=" * 70)
print("PALEY T_7: DIRECTED CYCLE ANALYSIS (DP)")
print("=" * 70)

T7 = build_paley(7)
h7 = hamiltonian_path_count(T7)
print(f"H(T_7) = {h7}")

dc = count_directed_cycles_all(T7)
# Group by size
by_size = defaultdict(list)
total_directed = 0
for vset, num in dc.items():
    k = len(vset)
    by_size[k].append((vset, num))
    total_directed += num

for k in sorted(by_size.keys()):
    n_sets = len(by_size[k])
    n_dir = sum(num for _, num in by_size[k])
    nums = sorted(set(num for _, num in by_size[k]))
    print(f"  {k}-cycles: {n_sets} vertex sets, {n_dir} directed cycles, "
          f"cycles/set: {nums}")

print(f"  Total vertex sets with cycles: {len(dc)}")
print(f"  Total directed cycles: {total_directed}")

# Now: in the OCF Omega(T), each DIRECTED cycle is a vertex.
# Two directed cycles are adjacent iff they share at least one vertex.
# For 3-cycles: each vertex set has 1 directed cycle, so vertex set = directed cycle.
# For 5-cycles: vertex set may have multiple directed cycles, ALL pairwise adjacent
#   (they share all 5 vertices!). So they form a clique in Omega.

# The IP of Omega counts independent sets of directed cycles.
# An independent set is a collection of pairwise vertex-disjoint directed cycles.
# At n=7: the only way to have 2 vertex-disjoint odd cycles is two 3-cycles (6 vertices).
# Each 3-cycle vertex set has exactly 1 directed cycle, so alpha_2 counts are the same.
# But alpha_1 = total directed cycles (not vertex sets)!

# So H = I(Omega, 2) = 1 + 2 * (total directed cycles) + 4 * (disjoint 3-cycle pairs)
c3 = find_3cycles(T7)
a2 = count_disjoint_pairs_c3(c3)
h_check = 1 + 2 * total_directed + 4 * a2
print(f"\nalpha_1 = {total_directed} (directed cycles)")
print(f"alpha_2 = {a2} (disjoint 3-cycle pairs)")
print(f"H check: 1 + 2*{total_directed} + 4*{a2} = {h_check}")
print(f"Match H(T_7) = {h7}? {h_check == h7}")

# ============================================================
# All regular n=7: directed cycle counts
# ============================================================
print(f"\n{'=' * 70}")
print("ALL REGULAR n=7: DIRECTED CYCLE COUNTS")
print("=" * 70)

n = 7
m = n*(n-1)//2

reg_data = []
count = 0
for bits in range(1 << m):
    T = tournament_from_bits(7, bits)
    scores = sorted(sum(T[i]) for i in range(7))
    if scores != [3]*7:
        continue
    count += 1
    h = hamiltonian_path_count(T)
    dc = count_directed_cycles_all(T)

    dc3_sets = sum(1 for v in dc if len(v) == 3)
    dc3_dir = sum(n for v, n in dc.items() if len(v) == 3)
    dc5_sets = sum(1 for v in dc if len(v) == 5)
    dc5_dir = sum(n for v, n in dc.items() if len(v) == 5)
    dc7_sets = sum(1 for v in dc if len(v) == 7)
    dc7_dir = sum(n for v, n in dc.items() if len(v) == 7)
    total_dir = dc3_dir + dc5_dir + dc7_dir

    c3 = find_3cycles(T)
    a2 = count_disjoint_pairs_c3(c3)

    reg_data.append({
        'bits': bits, 'h': h,
        'dc3_sets': dc3_sets, 'dc3_dir': dc3_dir,
        'dc5_sets': dc5_sets, 'dc5_dir': dc5_dir,
        'dc7_sets': dc7_sets, 'dc7_dir': dc7_dir,
        'total_dir': total_dir, 'a2': a2
    })

    if count % 100 == 0:
        print(f"  {count} done...", flush=True)

# Group by H
by_h = defaultdict(list)
for d in reg_data:
    by_h[d['h']].append(d)

print(f"\nTotal regular: {len(reg_data)}")
for h in sorted(by_h.keys(), reverse=True):
    entries = by_h[h]
    e = entries[0]
    # Check if all entries have same values
    dc5_dirs = set(d['dc5_dir'] for d in entries)
    dc7_dirs = set(d['dc7_dir'] for d in entries)
    total_dirs = set(d['total_dir'] for d in entries)
    a2s = set(d['a2'] for d in entries)
    print(f"  H={h}: n={len(entries)}, dc3_dir={e['dc3_dir']}, "
          f"dc5_dir={dc5_dirs}, dc7_dir={dc7_dirs}, total_dir={total_dirs}, a2={a2s}")
    # Verify H formula
    for d in entries[:1]:
        h_check = 1 + 2*d['total_dir'] + 4*d['a2']
        print(f"    check: 1 + 2*{d['total_dir']} + 4*{d['a2']} = {h_check} vs H={d['h']}")

print("\nDone.")
