#!/usr/bin/env python3
"""
Fast cycle counting for Paley T_7 vs other regular tournaments.
Uses direct Hamiltonian path/cycle counting via DP.

key question: Does Paley T_7 maximize TOTAL cycles (alpha_1)?

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

def count_cycles_by_length(T):
    """Count directed cycles of each odd length (as vertex sets, each counted once)."""
    n = len(T)
    counts = {}
    for k in range(3, n+1, 2):  # odd lengths only
        c = 0
        for combo in combinations(range(n), k):
            # Check if there's a directed Hamiltonian cycle on this subset
            verts = list(combo)
            sub = [[T[verts[i]][verts[j]] for j in range(k)] for i in range(k)]
            # DP for Ham cycles: fix vertex 0 as start, find paths visiting all, returning to 0
            full = (1 << k) - 1
            dp = {}
            dp[(1, 0)] = 1  # Start at vertex 0 (local index)
            for mask in range(1, 1 << k):
                if not (mask & 1):
                    continue  # must include vertex 0
                for v in range(k):
                    if not (mask & (1 << v)):
                        continue
                    val = dp.get((mask, v), 0)
                    if val == 0:
                        continue
                    for u in range(k):
                        if mask & (1 << u):
                            continue
                        if sub[v][u]:
                            dp[(mask | (1 << u), u)] = dp.get((mask | (1 << u), u), 0) + val
            # Count cycles: paths from 0 visiting all, ending at v with v->0
            ham_cycles = 0
            for v in range(1, k):
                if dp.get((full, v), 0) > 0 and sub[v][0]:
                    ham_cycles += dp[(full, v)]
            # Each directed cycle on k vertices is counted k times (once per starting vertex)
            # But we fixed vertex 0, so we're counting directed cycles through 0
            # Actually: we count cycles starting from vertex 0, visiting all, returning to 0
            # Each undirected cycle gives 2 directed cycles (both directions)
            # We count directed cycles starting at 0 = (total directed cycles) / k
            # Actually: each directed cycle has exactly one traversal starting at 0
            # So ham_cycles = number of directed Hamiltonian cycles on this k-subset
            # As VERTEX SETS: each vertex set supports either 0 or some number of directed cycles
            if ham_cycles > 0:
                c += 1  # Count vertex SET (not number of directed cycles)
        counts[k] = c
    return counts

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

def count_disjoint_pairs(cycles):
    m = len(cycles)
    count = 0
    for i in range(m):
        si = set(cycles[i])
        for j in range(i+1, m):
            if not (si & set(cycles[j])):
                count += 1
    return count

# ============================================================
# Paley T_7 cycle analysis
# ============================================================
print("=" * 70)
print("PALEY T_7 vs OTHER REGULAR: CYCLE COUNTS")
print("=" * 70)

T7 = build_paley(7)
h7 = hamiltonian_path_count(T7)
cycles_7 = count_cycles_by_length(T7)
c3_7_list = find_3cycles(T7)
a2_7 = count_disjoint_pairs(c3_7_list)

print(f"T_7: H={h7}")
for k, c in sorted(cycles_7.items()):
    print(f"  c_{k} = {c}")
alpha_1 = sum(cycles_7.values())
print(f"  alpha_1 = {alpha_1}")
print(f"  alpha_2 = {a2_7}")
print(f"  H check: 1 + 2*{alpha_1} + 4*{a2_7} = {1 + 2*alpha_1 + 4*a2_7}")

# ============================================================
# Sample regular tournaments by alpha_2 class
# ============================================================
print(f"\n{'=' * 70}")
print("REGULAR n=7 TOURNAMENTS: CYCLE COUNTS BY ALPHA_2 CLASS")
print("=" * 70)

n = 7
m = n*(n-1)//2

# Group by alpha_2, sample a few from each
by_a2 = defaultdict(list)
for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    scores = sorted(sum(T[i]) for i in range(n))
    if scores != [3]*7:
        continue
    c3_list = find_3cycles(T)
    a2 = count_disjoint_pairs(c3_list)
    if len(by_a2[a2]) < 10:  # Keep up to 10 per class
        by_a2[a2].append(bits)

for a2 in sorted(by_a2.keys()):
    print(f"\nalpha_2 = {a2}: {len(by_a2[a2])} samples")
    for bits in by_a2[a2][:3]:  # Analyze 3 per class
        T = tournament_from_bits(7, bits)
        h = hamiltonian_path_count(T)
        cc = count_cycles_by_length(T)
        a1 = sum(cc.values())
        h_check = 1 + 2*a1 + 4*a2
        print(f"  bits={bits}: H={h}, c3={cc.get(3,0)}, c5={cc.get(5,0)}, c7={cc.get(7,0)}, "
              f"alpha_1={a1}, check={h_check}")

# ============================================================
# Full analysis: alpha_1 by alpha_2 class
# ============================================================
print(f"\n{'=' * 70}")
print("ALPHA_1 vs ALPHA_2 FOR ALL REGULAR n=7")
print("=" * 70)

a1_by_a2 = defaultdict(list)
h_by_a2 = defaultdict(list)

for bits in range(1 << m):
    T = tournament_from_bits(7, bits)
    scores = sorted(sum(T[i]) for i in range(n))
    if scores != [3]*7:
        continue
    c3_list = find_3cycles(T)
    a2 = count_disjoint_pairs(c3_list)
    cc = count_cycles_by_length(T)
    a1 = sum(cc.values())
    h = hamiltonian_path_count(T)
    a1_by_a2[a2].append(a1)
    h_by_a2[a2].append(h)

for a2 in sorted(a1_by_a2.keys()):
    a1_vals = a1_by_a2[a2]
    h_vals = h_by_a2[a2]
    print(f"  alpha_2={a2:2d}: n={len(a1_vals):4d}, "
          f"alpha_1 range [{min(a1_vals)},{max(a1_vals)}], avg={sum(a1_vals)/len(a1_vals):.1f}, "
          f"H range [{min(h_vals)},{max(h_vals)}], avg H={sum(h_vals)/len(h_vals):.1f}")

print("\nDone.")
