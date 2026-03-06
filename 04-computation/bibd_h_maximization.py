#!/usr/bin/env python3
"""
BIBD STRUCTURE OF PALEY TOURNAMENTS AND H-MAXIMIZATION

Key insight: The cyclic 3-cycles of T_p form a 2-(p,3,(p+1)/4) BIBD.
This gives exact formulas for alpha_2 (disjoint 3-cycle pairs).

QUESTION: Does the BIBD structure maximize alpha_2 among all n=7 tournaments?
And does alpha_2 maximization drive H-maximization?

kind-pasteur-2026-03-06-S18h
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits, hamiltonian_path_count, find_odd_cycles, conflict_graph

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
# PHASE 1: Collect (c3, alpha_2, H) for all n=7 tournaments
# ============================================================
print("=" * 70)
print("PHASE 1: Exhaustive n=7 alpha_2 vs H analysis")
print("=" * 70)

n = 7
m = n * (n - 1) // 2

# Collect data
from collections import defaultdict
alpha2_by_c3 = defaultdict(list)  # c3 -> list of (alpha_2, H) pairs
h_maximizers = []  # (bits, H, c3, alpha_2)
max_alpha2 = 0

for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    h = hamiltonian_path_count(T)
    c3_list = find_3cycles(T)
    c3 = len(c3_list)
    a2 = count_disjoint_pairs(c3_list)

    alpha2_by_c3[c3].append((a2, h))
    if a2 > max_alpha2:
        max_alpha2 = a2

    if h == 189:
        h_maximizers.append((bits, h, c3, a2))

    if (bits + 1) % 500000 == 0:
        print(f"  {bits+1}/{1 << m} done...", flush=True)

# ============================================================
# ANALYSIS
# ============================================================
print(f"\n{'=' * 70}")
print("RESULTS")
print("=" * 70)

# Paley reference
T7 = build_paley(7)
c3_paley = len(find_3cycles(T7))
a2_paley = count_disjoint_pairs(find_3cycles(T7))
h_paley = hamiltonian_path_count(T7)
print(f"\nPaley T_7: c3={c3_paley}, alpha_2={a2_paley}, H={h_paley}")
print(f"Global max alpha_2 = {max_alpha2}")
print(f"Paley achieves global max alpha_2? {a2_paley == max_alpha2}")

# For each c3 count, show max alpha_2 and correlation with H
print(f"\nBy c3 count (max alpha_2 and correlation with H):")
for c3 in sorted(alpha2_by_c3.keys()):
    data = alpha2_by_c3[c3]
    a2_vals = [x[0] for x in data]
    h_vals = [x[1] for x in data]
    max_a2 = max(a2_vals)
    max_h = max(h_vals)
    # Find H at max alpha_2
    h_at_max_a2 = max(x[1] for x in data if x[0] == max_a2)
    marker = " <-- PALEY" if c3 == c3_paley else ""
    print(f"  c3={c3:2d}: n={len(data):6d}, max_a2={max_a2:2d}, "
          f"max_H={max_h:3d}, H@max_a2={h_at_max_a2:3d}{marker}")

# H-maximizers detail
print(f"\nH=189 maximizers ({len(h_maximizers)} total):")
c3_counts = defaultdict(int)
a2_counts = defaultdict(int)
for bits, h, c3, a2 in h_maximizers:
    c3_counts[c3] += 1
    a2_counts[a2] += 1
for c3 in sorted(c3_counts.keys()):
    print(f"  c3={c3}: {c3_counts[c3]} tournaments")
for a2 in sorted(a2_counts.keys()):
    print(f"  alpha_2={a2}: {a2_counts[a2]} tournaments")

# KEY: Among tournaments with c3=14, does max alpha_2 give max H?
print(f"\nFocused: c3=14 tournaments (same as Paley):")
c3_14 = alpha2_by_c3[14]
a2_dist = defaultdict(list)
for a2, h in c3_14:
    a2_dist[a2].append(h)
for a2 in sorted(a2_dist.keys(), reverse=True):
    h_vals = a2_dist[a2]
    marker = " <-- PALEY" if a2 == a2_paley else ""
    print(f"  alpha_2={a2}: n={len(h_vals)}, H range [{min(h_vals)},{max(h_vals)}], "
          f"avg H={sum(h_vals)/len(h_vals):.1f}{marker}")

# Correlation: alpha_2 vs H across all tournaments
all_a2 = []
all_h = []
for data in alpha2_by_c3.values():
    for a2, h in data:
        all_a2.append(a2)
        all_h.append(h)
n_data = len(all_a2)
mean_a2 = sum(all_a2) / n_data
mean_h = sum(all_h) / n_data
cov = sum((a2 - mean_a2) * (h - mean_h) for a2, h in zip(all_a2, all_h)) / n_data
std_a2 = (sum((a2 - mean_a2)**2 for a2 in all_a2) / n_data) ** 0.5
std_h = (sum((h - mean_h)**2 for h in all_h) / n_data) ** 0.5
corr = cov / (std_a2 * std_h) if std_a2 > 0 and std_h > 0 else 0
print(f"\nCorrelation(alpha_2, H) across all n=7: r = {corr:.4f}")

print("\nDone.")
