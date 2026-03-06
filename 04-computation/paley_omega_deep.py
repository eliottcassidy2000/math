#!/usr/bin/env python3
"""
Deep Omega analysis for Paley T_7.
Uses tournament_lib.find_odd_cycles (fast) instead of brute-force permutation search.

kind-pasteur-2026-03-06-S18h
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits, hamiltonian_path_count, find_odd_cycles, conflict_graph
from collections import defaultdict

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

def compute_ip_from_cg(cg):
    m = len(cg)
    if m == 0:
        return [1]
    nbr = [0] * m
    for i in range(m):
        for j in range(m):
            if cg[i][j]:
                nbr[i] |= 1 << j
    coeffs = [0] * (m + 1)
    for mask in range(1 << m):
        ok = True
        seen = 0
        temp = mask
        while temp:
            v = (temp & -temp).bit_length() - 1
            if nbr[v] & seen:
                ok = False
                break
            seen |= 1 << v
            temp &= temp - 1
        if ok:
            coeffs[bin(mask).count('1')] += 1
    while len(coeffs) > 1 and coeffs[-1] == 0:
        coeffs.pop()
    return coeffs

# ============================================================
# Paley T_7 full Omega analysis
# ============================================================
print("=" * 70)
print("PALEY T_7: FULL OMEGA ANALYSIS")
print("=" * 70)

T7 = build_paley(7)
h7 = hamiltonian_path_count(T7)
print(f"H(T_7) = {h7}")

cycles = find_odd_cycles(T7)
c3 = [c for c in cycles if len(c) == 3]
c5 = [c for c in cycles if len(c) == 5]
c7 = [c for c in cycles if len(c) == 7]
print(f"Cycles: c3={len(c3)}, c5={len(c5)}, c7={len(c7)}, total={len(cycles)}")

cg = conflict_graph(cycles)
ip = compute_ip_from_cg(cg)
print(f"I(Omega(T_7), x) = {ip}")
h_check = sum(c * (2**k) for k, c in enumerate(ip))
print(f"I(Omega(T_7), 2) = {h_check}")

# ============================================================
# Disjoint pair analysis
# ============================================================
print(f"\n{'=' * 70}")
print("DISJOINT PAIRS IN OMEGA(T_7)")
print("=" * 70)

cycle_sets = [set(c) for c in cycles]
m = len(cycles)
disjoint_pairs = []
for i in range(m):
    for j in range(i+1, m):
        if not (cycle_sets[i] & cycle_sets[j]):
            disjoint_pairs.append((i, j))

print(f"Total disjoint pairs (alpha_2): {len(disjoint_pairs)}")

type_counts = defaultdict(int)
for i, j in disjoint_pairs:
    t1 = len(cycles[i])
    t2 = len(cycles[j])
    pair_type = tuple(sorted([t1, t2]))
    type_counts[pair_type] += 1

for pair_type, count in sorted(type_counts.items()):
    print(f"  ({pair_type[0]}-cycle, {pair_type[1]}-cycle): {count} disjoint pairs")

# ============================================================
# H decomposition
# ============================================================
print(f"\n{'=' * 70}")
print("H DECOMPOSITION BY CYCLE TYPE")
print("=" * 70)

alpha_1 = len(cycles)
alpha_2 = ip[2] if len(ip) > 2 else 0
H_formula = sum(c * (2**k) for k, c in enumerate(ip))

print(f"alpha_0 = 1")
print(f"alpha_1 = {alpha_1} ({len(c3)} 3-cycles + {len(c5)} 5-cycles + {len(c7)} 7-cycles)")
print(f"alpha_2 = {alpha_2}")
if len(ip) > 3:
    print(f"alpha_3 = {ip[3]}")
print(f"H = 1 + 2*{alpha_1} + 4*{alpha_2} = {H_formula}")

# ============================================================
# 5-cycle BIBD check
# ============================================================
print(f"\n{'=' * 70}")
print("5-CYCLE STRUCTURE")
print("=" * 70)

# Each 5-cycle uses 5 of 7 vertices, complement is 2 vertices
complement_count = defaultdict(int)
for cyc in c5:
    comp = tuple(sorted(set(range(7)) - set(cyc)))
    complement_count[comp] += 1

print(f"5-cycle complements (pair -> count):")
for pair in sorted(complement_count.keys()):
    print(f"  omit {pair}: {complement_count[pair]} 5-cycles")

comp_values = sorted(set(complement_count.values()))
if len(comp_values) == 1:
    lam5 = comp_values[0]
    print(f"UNIFORM! lambda_5 = {lam5}")
    # Total 5-cycles = C(7,2) * lam5 / 1 ... no
    # Each 5-cycle has C(5,2)=10 pairs of included vertices and 1 pair of excluded
    # sum of complement_count = total 5-cycles = c5
    print(f"Check: C(7,2) * {lam5} = {21 * lam5}, actual c5 = {len(c5)}")
    # The complement gives 21 (pair) slots, each with lam5 5-cycles
    # But each 5-cycle has 1 complement pair, so total = C(7,2)*lam5 is wrong
    # Actually each 5-cycle maps to exactly 1 complement pair, so
    # sum(complement_count.values()) = len(c5) which checks out
else:
    print(f"NOT uniform! Values: {comp_values}")

# Per-vertex counts
v_counts = defaultdict(lambda: defaultdict(int))
for cyc in cycles:
    k = len(cyc)
    for v in cyc:
        v_counts[v][k] += 1

print(f"\nPer-vertex cycle counts:")
for v in range(7):
    counts = v_counts[v]
    print(f"  v={v}: c3={counts.get(3,0)}, c5={counts.get(5,0)}, c7={counts.get(7,0)}")

# ============================================================
# 7-cycle (Hamiltonian cycle) analysis
# ============================================================
print(f"\n{'=' * 70}")
print("HAMILTONIAN CYCLES")
print("=" * 70)

print(f"Number of directed Hamiltonian cycles: {len(c7)}")
# Each undirected Hamiltonian cycle gives 2 directed cycles (CW and CCW)
# But find_odd_cycles might count vertex SETS (each counted once)
# Let me check
for cyc in c7[:5]:
    print(f"  {cyc}")

# ============================================================
# Compare with other n=7 maximizers
# ============================================================
print(f"\n{'=' * 70}")
print("OTHER H=189 MAXIMIZERS: CYCLE STRUCTURE")
print("=" * 70)

n = 7
m_bits = n*(n-1)//2
ip_dist = defaultdict(int)
cycle_dist = defaultdict(int)

for bits in range(1 << m_bits):
    T = tournament_from_bits(n, bits)
    h = hamiltonian_path_count(T)
    if h != 189:
        continue
    cycles_t = find_odd_cycles(T)
    ck = defaultdict(int)
    for c in cycles_t:
        ck[len(c)] += 1
    cycle_key = tuple(sorted(ck.items()))
    cycle_dist[cycle_key] += 1
    cg_t = conflict_graph(cycles_t)
    ip_t = compute_ip_from_cg(cg_t)
    ip_dist[tuple(ip_t)] += 1

print(f"IP distribution among H=189 maximizers:")
for ip_key, count in sorted(ip_dist.items()):
    h_check = sum(c * (2**k) for k, c in enumerate(ip_key))
    print(f"  IP={list(ip_key)}: {count} tournaments, I(2)={h_check}")

print(f"\nCycle structure distribution:")
for ck, count in sorted(cycle_dist.items()):
    print(f"  {dict(ck)}: {count} tournaments")

print("\nDone.")
