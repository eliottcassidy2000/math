#!/usr/bin/env python3
"""
INV-004b: Deep analysis of cut-flip delta preservation.

Key finding from rcone_ocf_analysis.py: delta_H = delta_I for ALL cut-flips
at n=4,5. This means cut-flips preserve E(T) = H(T) - I(Omega(T), 2).

This script investigates:
1. Does delta_H = delta_I hold exhaustively at n=6?
2. Are all tournaments connected via cut-flips? (cut-flip graph connectivity)
3. What is the structure of delta_H under single-vertex cuts (S = {v})?
4. Can we express delta_H algebraically in terms of the tournament structure?

Author: opus-2026-03-06-S18
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import (tournament_from_bits, hamiltonian_path_count,
                             find_odd_cycles, conflict_graph)
from collections import defaultdict

def indep_poly_at_2(adj):
    m = len(adj)
    if m == 0:
        return 1
    nbr = [0] * m
    for i in range(m):
        for j in range(m):
            if adj[i][j]:
                nbr[i] |= 1 << j
    total = 0
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
            total += 2**bin(mask).count('1')
    return total

def compute_E(T):
    """Compute E(T) = H(T) - I(Omega(T), 2)."""
    H = hamiltonian_path_count(T)
    cycles = find_odd_cycles(T)
    if cycles:
        cg = conflict_graph(cycles)
        I2 = indep_poly_at_2(cg)
    else:
        I2 = 1
    return H, I2, H - I2

def cut_flip(T, S):
    """Apply cut-flip phi_S: reverse all arcs between S and V\\S."""
    n = len(T)
    T_new = [row[:] for row in T]
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if (i in S) != (j in S):
                T_new[i][j] = 1 - T[i][j]
    return T_new

def tournament_to_bits(T):
    """Convert tournament matrix to canonical bit representation."""
    n = len(T)
    bits = 0
    k = 0
    for i in range(n):
        for j in range(i+1, n):
            if T[i][j]:
                bits |= 1 << k
            k += 1
    return bits

print("=" * 70)
print("CUT-FLIP DELTA PRESERVATION — DEEP ANALYSIS")
print("=" * 70)

# Part 1: Exhaustive verification at n=5
print("\n--- Part 1: Exhaustive delta_H = delta_I at n=5 ---")
n = 5
m = n * (n - 1) // 2
total_flips = 0
mismatches = 0
delta_values = defaultdict(int)

for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    H, I2, E = compute_E(T)

    # Test all single-vertex cuts S = {v}
    for v in range(n):
        S = {v}
        T_flip = cut_flip(T, S)
        H_f, I2_f, E_f = compute_E(T_flip)
        dH = H_f - H
        dI = I2_f - I2
        total_flips += 1
        if dH != dI:
            mismatches += 1
        delta_values[dH] += 1

print(f"  n={n}: {total_flips} single-vertex cut-flips, mismatches: {mismatches}")
print(f"  delta_H distribution: {dict(sorted(delta_values.items()))}")

# Part 2: Single-vertex cut-flip = score change. What determines delta_H?
print("\n--- Part 2: Single-vertex cut-flip structure ---")
print("  phi_{v} reverses all arcs incident to v.")
print("  If d+(v) = out-degree of v, then after flip d+(v) -> (n-1) - d+(v).")
print("  delta_H depends on the local structure around v.")

n = 5
m = n * (n - 1) // 2
# Group delta_H by (out-degree of v, number of 3-cycles through v)
delta_by_structure = defaultdict(list)

for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    H = hamiltonian_path_count(T)

    for v in range(n):
        out_deg = sum(T[v][j] for j in range(n) if j != v)
        # Count 3-cycles through v
        c3v = 0
        for i in range(n):
            for j in range(i+1, n):
                if i == v or j == v:
                    continue
                # Check (v, i, j) forms a 3-cycle
                verts = [v, i, j]
                if T[verts[0]][verts[1]] and T[verts[1]][verts[2]] and T[verts[2]][verts[0]]:
                    c3v += 1
                elif T[verts[0]][verts[2]] and T[verts[2]][verts[1]] and T[verts[1]][verts[0]]:
                    c3v += 1

        S = {v}
        T_flip = cut_flip(T, S)
        H_f = hamiltonian_path_count(T_flip)
        dH = H_f - H
        delta_by_structure[(out_deg, c3v)].append(dH)

print(f"  (out_deg, c3_through_v) -> delta_H statistics:")
for key in sorted(delta_by_structure.keys()):
    vals = delta_by_structure[key]
    avg = sum(vals) / len(vals)
    mn, mx = min(vals), max(vals)
    unique = len(set(vals))
    print(f"    {key}: count={len(vals)}, avg={avg:.2f}, range=[{mn},{mx}], unique_values={unique}")

# Part 3: Cut-flip orbit connectivity
print("\n--- Part 3: Cut-flip orbit connectivity ---")
print("  Are all n-tournaments connected via single-vertex cut-flips?")

for n in [4, 5]:
    m = n * (n - 1) // 2
    total = 1 << m

    # BFS from tournament 0
    visited = set()
    queue = [0]
    visited.add(0)

    while queue:
        current = queue.pop(0)
        T = tournament_from_bits(n, current)

        # Try all single-vertex cuts
        for v in range(n):
            S = {v}
            T_flip = cut_flip(T, S)
            flip_bits = tournament_to_bits(T_flip)
            if flip_bits not in visited:
                visited.add(flip_bits)
                queue.append(flip_bits)

    print(f"  n={n}: reachable from T_0 via single-vertex flips: {len(visited)}/{total}")

# Part 4: Try general cut-flips for connectivity
print("\n--- Part 4: General cut-flip connectivity ---")

for n in [4, 5]:
    m = n * (n - 1) // 2
    total = 1 << m

    visited = set()
    queue = [0]
    visited.add(0)

    while queue:
        current = queue.pop(0)
        T = tournament_from_bits(n, current)

        # Try ALL non-trivial cuts
        for s_mask in range(1, (1 << n) - 1):
            S = {v for v in range(n) if (s_mask >> v) & 1}
            T_flip = cut_flip(T, S)
            flip_bits = tournament_to_bits(T_flip)
            if flip_bits not in visited:
                visited.add(flip_bits)
                queue.append(flip_bits)

    print(f"  n={n}: reachable via all cut-flips: {len(visited)}/{total}")

# Part 5: Formula for delta_H under single-vertex flip
print("\n--- Part 5: Algebraic structure of delta_H ---")
print("  phi_{v} reverses all arcs through v.")
print("  Key observation: H(T) counts Hamiltonian paths.")
print("  Flipping v changes which paths through v are valid.")

n = 4
m = n * (n - 1) // 2
print(f"\n  Detailed n={n} analysis:")
for bits in range(min(1 << m, 20)):
    T = tournament_from_bits(n, bits)
    H, I2, E = compute_E(T)
    if E != 0:
        print(f"  *** OCF FAILURE at bits={bits}!")
        continue

    deltas = []
    for v in range(n):
        S = {v}
        T_flip = cut_flip(T, S)
        H_f, I2_f, E_f = compute_E(T_flip)
        deltas.append(H_f - H)

    scores = [sum(T[i][j] for j in range(n) if j != i) for i in range(n)]
    print(f"  bits={bits:3d}: H={H:2d}, scores={scores}, delta_H={deltas}")

# Part 6: Check if delta_H for phi_{v} depends only on local info
print("\n--- Part 6: Is delta_H(phi_v) determined by the link of v? ---")
print("  The 'link' of v = sub-tournament on N+(v) union N-(v) with")
print("  the partition into out-neighbors and in-neighbors.")

n = 5
m = n * (n - 1) // 2
link_to_delta = defaultdict(set)

for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    H = hamiltonian_path_count(T)

    for v in range(n):
        # Compute link signature: encode the sub-tournament on V\{v}
        # together with the arc directions from v
        others = [u for u in range(n) if u != v]
        # Arc pattern from v
        arc_pattern = tuple(T[v][u] for u in others)
        # Sub-tournament on others
        sub_bits = 0
        k = 0
        for i in range(len(others)):
            for j in range(i+1, len(others)):
                if T[others[i]][others[j]]:
                    sub_bits |= 1 << k
                k += 1

        link_sig = (arc_pattern, sub_bits)

        S = {v}
        T_flip = cut_flip(T, S)
        H_f = hamiltonian_path_count(T_flip)
        dH = H_f - H

        link_to_delta[link_sig].add(dH)

ambiguous = sum(1 for s in link_to_delta.values() if len(s) > 1)
print(f"  n={n}: {len(link_to_delta)} distinct link signatures")
print(f"  Signatures with multiple delta_H values: {ambiguous}")
if ambiguous == 0:
    print(f"  => delta_H is FULLY DETERMINED by the link of v!")
else:
    print(f"  => delta_H is NOT determined by link alone")

print(f"\n{'='*70}")
print("ANALYSIS COMPLETE")
print("=" * 70)
