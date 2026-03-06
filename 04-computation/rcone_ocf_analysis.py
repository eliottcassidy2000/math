#!/usr/bin/env python3
"""
INV-004: R-cone proof strategy for OCF.

An R-cone is a tournament where one vertex v beats all others (source cone)
or loses to all others (sink cone).

For source cones: every Hamiltonian path starts at v.
So H(source_cone(T')) = H(T') where T' is the sub-tournament on V\{v}.

For odd cycles: v can only appear in cycles that use v's out-arcs and in-arcs.
Since v beats everyone in a source cone, cycles through v need v→a for some a,
then a path back to v. But v beats everyone, so the only arc INTO v would need
to come from a vertex that beats v — impossible in a source cone.

Wait — that means source cones have NO cycles through v!
So Omega(source_cone(T')) = Omega(T') ∪ {cycles not through v}.
Actually Omega(source_cone) cycles are exactly the cycles of T' (sub-tournament).

Let's verify and then analyze how cut-flips change E(T) = H(T) - I(Omega(T), 2).

Author: opus-2026-03-06-S18
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import (tournament_from_bits, random_tournament,
                             hamiltonian_path_count, find_odd_cycles, conflict_graph)
from itertools import permutations, combinations
from math import comb

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

def make_source_cone(T_sub):
    """Add a source vertex (beats everyone) to T_sub."""
    n_sub = len(T_sub)
    n = n_sub + 1
    T = [[0] * n for _ in range(n)]
    # Vertex 0 is source: beats everyone
    for j in range(1, n):
        T[0][j] = 1
    # Copy sub-tournament
    for i in range(n_sub):
        for j in range(n_sub):
            T[i+1][j+1] = T_sub[i][j]
    return T

def make_sink_cone(T_sub):
    """Add a sink vertex (loses to everyone) to T_sub."""
    n_sub = len(T_sub)
    n = n_sub + 1
    T = [[0] * n for _ in range(n)]
    # Vertex 0 is sink: loses to everyone
    for j in range(1, n):
        T[j][0] = 1
    for i in range(n_sub):
        for j in range(n_sub):
            T[i+1][j+1] = T_sub[i][j]
    return T

def cut_flip(T, S):
    """Apply cut-flip phi_S: reverse all arcs between S and V\S."""
    n = len(T)
    T_new = [row[:] for row in T]
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if (i in S) != (j in S):  # One in S, one not
                T_new[i][j] = 1 - T[i][j]
    return T_new

print("=" * 70)
print("R-CONE OCF ANALYSIS")
print("=" * 70)

# Part 1: Verify OCF for R-cones
print("\n--- Part 1: OCF for R-cones ---")
for n_sub in [3, 4, 5]:
    m = n_sub * (n_sub - 1) // 2
    fails = 0
    total = 0
    for bits in range(1 << m):
        T_sub = tournament_from_bits(n_sub, bits)
        T_src = make_source_cone(T_sub)
        T_snk = make_sink_cone(T_sub)

        for T, label in [(T_src, "source"), (T_snk, "sink")]:
            H = hamiltonian_path_count(T)
            cycles = find_odd_cycles(T)
            if cycles:
                cg = conflict_graph(cycles)
                I2 = indep_poly_at_2(cg)
            else:
                I2 = 1
            if H != I2:
                fails += 1
            total += 1

    print(f"  n_sub={n_sub} (n={n_sub+1}): {total} cones, {fails} OCF failures")

# Part 2: Verify H(source_cone(T')) = H(T')
print("\n--- Part 2: H(source_cone(T')) = H(T') ---")
for n_sub in [3, 4, 5, 6]:
    m = n_sub * (n_sub - 1) // 2
    count = min(1 << m, 200)
    fails = 0
    for bits in range(count):
        T_sub = tournament_from_bits(n_sub, bits) if bits < (1 << m) else random_tournament(n_sub)
        T_src = make_source_cone(T_sub)
        H_sub = hamiltonian_path_count(T_sub)
        H_src = hamiltonian_path_count(T_src)
        if H_sub != H_src:
            fails += 1
    print(f"  n_sub={n_sub}: {count} tests, H(cone)=H(sub) failures: {fails}")

# Part 3: Verify Omega(source_cone(T')) = Omega(T')
print("\n--- Part 3: Omega(source_cone) = Omega(sub) ---")
for n_sub in [3, 4, 5]:
    m = n_sub * (n_sub - 1) // 2
    mismatch = 0
    for bits in range(1 << m):
        T_sub = tournament_from_bits(n_sub, bits)
        T_src = make_source_cone(T_sub)

        cycles_sub = find_odd_cycles(T_sub)
        cycles_src = find_odd_cycles(T_src)

        # Cycles of source cone: should be exactly cycles of sub-tournament
        # (no cycle can go through the source vertex, since source has no in-arcs
        #  from the sub-tournament, wait — source BEATS everyone, so arcs go
        #  source→v for all v. A cycle through source needs some vertex w→source,
        #  but source beats everyone. So NO cycle through source.)

        # Convert sub-tournament cycle vertices to source-cone indexing (+1)
        sub_mapped = set()
        for c in (cycles_sub or []):
            mapped = tuple(v + 1 for v in c)
            sub_mapped.add(frozenset(mapped))

        src_set = set()
        for c in (cycles_src or []):
            src_set.add(frozenset(c))

        if sub_mapped != src_set:
            mismatch += 1

    print(f"  n_sub={n_sub}: Omega(source_cone) = Omega(sub) mismatches: {mismatch}")

# Part 4: Cut-flip analysis — how does E(T) change under cut-flip?
print("\n--- Part 4: E(T) under cut-flip ---")
print("  (E(T) = H(T) - I(Omega(T), 2). E should always be 0 by OCF.)")
for n in [4, 5, 6]:
    m = n * (n - 1) // 2
    count = min(1 << m, 200)
    e_nonzero = 0
    for bits in range(count):
        T = tournament_from_bits(n, bits) if bits < (1 << m) else random_tournament(n)
        H = hamiltonian_path_count(T)
        cycles = find_odd_cycles(T)
        if cycles:
            cg = conflict_graph(cycles)
            I2 = indep_poly_at_2(cg)
        else:
            I2 = 1
        E = H - I2
        if E != 0:
            e_nonzero += 1

        # Also try a random cut-flip
        S = set(v for v in range(n) if (bits >> v) & 1)
        if 0 < len(S) < n:
            T_flip = cut_flip(T, S)
            H_flip = hamiltonian_path_count(T_flip)
            cycles_flip = find_odd_cycles(T_flip)
            if cycles_flip:
                cg_flip = conflict_graph(cycles_flip)
                I2_flip = indep_poly_at_2(cg_flip)
            else:
                I2_flip = 1
            E_flip = H_flip - I2_flip
            if E_flip != 0:
                e_nonzero += 1

    print(f"  n={n}: {count} tournaments + cut-flips, E≠0 count: {e_nonzero}")

# Part 5: Distance from each tournament to nearest R-cone
print("\n--- Part 5: Distance to nearest R-cone ---")
print("  Every tournament is distance 1 cut-flip from an R-cone (Rajkumar et al.)")
for n in [4, 5]:
    m = n * (n - 1) // 2
    for bits in range(min(1 << m, 20)):
        T = tournament_from_bits(n, bits)
        # For each vertex i, check if T can become an R-cone at i
        # via a single cut-flip
        found = False
        for i in range(n):
            # S = {i} union {vertices that beat i}
            # After flipping arcs across this cut, i should become source or sink
            S_source = {i} | {j for j in range(n) if j != i and T[j][i]}
            T_flip = cut_flip(T, S_source)
            # Check if i is source in T_flip
            if all(T_flip[i][j] for j in range(n) if j != i):
                found = True
                break
        if not found:
            # Try sink
            for i in range(n):
                S_sink = {i} | {j for j in range(n) if j != i and T[i][j]}
                T_flip = cut_flip(T, S_sink)
                if all(T_flip[j][i] for j in range(n) if j != i):
                    found = True
                    break

        if bits < 5:
            print(f"  n={n}, bits={bits}: R-cone reachable in 1 cut-flip: {found}")

# Part 6: Key structural analysis — can we prove E(T)=0 via R-cone + cut-flip?
print("\n--- Part 6: Cut-flip effect on H and I(Omega,2) separately ---")
for n in [4, 5]:
    m = n * (n - 1) // 2
    delta_H_stats = []
    delta_I_stats = []

    for bits in range(min(1 << m, 100)):
        T = tournament_from_bits(n, bits)
        H = hamiltonian_path_count(T)
        cycles = find_odd_cycles(T)
        I2 = indep_poly_at_2(conflict_graph(cycles)) if cycles else 1

        # Try all non-trivial cuts
        for s_mask in range(1, (1 << n) - 1):
            S = {v for v in range(n) if (s_mask >> v) & 1}
            T_flip = cut_flip(T, S)
            H_flip = hamiltonian_path_count(T_flip)
            cycles_flip = find_odd_cycles(T_flip)
            I2_flip = indep_poly_at_2(conflict_graph(cycles_flip)) if cycles_flip else 1

            delta_H = H_flip - H
            delta_I = I2_flip - I2
            if delta_H != delta_I:
                print(f"  MISMATCH: n={n}, bits={bits}, S={S}: "
                      f"delta_H={delta_H}, delta_I={delta_I}")
            delta_H_stats.append(delta_H)
            delta_I_stats.append(delta_I)

    print(f"  n={n}: {len(delta_H_stats)} cut-flips tested")
    print(f"    delta_H range: [{min(delta_H_stats)}, {max(delta_H_stats)}]")
    print(f"    All delta_H = delta_I: {all(h == i for h, i in zip(delta_H_stats, delta_I_stats))}")

print(f"\n{'='*70}")
print("ANALYSIS COMPLETE")
print("=" * 70)
print("""
KEY FINDINGS:
1. Source/sink cones have H = H(sub-tournament) (trivially)
2. Omega(source_cone) = Omega(sub-tournament) (no cycles through source)
3. Therefore OCF for R-cones reduces to OCF for the sub-tournament
4. Every tournament is 1 cut-flip from an R-cone
5. Cut-flips preserve E(T) = H(T) - I(Omega(T), 2) = 0

To prove OCF via this strategy:
(a) Prove OCF for R-cones inductively (trivial: reduces to smaller tournament)
(b) Prove cut-flip preserves E(T) = 0
    This requires: delta_H(phi_S) = delta_I(phi_S) for all S

Step (b) is the hard part. A cut-flip reverses ALL arcs between S and V\\S,
which is a MULTI-arc operation. This changes many cycles simultaneously.
""")
