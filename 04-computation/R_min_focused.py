#!/usr/bin/env python3
"""
Focused R-minimization analysis: just maximizers and their U_sum breakdown.

kind-pasteur-2026-03-06-S18g
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits, hamiltonian_path_count, find_odd_cycles, conflict_graph
from collections import defaultdict

MAX_H = {1: 1, 2: 1, 3: 3, 4: 5, 5: 15, 6: 45, 7: 189}

def score_seq(T):
    return tuple(sorted(sum(T[i]) for i in range(len(T))))

def delete_vertex(T, v):
    n = len(T)
    verts = [i for i in range(n) if i != v]
    return [[T[verts[i]][verts[j]] for j in range(n-1)] for i in range(n-1)]

# ============================================================
# For each n, analyze the IP and U_sum breakdown of maximizers
# ============================================================
print("=" * 70)
print("MAXIMIZER IP AND U_sum BREAKDOWN")
print("=" * 70)

for n in [5, 6]:
    m = n * (n - 1) // 2
    max_h = MAX_H[n]

    seen_configs = set()
    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        h = hamiltonian_path_count(T)
        if h != max_h:
            continue

        s = score_seq(T)
        cycles = find_odd_cycles(T)
        cg = conflict_graph(cycles)
        m_cycles = len(cycles)

        config_key = (s, m_cycles)
        if config_key in seen_configs:
            continue
        seen_configs.add(config_key)

        # Build neighbor bitmasks
        nbr = [0] * m_cycles
        for i in range(m_cycles):
            for j in range(m_cycles):
                if cg[i][j]:
                    nbr[i] |= 1 << j

        # Enumerate independent sets
        h_by_k = defaultdict(int)
        u_by_k = defaultdict(int)
        avg_cycle_size_by_k = defaultdict(float)
        count_by_k = defaultdict(int)

        for mask in range(1 << m_cycles):
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
            if not ok:
                continue

            k = bin(mask).count('1')
            U = set()
            total_cycle_verts = 0
            temp = mask
            while temp:
                v = (temp & -temp).bit_length() - 1
                U.update(cycles[v])
                total_cycle_verts += len(cycles[v])
                temp &= temp - 1

            h_by_k[k] += 2**k
            u_by_k[k] += 2**k * len(U)
            count_by_k[k] += 1

        c3 = len([c for c in cycles if len(c) == 3])
        c5 = len([c for c in cycles if len(c) == 5])
        c7 = len([c for c in cycles if len(c) == 7])

        print(f"\nn={n}, score={s}, H={h}")
        print(f"  Cycles: c3={c3}, c5={c5}, c7={c7}, total={m_cycles}")
        print(f"  |S| | alpha_k | 2^k contr | U contr | avg |U|")
        for k in sorted(h_by_k.keys()):
            hk = h_by_k[k]
            uk = u_by_k[k]
            alpha_k = count_by_k[k]
            avg_u = uk / hk if hk > 0 else 0
            print(f"  {k:>3} | {alpha_k:>7} | {hk:>9} | {uk:>7} | {avg_u:>7.3f}")

        total_h = sum(h_by_k.values())
        total_u = sum(u_by_k.values())
        R = n - total_u / total_h
        print(f"  Total: H={total_h}, U_sum={total_u}, R={R:.6f}")

        # Compute deletion H-values
        del_hs = [hamiltonian_path_count(delete_vertex(T, v)) for v in range(n)]
        print(f"  del_hs: {del_hs}")
        print(f"  U_sum check: n*H - sum(del_hs) = {n*h - sum(del_hs)} (should be {total_u})")

# ============================================================
# Derivative identity: d/dx I(G, x)|_{x=2} relates to U_sum
# ============================================================
print(f"\n{'='*70}")
print("DERIVATIVE IDENTITY")
print("=" * 70)

# I(G, x) = sum_k alpha_k x^k
# I'(G, x) = sum_k k * alpha_k x^{k-1}
# I'(G, 2) = sum_k k * alpha_k * 2^{k-1}
# 2*I'(G,2) = sum_k k * alpha_k * 2^k
#
# But U_sum = sum_S 2^|S| * |U(S)|. This is NOT the same as
# sum_k k * alpha_k * 2^k because |U(S)| != |S|.
# |U(S)| = sum of cycle SIZES in S, not the count of cycles.
#
# For 3-cycles only: |U(S)| = 3*|S|. So U_sum = 3 * sum_k k * alpha_k * 2^k
# = 3 * 2 * I'(G, 2) = 6 * I'(G, 2).
#
# For mixed cycle sizes, we need a WEIGHTED derivative.
# Define I_w(G, x) = sum_S x^|S| * sum_{C in S} |V(C)|
# Then U_sum = I_w(G, 2).

# At n=5: only 3 and 5-cycles. But in independent sets, the cycles
# are vertex-disjoint. Two 3-cycles use 6 vertices (impossible for n=5).
# So max independent set size is:
#   - One 3-cycle (|U|=3)
#   - One 5-cycle (|U|=5)
#   - No pairs possible (n=5, any two cycles share a vertex)
# Wait, at n=5 the conflict graph is a clique? Let me check.

for n in [5]:
    m = n * (n - 1) // 2
    max_h = MAX_H[n]
    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        h = hamiltonian_path_count(T)
        if h != max_h:
            continue

        cycles = find_odd_cycles(T)
        cg = conflict_graph(cycles)
        s = score_seq(T)
        if s != (2,2,2,2,2):
            continue

        # Check if clique
        m_cycles = len(cycles)
        is_clique = all(cg[i][j] for i in range(m_cycles) for j in range(m_cycles) if i != j)
        print(f"n=5, regular: {m_cycles} cycles, clique={is_clique}")
        print(f"  IP = [1, {m_cycles}], H = 1 + 2*{m_cycles} = {1 + 2*m_cycles}")

        # For clique: only independent sets of size 0 (|U|=0) and size 1 (|U|=|V(C)|)
        # U_sum = sum_C 2 * |V(C)|
        u_sum_formula = sum(2 * len(c) for c in cycles)
        print(f"  U_sum = sum_C 2*|V(C)| = {u_sum_formula}")
        print(f"  R = {n} - {u_sum_formula}/{1+2*m_cycles} = {n - u_sum_formula/(1+2*m_cycles):.4f}")

        # For non-regular n=5 maximizer:
        break

    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        h = hamiltonian_path_count(T)
        if h != max_h:
            continue

        s = score_seq(T)
        if s == (2,2,2,2,2):
            continue

        cycles = find_odd_cycles(T)
        cg = conflict_graph(cycles)
        m_cycles = len(cycles)
        is_clique = all(cg[i][j] for i in range(m_cycles) for j in range(m_cycles) if i != j)
        print(f"\nn=5, non-regular {s}: {m_cycles} cycles, clique={is_clique}")

        if not is_clique:
            # Find independent sets of size > 1
            for i in range(m_cycles):
                for j in range(i+1, m_cycles):
                    if not cg[i][j]:
                        print(f"  Disjoint: {cycles[i]} + {cycles[j]}, |U|={len(set(cycles[i]) | set(cycles[j]))}")
        break

# ============================================================
# GENERAL FORMULA for R when Omega is a clique
# ============================================================
print(f"\n{'='*70}")
print("R FORMULA FOR CLIQUE OMEGA")
print("=" * 70)

# When Omega is a clique (all cycles pairwise share vertices):
# Only indep sets of size 0 and 1.
# H = 1 + 2*m  (m = number of cycles)
# U_sum = sum_C 2*|V(C)|
# R = n - U_sum/H = n - 2*sum|V(C)| / (1 + 2*m)
#
# For regular tournaments where every cycle is a 3-cycle:
# R = n - 2*3*m / (1 + 2*m) = n - 6m/(1+2m)
# As m increases, R -> n - 3.
#
# For non-clique Omega:
# Independent sets of size >= 2 contribute to both H and U_sum.
# The size-2 contribution: 4*alpha_2 to H, 4*sum_{disjoint pairs} |U(pair)| to U_sum.
# Since disjoint pairs have |U| = |V(C1)| + |V(C2)| >= 6 (two 3-cycles),
# their U contribution is relatively LARGE.
# The alpha_2 term boosts E[|U|] when disjoint pairs cover many vertices.

print("For clique Omega with only 3-cycles:")
for m in [5, 10, 15, 20]:
    h = 1 + 2*m
    u = 6*m
    r = 5 - u/h  # assuming n=5
    print(f"  m={m}: H={h}, U_sum={u}, R(n=5)={r:.4f}")

print("\nFor non-clique: adding alpha_2 independent pairs of 3-cycles (|U|=6 each):")
for m, a2 in [(10, 1), (10, 3), (10, 5), (14, 4), (20, 1)]:
    # H = 1 + 2m + 4*a2
    h = 1 + 2*m + 4*a2
    # U_sum: size-0 contributes 0, size-1 contributes 2*3m=6m,
    # size-2 contributes 4*a2*6 = 24*a2
    u = 6*m + 24*a2
    r = 6 - u/h  # n=6
    print(f"  m={m}, a2={a2}: H={h}, U_sum~={u}, R(n=6)={r:.4f}")

print("\nDone.")
