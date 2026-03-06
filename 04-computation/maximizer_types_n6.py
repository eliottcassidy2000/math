#!/usr/bin/env python3
"""
Investigate the two types of n=6 maximizers:
- Type A: 240 with deletion spectrum (11,11,11,11,11,11)
- Type B: 240 with deletion spectrum (13,13,13,13,13,13)

What distinguishes them? Cycle counts, SC structure, sigma-orbits?

Also: test the vertex-deletion H-sum identity more carefully.

kind-pasteur-2026-03-06-S18f
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import (tournament_from_bits, hamiltonian_path_count,
                             find_odd_cycles, conflict_graph, opposite_tournament)
from itertools import permutations

def indep_poly_coeffs(adj):
    m = len(adj)
    if m == 0:
        return [1]
    nbr = [0] * m
    for i in range(m):
        for j in range(m):
            if adj[i][j]:
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

def canonical_form(T):
    n = len(T)
    best = None
    for perm in permutations(range(n)):
        form = tuple(T[perm[i]][perm[j]] for i in range(n) for j in range(n))
        if best is None or form < best:
            best = form
    return best

def is_self_converse(T):
    n = len(T)
    Top = opposite_tournament(T)
    return canonical_form(T) == canonical_form(Top)

def score_sequence(T):
    n = len(T)
    return tuple(sorted([sum(T[i]) for i in range(n)]))

# ============================================================
# Find all n=6 maximizers and classify them
# ============================================================
print("=" * 70)
print("n=6 MAXIMIZER TYPE ANALYSIS")
print("=" * 70)

n = 6
m = n * (n - 1) // 2
total = 1 << m

type_A = []  # deletion spectrum all 11
type_B = []  # deletion spectrum all 13

for bits in range(total):
    T = tournament_from_bits(n, bits)
    h = hamiltonian_path_count(T)
    if h != 45:
        continue

    # Compute deletion spectrum
    spectrum = []
    for v in range(n):
        verts = [i for i in range(n) if i != v]
        sub = [[T[verts[i]][verts[j]] for j in range(n-1)] for i in range(n-1)]
        spectrum.append(hamiltonian_path_count(sub))

    if all(s == 11 for s in spectrum):
        type_A.append(bits)
    elif all(s == 13 for s in spectrum):
        type_B.append(bits)
    else:
        print(f"  UNEXPECTED: bits={bits}, spectrum={spectrum}")

print(f"Type A (all 11): {len(type_A)} tournaments")
print(f"Type B (all 13): {len(type_B)} tournaments")

# ============================================================
# Analyze cycle structure of each type
# ============================================================
print(f"\n{'='*70}")
print("CYCLE STRUCTURE BY TYPE")
print("=" * 70)

for label, group in [("Type A (del=11)", type_A), ("Type B (del=13)", type_B)]:
    cycles_data = {}
    ip_data = {}
    sc_count = 0
    scores_seen = set()

    for bits in group[:60]:  # check a subset
        T = tournament_from_bits(n, bits)
        sc = is_self_converse(T)
        if sc:
            sc_count += 1

        cycles = find_odd_cycles(T)
        c3 = sum(1 for c in cycles if len(c) == 3)
        c5 = sum(1 for c in cycles if len(c) == 5)
        key = (c3, c5)
        cycles_data[key] = cycles_data.get(key, 0) + 1

        cg = conflict_graph(cycles)
        ip = tuple(indep_poly_coeffs(cg))
        ip_data[ip] = ip_data.get(ip, 0) + 1

        scores_seen.add(score_sequence(T))

    print(f"\n{label}:")
    print(f"  SC count (of {min(len(group), 60)} checked): {sc_count}")
    print(f"  Score sequences: {scores_seen}")
    print(f"  Cycle distributions: {cycles_data}")
    print(f"  Independence polynomials: {ip_data}")

# ============================================================
# The vertex-deletion H-sum identity
# For any tournament T: sum_v H(T-v) = ?
# Claim: sum_v H(T-v) = sum over Ham paths P of (2 + #{internal v where v_{k-1}->v_{k+1}})
# ============================================================
print(f"\n{'='*70}")
print("VERTEX-DELETION H-SUM IDENTITY")
print("=" * 70)

# Let's verify: sum_v H(T-v) by direct counting from paths
for n in [4, 5]:
    m = n * (n - 1) // 2
    print(f"\nn={n}:")

    for bits in range(min(1 << m, 20)):
        T = tournament_from_bits(n, bits)
        h = hamiltonian_path_count(T)

        # Method 1: direct sum of H(T-v)
        h_sum = 0
        for v in range(n):
            verts = [i for i in range(n) if i != v]
            sub = [[T[verts[i]][verts[j]] for j in range(n-1)] for i in range(n-1)]
            h_sum += hamiltonian_path_count(sub)

        ratio = h_sum / h if h > 0 else 0
        print(f"  bits={bits}: H={h}, sum H(T-v)={h_sum}, ratio={ratio:.4f}")

# ============================================================
# The KEY identity: sum_v H(T-v) via OCF
# H(T) = I(Omega(T), 2) = 1 + 2*alpha_1 + 4*alpha_2 + ...
# H(T-v) = I(Omega(T-v), 2) = 1 + 2*alpha_1(T-v) + ...
# sum_v H(T-v) = n + 2*sum_v alpha_1(T-v) + 4*sum_v alpha_2(T-v) + ...
#
# Now: alpha_1(T-v) = number of odd cycles in T-v (not through v)
# sum_v alpha_1(T-v) = sum_v #{odd cycles not through v}
# = sum_C (n - |C|) where sum is over odd cycles C of T
# (each cycle C misses |C| vertices, so is counted n - |C| times)
#
# Wait: sum_v alpha_1(T-v) = sum over cycles C in T of (n - |C|)
# because cycle C is NOT through v iff v not in C, so C is counted
# for each v not in C, which is n - |C| vertices.
#
# Similarly: alpha_2(T-v) = #{vertex-disjoint pairs of cycles in T-v}
# ============================================================
print(f"\n{'='*70}")
print("OCF DECOMPOSITION OF sum_v H(T-v)")
print("=" * 70)

for n in [5, 6]:
    m = n * (n - 1) // 2
    print(f"\nn={n}:")

    for bits in range(min(1 << m, 30)):
        T = tournament_from_bits(n, bits)
        h = hamiltonian_path_count(T)
        if h <= 1:
            continue

        cycles = find_odd_cycles(T)
        if not cycles:
            continue
        cg = conflict_graph(cycles)
        ip = indep_poly_coeffs(cg)

        # sum_v alpha_1(T-v) = sum_C (n - |C|)
        sum_alpha1 = sum(n - len(c) for c in cycles)

        # Verify against direct computation
        direct_sum_alpha1 = 0
        for v in range(n):
            cycles_v = [c for c in cycles if v not in c]
            direct_sum_alpha1 += len(cycles_v)

        assert sum_alpha1 == direct_sum_alpha1

        # Also compute sum_v alpha_2(T-v) directly
        direct_sum_alpha2 = 0
        for v in range(n):
            cycles_v = [c for c in cycles if v not in c]
            # Count independent pairs
            if len(cycles_v) >= 2:
                # Build conflict subgraph for cycles not through v
                cycle_idx = {id(c): i for i, c in enumerate(cycles)}
                # Actually need to reindex
                sub_cg = [[0]*len(cycles_v) for _ in range(len(cycles_v))]
                for i in range(len(cycles_v)):
                    for j in range(len(cycles_v)):
                        ci = cycles.index(cycles_v[i])
                        cj = cycles.index(cycles_v[j])
                        sub_cg[i][j] = cg[ci][cj]
                sub_ip = indep_poly_coeffs(sub_cg)
                direct_sum_alpha2 += sub_ip[2] if len(sub_ip) > 2 else 0

        # Predicted sum: n + 2*sum_alpha1 + 4*sum_alpha2 + ...
        pred_h_sum = n + 2 * sum_alpha1 + 4 * direct_sum_alpha2

        # Direct h_sum
        h_sum = 0
        for v in range(n):
            verts = [i for i in range(n) if i != v]
            sub = [[T[verts[i]][verts[j]] for j in range(n-1)] for i in range(n-1)]
            h_sum += hamiltonian_path_count(sub)

        if bits < 15:
            print(f"  bits={bits}: H={h}, sum_H(T-v)={h_sum}, "
                  f"n+2*sum_a1+4*sum_a2={pred_h_sum}, "
                  f"match={h_sum==pred_h_sum}, "
                  f"sum_a1={sum_alpha1}, sum_a2={direct_sum_alpha2}")

# ============================================================
# For the maximizers: what does the sum identity tell us?
# H(T_max) - H(T-v) = 2 * sum_{C through v} mu(C) [Claim A]
# For n=7 maximizer: H(T-v)=45 for all v, H(T)=189
# So 189 - 45 = 144 = 2 * sum mu(C through v) for each v
# => sum mu(C through v) = 72 for EVERY vertex v
#
# For n=6 type A: H(T)=45, H(T-v)=11 for all v
# 45 - 11 = 34 = 2 * sum mu(C through v)
# => sum mu(C through v) = 17 for every v
#
# For n=6 type B: H(T)=45, H(T-v)=13 for all v
# 45 - 13 = 32 = 2 * sum mu(C through v)
# => sum mu(C through v) = 16 for every v
# ============================================================
print(f"\n{'='*70}")
print("CLAIM A DECOMPOSITION FOR MAXIMIZERS")
print("=" * 70)

# For n=7, the maximizer has sum mu = 72 for every vertex
# What cycles contribute?
for bits_ex in type_A[:1]:
    T = tournament_from_bits(6, bits_ex)
    cycles = find_odd_cycles(T)
    print(f"\nType A example (bits={bits_ex}, H=45, del H=11):")
    c3 = [c for c in cycles if len(c) == 3]
    c5 = [c for c in cycles if len(c) == 5]
    print(f"  c3={len(c3)}, c5={len(c5)}")
    print(f"  Through each vertex:")
    for v in range(6):
        c3v = [c for c in c3 if v in c]
        c5v = [c for c in c5 if v in c]
        print(f"    v={v}: {len(c3v)} 3-cycles, {len(c5v)} 5-cycles through v")

for bits_ex in type_B[:1]:
    T = tournament_from_bits(6, bits_ex)
    cycles = find_odd_cycles(T)
    print(f"\nType B example (bits={bits_ex}, H=45, del H=13):")
    c3 = [c for c in cycles if len(c) == 3]
    c5 = [c for c in cycles if len(c) == 5]
    print(f"  c3={len(c3)}, c5={len(c5)}")
    print(f"  Through each vertex:")
    for v in range(6):
        c3v = [c for c in c3 if v in c]
        c5v = [c for c in c5 if v in c]
        print(f"    v={v}: {len(c3v)} 3-cycles, {len(c5v)} 5-cycles through v")

# ============================================================
# Check: is the H-sum identity exact via OCF?
# H = 1 + 2*a1 + 4*a2 + 8*a3 + ...
# sum H(T-v) = n + 2*S1 + 4*S2 + 8*S3 + ...
# where Sk = sum_v alpha_k(T-v)
#
# For maximizer: each H(T-v) is constant c
# => sum H(T-v) = n*c
# => n*c = n + 2*S1 + 4*S2 + ...
# => c = 1 + 2*S1/n + 4*S2/n + ...
# So S1/n, S2/n, etc. must be integers!
# ============================================================
print(f"\n{'='*70}")
print("INTEGRALITY CHECK FOR MAXIMIZER DELETION")
print("=" * 70)

for n_val, max_h, del_h_A, del_h_B in [(6, 45, 11, 13), (7, 189, 45, 45)]:
    print(f"\nn={n_val}: max H={max_h}")
    for del_h, label in [(del_h_A, "A"), (del_h_B, "B")]:
        diff = max_h - del_h
        print(f"  Type {label}: H(T-v)={del_h}, diff={diff} = 2*{diff//2}")
        # diff = 2 * sum mu(C through v)
        # For regular tournament: each vertex sees same number of cycles

print("\nDone.")
