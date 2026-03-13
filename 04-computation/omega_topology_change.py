"""
omega_topology_change.py -- kind-pasteur-2026-03-13-S61
The Vitali atom changes Omega's TOPOLOGY, not its vertex count.
Find examples with small Omega for full IP decomposition.

Key finding: at n=8, there exist lambda-preserving (1,1,2,2) reversals where
c3=c3', c5=c5', c7=c7' (ALL cycle counts preserved) yet H changes.
This means the CONFLICT GRAPH (which cycles share vertices) changes.
"""

import numpy as np
from itertools import combinations, permutations
from collections import Counter

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

def reverse_subtournament(A, n, subset):
    B = A.copy()
    for i in subset:
        for j in subset:
            if i != j:
                B[i][j] = A[j][i]
    return B

def lambda_graph(A, n):
    L = np.zeros((n, n), dtype=int)
    for u in range(n):
        for v in range(u+1, n):
            for w in range(n):
                if w == u or w == v:
                    continue
                if (A[u][v] and A[v][w] and A[w][u]) or (A[v][u] and A[u][w] and A[w][v]):
                    L[u][v] += 1
                    L[v][u] += 1
    return L

def sub_scores(A, n, subset):
    k = len(subset)
    return tuple(sorted([sum(A[subset[i]][subset[j]] for j in range(k) if i != j) for i in range(k)]))

def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask_size in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != mask_size:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                prev_mask = mask ^ (1 << v)
                total = 0
                for u in range(n):
                    if (prev_mask & (1 << u)) and A[u][v]:
                        total += dp.get((prev_mask, u), 0)
                if total:
                    dp[(mask, v)] = total
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def find_all_directed_odd_cycles(A, n):
    """Find ALL directed odd cycles, canonical form (smallest vertex first)."""
    cycles = []
    for k in range(3, n+1, 2):
        for combo in combinations(range(n), k):
            for perm in permutations(combo[1:]):
                path = (combo[0],) + perm
                valid = True
                for i in range(k):
                    if A[path[i]][path[(i+1) % k]] != 1:
                        valid = False
                        break
                if valid:
                    cycles.append(path)
    return cycles

def build_omega_adj(cycles):
    """Build conflict graph adjacency."""
    n_c = len(cycles)
    adj = [set() for _ in range(n_c)]
    vsets = [frozenset(c) for c in cycles]
    for i in range(n_c):
        for j in range(i+1, n_c):
            if len(vsets[i] & vsets[j]) > 0:
                adj[i].add(j)
                adj[j].add(i)
    return adj

def ip_coefficients(adj, n_v):
    """Compute independence polynomial coefficients."""
    coeffs = Counter()
    for mask in range(1 << n_v):
        vertices = [i for i in range(n_v) if mask & (1 << i)]
        is_independent = True
        for i in range(len(vertices)):
            for j in range(i+1, len(vertices)):
                if vertices[j] in adj[vertices[i]]:
                    is_independent = False
                    break
            if not is_independent:
                break
        if is_independent:
            coeffs[len(vertices)] += 1
    return coeffs

# First: find a topology-changing example at n=7 to verify the mechanism works
print("=" * 70)
print("STEP 1: Check for topology changes at n=7")
print("=" * 70)

n = 7
total_bits = n * (n - 1) // 2
np.random.seed(555)

for trial in range(2000):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    lam = lambda_graph(A, n)

    for subset in combinations(range(n), 4):
        ss = sub_scores(A, n, list(subset))
        if ss != (1, 1, 2, 2):
            continue
        B = reverse_subtournament(A, n, list(subset))
        if not np.array_equal(lam, lambda_graph(B, n)):
            continue

        H_A = count_ham_paths(A, n)
        H_B = count_ham_paths(B, n)
        if H_A == H_B:
            continue

        cycles_A = find_all_directed_odd_cycles(A, n)
        cycles_B = find_all_directed_odd_cycles(B, n)

        by_len_A = Counter(len(c) for c in cycles_A)
        by_len_B = Counter(len(c) for c in cycles_B)

        print(f"n=7 example: H={H_A}->{H_B}, cycles: {dict(by_len_A)} -> {dict(by_len_B)}")

        # Compute IP
        omega_A = build_omega_adj(cycles_A)
        omega_B = build_omega_adj(cycles_B)

        coeffs_A = ip_coefficients(omega_A, len(cycles_A))
        coeffs_B = ip_coefficients(omega_B, len(cycles_B))

        print(f"  IP before: {dict(sorted(coeffs_A.items()))}")
        print(f"  IP after:  {dict(sorted(coeffs_B.items()))}")

        for k in sorted(set(list(coeffs_A.keys()) + list(coeffs_B.keys()))):
            delta = coeffs_B.get(k, 0) - coeffs_A.get(k, 0)
            if delta != 0:
                print(f"  delta_i_{k} = {delta} (contributes {delta * 2**k})")

        I_A = sum(coeffs_A[k] * 2**k for k in coeffs_A)
        I_B = sum(coeffs_B[k] * 2**k for k in coeffs_B)
        print(f"  I_A={I_A}, I_B={I_B}, delta_I={I_B-I_A}")

        # Number of edges in Omega
        edges_A = sum(len(adj) for adj in omega_A) // 2
        edges_B = sum(len(adj) for adj in omega_B) // 2
        print(f"  Omega edges: {edges_A} -> {edges_B}, delta={edges_B-edges_A}")
        break
    else:
        continue
    break

# Now at n=8: search for topology-changing examples with small Omega
print("\n" + "=" * 70)
print("STEP 2: Find topology-changing example at n=8 with small Omega")
print("=" * 70)

n = 8
total_bits = n * (n - 1) // 2
np.random.seed(9876)

found_topo = False
for trial in range(3000):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    lam = lambda_graph(A, n)

    for subset in combinations(range(n), 4):
        ss = sub_scores(A, n, list(subset))
        if ss != (1, 1, 2, 2):
            continue
        B = reverse_subtournament(A, n, list(subset))
        if not np.array_equal(lam, lambda_graph(B, n)):
            continue

        H_A = count_ham_paths(A, n)
        H_B = count_ham_paths(B, n)
        if H_A == H_B:
            continue

        cycles_A = find_all_directed_odd_cycles(A, n)
        cycles_B = find_all_directed_odd_cycles(B, n)

        by_len_A = Counter(len(c) for c in cycles_A)
        by_len_B = Counter(len(c) for c in cycles_B)

        # Check if ALL counts preserved (topology change only)
        if by_len_A == by_len_B:
            total_cycles = len(cycles_A)
            if total_cycles <= 20:
                print(f"\nPURE TOPOLOGY CHANGE at n=8 with {total_cycles} cycles!")
                print(f"  trial={trial}, subset={subset}")
                print(f"  H={H_A}->{H_B}, delta={H_B-H_A}")
                print(f"  Cycles: {dict(sorted(by_len_A.items()))}")

                omega_A = build_omega_adj(cycles_A)
                omega_B = build_omega_adj(cycles_B)

                coeffs_A = ip_coefficients(omega_A, len(cycles_A))
                coeffs_B = ip_coefficients(omega_B, len(cycles_B))

                print(f"  IP before: {dict(sorted(coeffs_A.items()))}")
                print(f"  IP after:  {dict(sorted(coeffs_B.items()))}")

                for k in sorted(set(list(coeffs_A.keys()) + list(coeffs_B.keys()))):
                    delta = coeffs_B.get(k, 0) - coeffs_A.get(k, 0)
                    if delta != 0:
                        print(f"  delta_i_{k} = {delta} (contributes {delta * 2**k})")

                edges_A = sum(len(adj) for adj in omega_A) // 2
                edges_B = sum(len(adj) for adj in omega_B) // 2
                print(f"  Omega edges: {edges_A} -> {edges_B}, delta={edges_B-edges_A}")

                # Which specific cycles changed?
                # Compare cycle SETS
                set_A = set(cycles_A)
                set_B = set(cycles_B)
                lost = set_A - set_B
                gained = set_B - set_A
                common = set_A & set_B

                print(f"\n  Cycles: {len(common)} common, {len(lost)} lost, {len(gained)} gained")
                print(f"  Lost cycles:")
                for c in sorted(lost, key=lambda x: (len(x), x)):
                    print(f"    len={len(c)}: {c}")
                print(f"  Gained cycles:")
                for c in sorted(gained, key=lambda x: (len(x), x)):
                    print(f"    len={len(c)}: {c}")

                # Key: how do the lost/gained cycles interact with common cycles?
                # Count edges from lost cycles to common cycles
                lost_list = list(lost)
                gained_list = list(gained)
                common_list = list(common)

                def count_conflicts(cycle_list_1, cycle_list_2):
                    """Count pairs that share a vertex."""
                    count = 0
                    for c1 in cycle_list_1:
                        for c2 in cycle_list_2:
                            if len(frozenset(c1) & frozenset(c2)) > 0:
                                count += 1
                    return count

                lost_common_conflicts = count_conflicts(lost_list, common_list)
                gained_common_conflicts = count_conflicts(gained_list, common_list)
                lost_self_conflicts = count_conflicts(lost_list, lost_list) - len(lost_list)  # subtract self
                gained_self_conflicts = count_conflicts(gained_list, gained_list) - len(gained_list)

                print(f"\n  Lost-to-common conflicts: {lost_common_conflicts}")
                print(f"  Gained-to-common conflicts: {gained_common_conflicts}")
                print(f"  Delta conflicts with common: {gained_common_conflicts - lost_common_conflicts}")
                print(f"  Lost-to-lost conflicts: {lost_self_conflicts//2}")
                print(f"  Gained-to-gained conflicts: {gained_self_conflicts//2}")

                found_topo = True
                break

            elif total_cycles <= 25:
                print(f"\n  Topology change found at n=8 with {total_cycles} cycles (too large for full IP)")
                print(f"  H={H_A}->{H_B}, delta={H_B-H_A}")
                print(f"  Cycles: {dict(sorted(by_len_A.items()))}")

                # Still compute lost/gained
                set_A = set(cycles_A)
                set_B = set(cycles_B)
                lost = set_A - set_B
                gained = set_B - set_A
                common = set_A & set_B
                print(f"  {len(common)} common, {len(lost)} lost, {len(gained)} gained")

                # Show by length
                lost_by_len = Counter(len(c) for c in lost)
                gained_by_len = Counter(len(c) for c in gained)
                print(f"  Lost by length: {dict(sorted(lost_by_len.items()))}")
                print(f"  Gained by length: {dict(sorted(gained_by_len.items()))}")

    if found_topo:
        break

    if trial % 500 == 0 and trial > 0:
        print(f"  ... {trial}/3000 searched")

if not found_topo:
    print("  No small enough example found. Showing largest found.")

print("\n" + "=" * 70)
print("STEP 3: WHAT TYPE OF TOPOLOGY CHANGE OCCURS?")
print("=" * 70)
print("""
The key insight: the Vitali atom preserves cycle COUNTS but changes
cycle IDENTITIES. This changes the conflict graph's edge set:

- Lost cycles had certain conflict relationships with common cycles
- Gained cycles have DIFFERENT conflict relationships
- The net change in edges changes the independence polynomial

This is a SECOND-ORDER effect, invisible to cycle counting but visible
to the independence polynomial. It's the "hidden structure" that the
overlap weight W=0,1,2 system was designed to capture:

- W=2 (share edge): very likely to conflict
- W=1 (share vertex): conflict (by Omega definition)
- W=0 (disjoint): independent (no conflict)

The Vitali atom doesn't change the W=0 COUNTS (we proved this).
But it changes WHICH cycles are at W=0 vs W=1 with each other.
The reshuffling of conflict relationships changes i_2 (and possibly higher)
even when i_1 is preserved.

AT N=7: This effect doesn't appear because:
- 7-cycles use ALL vertices, so ALL cycles conflict with every 7-cycle
- The ONLY independent pairs are disjoint 3-cycle pairs
- The disjoint 3-cycle pair count is preserved (proved)
- So i_2 is preserved and delta_H = 2*delta_c7 exactly

AT N=8: This effect DOES appear because:
- 7-cycles use 7 of 8 vertices, leaving 1 vertex out
- A 3-cycle can be "near" that excluded vertex
- The reshuffling of 3-cycles changes which 3-cycles are
  independent of which 7-cycles (or 5-cycles)
- This changes i_2 without changing i_1
""")

print("Done.")
