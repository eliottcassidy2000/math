"""
omega_small_ip_n8.py -- kind-pasteur-2026-03-13-S61
Find smallest directed-cycle count at n=8 for full IP decomposition.

Key discovery from omega_ip_decomposition_n8.py:
- ALL H-changing reversals at n=8 preserve vertex set counts
- But DIRECTED cycle multiplicities on vertex sets change
- The overlap weight spectrum at vertex-set level is PERFECTLY preserved
- Need a small enough example (<=25 directed cycles) for full IP

Strategy: search with different random seeds for tournaments with
fewer total directed odd cycles.
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

def count_all_directed_odd_cycles(A, n):
    """Count total directed odd cycles (fast, no enumeration)."""
    total = 0
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
                    total += 1
    return total

def enumerate_all_directed_cycles(A, n):
    """Enumerate all directed odd cycles as tuples."""
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

def build_omega(cycles):
    nc = len(cycles)
    adj = [set() for _ in range(nc)]
    vsets = [frozenset(c) for c in cycles]
    for i in range(nc):
        for j in range(i+1, nc):
            if len(vsets[i] & vsets[j]) > 0:
                adj[i].add(j)
                adj[j].add(i)
    return adj

def ip_coefficients(adj, n_v):
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

n = 8
total_bits = n * (n - 1) // 2

print("=" * 70)
print("SEARCHING FOR SMALL-OMEGA H-CHANGING EXAMPLE AT n=8")
print("=" * 70)

best_example = None
best_total = 999

for seed in range(100):
    np.random.seed(seed * 137 + 42)
    for trial in range(500):
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

            # Quick count of directed cycles
            total_dir_A = count_all_directed_odd_cycles(A, n)

            if total_dir_A < best_total:
                best_total = total_dir_A
                best_example = (bits, subset, H_A, H_B, total_dir_A)
                print(f"  NEW BEST: seed={seed}, trial={trial}, total_dir={total_dir_A}, H={H_A}->{H_B}")

                if total_dir_A <= 25:
                    print(f"  SMALL ENOUGH for full IP! Breaking.")
                    break
            break  # one per tournament

        if best_example and best_example[4] <= 25:
            break
    if best_example and best_example[4] <= 25:
        break

    if seed % 10 == 9:
        print(f"  ... seed {seed}, best so far: {best_total}")

if best_example is None:
    print("No H-changing example found!")
else:
    bits, subset, H_A, H_B, total_dir = best_example
    print(f"\nBest example: {total_dir} directed cycles")
    print(f"  bits={bits}, subset={subset}")
    print(f"  H: {H_A} -> {H_B}, delta={H_B - H_A}")

    A = bits_to_adj(bits, n)
    B = reverse_subtournament(A, n, list(subset))

    dir_A = enumerate_all_directed_cycles(A, n)
    dir_B = enumerate_all_directed_cycles(B, n)

    by_len_A = Counter(len(c) for c in dir_A)
    by_len_B = Counter(len(c) for c in dir_B)
    print(f"\n  Directed cycles A: {len(dir_A)} = {dict(sorted(by_len_A.items()))}")
    print(f"  Directed cycles B: {len(dir_B)} = {dict(sorted(by_len_B.items()))}")

    if len(dir_A) <= 25 and len(dir_B) <= 25:
        print(f"\n{'='*70}")
        print(f"FULL INDEPENDENCE POLYNOMIAL DECOMPOSITION")
        print(f"{'='*70}")

        omega_A = build_omega(dir_A)
        omega_B = build_omega(dir_B)

        edges_A = sum(len(a) for a in omega_A) // 2
        edges_B = sum(len(a) for a in omega_B) // 2
        print(f"\n  Omega edges: A={edges_A}, B={edges_B}, delta={edges_B - edges_A}")

        coeffs_A = ip_coefficients(omega_A, len(dir_A))
        coeffs_B = ip_coefficients(omega_B, len(dir_B))

        print(f"\n  IP coefficients A: {dict(sorted(coeffs_A.items()))}")
        print(f"  IP coefficients B: {dict(sorted(coeffs_B.items()))}")

        I_A = sum(coeffs_A[k] * 2**k for k in coeffs_A)
        I_B = sum(coeffs_B[k] * 2**k for k in coeffs_B)
        print(f"\n  I(Omega_A, 2) = {I_A} (should be H_A={H_A})")
        print(f"  I(Omega_B, 2) = {I_B} (should be H_B={H_B})")
        print(f"  OCF check: {'PASS' if I_A == H_A and I_B == H_B else 'FAIL'}")

        for k in sorted(set(list(coeffs_A.keys()) + list(coeffs_B.keys()))):
            delta = coeffs_B.get(k, 0) - coeffs_A.get(k, 0)
            if delta != 0:
                print(f"\n  delta_i_{k} = {delta} (contributes {delta * 2**k} to delta_H)")

        # List the actual cycles
        print(f"\n  DIRECTED CYCLES IN A:")
        for i, c in enumerate(dir_A):
            print(f"    [{i}] c{len(c)}: {c} -> verts {sorted(frozenset(c))}")

        print(f"\n  DIRECTED CYCLES IN B:")
        for i, c in enumerate(dir_B):
            print(f"    [{i}] c{len(c)}: {c} -> verts {sorted(frozenset(c))}")

        # Which cycles are shared vs different?
        set_dir_A = set(dir_A)
        set_dir_B = set(dir_B)
        common_dir = set_dir_A & set_dir_B
        lost_dir = set_dir_A - set_dir_B
        gained_dir = set_dir_B - set_dir_A

        print(f"\n  Common directed cycles: {len(common_dir)}")
        print(f"  Lost directed cycles: {len(lost_dir)}")
        for c in sorted(lost_dir, key=lambda x: (len(x), x)):
            print(f"    c{len(c)}: {c}")
        print(f"  Gained directed cycles: {len(gained_dir)}")
        for c in sorted(gained_dir, key=lambda x: (len(x), x)):
            print(f"    c{len(c)}: {c}")

        # Conflict analysis for lost vs gained
        print(f"\n  CONFLICT ANALYSIS: How lost/gained cycles interact with common ones")
        common_list = sorted(common_dir, key=lambda x: (len(x), x))
        lost_list = sorted(lost_dir, key=lambda x: (len(x), x))
        gained_list = sorted(gained_dir, key=lambda x: (len(x), x))

        # For each lost cycle, count how many common cycles it conflicts with
        for c in lost_list:
            conflicts = sum(1 for c2 in common_list if frozenset(c) & frozenset(c2))
            print(f"    LOST c{len(c)}: {c} conflicts with {conflicts}/{len(common_list)} common cycles")

        for c in gained_list:
            conflicts = sum(1 for c2 in common_list if frozenset(c) & frozenset(c2))
            print(f"    GAINED c{len(c)}: {c} conflicts with {conflicts}/{len(common_list)} common cycles")

        # Net change in conflicts
        total_lost_conflicts = sum(sum(1 for c2 in common_list if frozenset(c) & frozenset(c2)) for c in lost_list)
        total_gained_conflicts = sum(sum(1 for c2 in common_list if frozenset(c) & frozenset(c2)) for c in gained_list)
        print(f"\n  Total lost-to-common conflicts: {total_lost_conflicts}")
        print(f"  Total gained-to-common conflicts: {total_gained_conflicts}")
        print(f"  Net conflict change: {total_gained_conflicts - total_lost_conflicts}")
        print(f"  More conflicts => fewer independent sets => LOWER H")

    else:
        print(f"\n  Still too many cycles for full IP. Best found: {total_dir}")
        print(f"  The n=8 topology change needs special analysis techniques.")

        # Even without full IP, we can compute the key quantity:
        # How many vertex-DISJOINT cycle pairs exist?
        def count_disjoint_pairs(cycles):
            vsets = [frozenset(c) for c in cycles]
            count = 0
            for i in range(len(vsets)):
                for j in range(i+1, len(vsets)):
                    if len(vsets[i] & vsets[j]) == 0:
                        count += 1
            return count

        dp_A = count_disjoint_pairs(dir_A)
        dp_B = count_disjoint_pairs(dir_B)
        print(f"\n  Disjoint directed cycle pairs: A={dp_A}, B={dp_B}, delta={dp_B - dp_A}")
        print(f"  delta_H from i_1 change: {2 * (len(dir_B) - len(dir_A))}")
        print(f"  delta_H from i_2 change (if only i_1,i_2): {4 * (dp_B - dp_A)}")
        residual = (H_B - H_A) - 2 * (len(dir_B) - len(dir_A)) - 4 * (dp_B - dp_A)
        print(f"  Residual (higher order): {residual}")
        print(f"  Actual delta_H = {H_B - H_A}")
        print(f"  2*di_1 + 4*di_2 = {2 * (len(dir_B) - len(dir_A)) + 4 * (dp_B - dp_A)}")

        # Also count disjoint triples
        def count_disjoint_triples(cycles):
            vsets = [frozenset(c) for c in cycles]
            count = 0
            nc = len(vsets)
            for i in range(nc):
                for j in range(i+1, nc):
                    if vsets[i] & vsets[j]:
                        continue
                    for k in range(j+1, nc):
                        if not (vsets[i] & vsets[k]) and not (vsets[j] & vsets[k]):
                            count += 1
            return count

        dt_A = count_disjoint_triples(dir_A)
        dt_B = count_disjoint_triples(dir_B)
        print(f"\n  Disjoint directed cycle triples: A={dt_A}, B={dt_B}, delta={dt_B - dt_A}")
        full_approx = 2*(len(dir_B)-len(dir_A)) + 4*(dp_B-dp_A) + 8*(dt_B-dt_A)
        print(f"  2*di_1 + 4*di_2 + 8*di_3 = {full_approx}")
        print(f"  Remaining residual: {(H_B - H_A) - full_approx}")

print("\nDone.")
