"""
ocf_direct_verify_n8.py -- kind-pasteur-2026-03-13-S61
Direct verification of OCF at n=8 and the decomposition formula.
Debug: why delta_H != 2*delta_c7?

Step 1: Verify I(Omega, 2) = H for a specific tournament
Step 2: Compute EVERY term of I(Omega, 2) before and after reversal
Step 3: Find where the discrepancy lies
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
    """Find ALL directed odd cycles as tuples (canonical form).
    Returns list of tuples, each being a directed cycle in canonical form
    (smallest vertex first, followed by next vertex in cycle direction).
    """
    cycles = []
    for k in range(3, n+1, 2):  # odd lengths only
        for combo in combinations(range(n), k):
            # Fix combo[0] as start
            for perm in permutations(combo[1:]):
                path = (combo[0],) + perm
                valid = True
                for i in range(k):
                    if A[path[i]][path[(i+1) % k]] != 1:
                        valid = False
                        break
                if valid:
                    # This is a valid directed k-cycle
                    # Canonical form: rotate so smallest vertex is first
                    # Then pick the direction where the next vertex is smaller
                    # But since we fixed combo[0] as start, and combo[0] is the smallest,
                    # this IS canonical.
                    cycles.append(path)
    return cycles

def cycle_vertex_set(cycle):
    return frozenset(cycle)

def build_omega(cycles):
    """Build conflict graph Omega. Returns adjacency as dict of sets."""
    n_cycles = len(cycles)
    adj = {i: set() for i in range(n_cycles)}
    vsets = [cycle_vertex_set(c) for c in cycles]
    for i in range(n_cycles):
        for j in range(i+1, n_cycles):
            if len(vsets[i] & vsets[j]) > 0:  # share at least one vertex
                adj[i].add(j)
                adj[j].add(i)
    return adj

def independence_polynomial_at_2(adj, n_vertices):
    """Compute I(G, 2) by enumerating all independent sets."""
    # For small graphs, brute force
    if n_vertices > 25:
        print(f"  WARNING: {n_vertices} vertices, using approximation")
        return -1

    total = 0
    for mask in range(1 << n_vertices):
        # Check if mask is an independent set
        vertices = [i for i in range(n_vertices) if mask & (1 << i)]
        is_independent = True
        for i in range(len(vertices)):
            for j in range(i+1, len(vertices)):
                if vertices[j] in adj[vertices[i]]:
                    is_independent = False
                    break
            if not is_independent:
                break
        if is_independent:
            total += 2 ** len(vertices)
    return total

print("=" * 70)
print("STEP 1: Verify I(Omega, 2) = H at n=8 for a random tournament")
print("=" * 70)

n = 8
np.random.seed(42)
bits = np.random.randint(0, 1 << (n*(n-1)//2))
A = bits_to_adj(bits, n)
H = count_ham_paths(A, n)

cycles = find_all_directed_odd_cycles(A, n)
print(f"Tournament bits={bits}")
print(f"H(T) = {H}")
print(f"Number of directed odd cycles: {len(cycles)}")
print(f"  By length: ", end="")
by_len = Counter(len(c) for c in cycles)
for k in sorted(by_len.keys()):
    print(f"c{k}={by_len[k]} ", end="")
print()

omega_adj = build_omega(cycles)
I_val = independence_polynomial_at_2(omega_adj, len(cycles))
print(f"I(Omega(T), 2) = {I_val}")
print(f"Match: {I_val == H}")

print("\n" + "=" * 70)
print("STEP 2: Find an H-changing reversal and compute I before/after")
print("=" * 70)

# Find a lambda-preserving (1,1,2,2) reversal that changes H
found = False
for trial in range(2000):
    bits = np.random.randint(0, 1 << (n*(n-1)//2))
    A = bits_to_adj(bits, n)
    lam_orig = lambda_graph(A, n)

    for subset in combinations(range(n), 4):
        ss = sub_scores(A, n, list(subset))
        if ss != (1, 1, 2, 2):
            continue
        B = reverse_subtournament(A, n, list(subset))
        lam_new = lambda_graph(B, n)
        if not np.array_equal(lam_orig, lam_new):
            continue

        H_A = count_ham_paths(A, n)
        H_B = count_ham_paths(B, n)
        if H_A != H_B:
            print(f"Found! trial={trial}, subset={subset}, H_A={H_A}, H_B={H_B}, delta={H_B-H_A}")

            # Compute all directed odd cycles before and after
            cycles_A = find_all_directed_odd_cycles(A, n)
            cycles_B = find_all_directed_odd_cycles(B, n)

            by_len_A = Counter(len(c) for c in cycles_A)
            by_len_B = Counter(len(c) for c in cycles_B)

            print(f"\nCycles BEFORE:")
            for k in sorted(by_len_A.keys()):
                print(f"  c{k} = {by_len_A[k]}")

            print(f"\nCycles AFTER:")
            for k in sorted(by_len_B.keys()):
                print(f"  c{k} = {by_len_B[k]}")

            for k in sorted(set(list(by_len_A.keys()) + list(by_len_B.keys()))):
                delta = by_len_B.get(k, 0) - by_len_A.get(k, 0)
                if delta != 0:
                    print(f"  delta_c{k} = {delta}")

            # Build Omega and compute I for both
            omega_A = build_omega(cycles_A)
            omega_B = build_omega(cycles_B)

            # Independence polynomial coefficients
            # i_k = number of independent sets of size k
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

            if len(cycles_A) <= 20 and len(cycles_B) <= 20:
                coeffs_A = ip_coefficients(omega_A, len(cycles_A))
                coeffs_B = ip_coefficients(omega_B, len(cycles_B))

                print(f"\nI.P. coefficients BEFORE: {dict(sorted(coeffs_A.items()))}")
                print(f"I.P. coefficients AFTER:  {dict(sorted(coeffs_B.items()))}")

                for k in sorted(set(list(coeffs_A.keys()) + list(coeffs_B.keys()))):
                    delta = coeffs_B.get(k, 0) - coeffs_A.get(k, 0)
                    if delta != 0:
                        print(f"  delta_i_{k} = {delta} (contributes {delta * 2**k} to delta_H)")

                I_A = sum(coeffs_A[k] * 2**k for k in coeffs_A)
                I_B = sum(coeffs_B[k] * 2**k for k in coeffs_B)
                print(f"\nI(Omega_A, 2) = {I_A} (should be {H_A})")
                print(f"I(Omega_B, 2) = {I_B} (should be {H_B})")
                print(f"delta_I = {I_B - I_A} (should be {H_B - H_A})")
            else:
                print(f"\n  Too many cycles ({len(cycles_A)}, {len(cycles_B)}) for full IP computation")

            found = True
            break
    if found:
        break

if not found:
    print("No H-changing reversal found in 2000 trials")

# Try to find a case where dH != 2*dc7
print("\n" + "=" * 70)
print("STEP 3: Find case where delta_H != 2*delta_c7")
print("=" * 70)

np.random.seed(1234)
for trial in range(3000):
    bits = np.random.randint(0, 1 << (n*(n-1)//2))
    A = bits_to_adj(bits, n)
    lam_orig = lambda_graph(A, n)

    for subset in combinations(range(n), 4):
        ss = sub_scores(A, n, list(subset))
        if ss != (1, 1, 2, 2):
            continue
        B = reverse_subtournament(A, n, list(subset))
        lam_new = lambda_graph(B, n)
        if not np.array_equal(lam_orig, lam_new):
            continue

        H_A = count_ham_paths(A, n)
        H_B = count_ham_paths(B, n)
        delta_H = H_B - H_A

        if delta_H == 0:
            continue

        # Quick c7 count (just the total)
        cycles_A = find_all_directed_odd_cycles(A, n)
        cycles_B = find_all_directed_odd_cycles(B, n)

        c7_A = sum(1 for c in cycles_A if len(c) == 7)
        c7_B = sum(1 for c in cycles_B if len(c) == 7)
        delta_c7 = c7_B - c7_A

        if delta_H != 2 * delta_c7:
            print(f"\nFOUND DISCREPANCY: trial={trial}, subset={subset}")
            print(f"  delta_H = {delta_H}, delta_c7 = {delta_c7}, 2*dc7 = {2*delta_c7}")

            by_len_A = Counter(len(c) for c in cycles_A)
            by_len_B = Counter(len(c) for c in cycles_B)

            print(f"\n  Cycles BEFORE: {dict(sorted(by_len_A.items()))}")
            print(f"  Cycles AFTER:  {dict(sorted(by_len_B.items()))}")

            for k in sorted(set(list(by_len_A.keys()) + list(by_len_B.keys()))):
                delta = by_len_B.get(k, 0) - by_len_A.get(k, 0)
                if delta != 0:
                    print(f"    delta_c{k} = {delta}")

            # Full IP decomposition
            if len(cycles_A) <= 20 and len(cycles_B) <= 20:
                omega_A = build_omega(cycles_A)
                omega_B = build_omega(cycles_B)

                coeffs_A = ip_coefficients(omega_A, len(cycles_A))
                coeffs_B = ip_coefficients(omega_B, len(cycles_B))

                print(f"\n  IP coefficients BEFORE: {dict(sorted(coeffs_A.items()))}")
                print(f"  IP coefficients AFTER:  {dict(sorted(coeffs_B.items()))}")

                for k in sorted(set(list(coeffs_A.keys()) + list(coeffs_B.keys()))):
                    delta = coeffs_B.get(k, 0) - coeffs_A.get(k, 0)
                    if delta != 0:
                        print(f"    delta_i_{k} = {delta} (contributes {delta * 2**k} to delta_H)")
            break
    else:
        continue
    break

print("\nDone.")
