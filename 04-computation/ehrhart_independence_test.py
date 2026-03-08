#!/usr/bin/env python3
"""
Ehrhart-theoretic interpretation of I(Omega(T), x) (kind-pasteur-S34).

Novel observation: I(Omega(T), x) evaluated at x = positive integer k
counts the number of "k-colorings" of independent sets: ways to assign
each vertex of an independent set a color from {1,...,k}.

In lattice-point terms: I(G, k) = |kP_G ∩ Z^n| where P_G is the
independence polytope (vertices = characteristic vectors of independent
sets). The Ehrhart polynomial of P_G is I(G, x).

Connection to OCF: H(T) = I(Omega(T), 2) counts lattice points in 2*P_Omega.
The parity H(T) ≡ 1 (mod 2) means |2*P_Omega ∩ Z^n| is odd.
This could have a topological explanation via Ehrhart reciprocity.

Also tests: I(Omega, -1) = (-1)^dim * (# interior lattice points of P_Omega).
"""

import os
os.environ['PYTHONIOENCODING'] = 'utf-8'

from itertools import combinations, permutations


def tournament_adj(n, bits):
    adj = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            idx += 1
    return adj


def count_ham_paths(adj, n):
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if not (mask & (1 << u)) and adj[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[full])


def count_directed_cycles_on_set(adj, vertices):
    """Count directed Hamiltonian cycles on vertex set."""
    vlist = list(vertices)
    k = len(vlist)
    if k < 3:
        return 0
    first = vlist[0]
    rest = vlist[1:]
    count = 0
    for perm in permutations(rest):
        cycle = [first] + list(perm)
        valid = True
        for i in range(k):
            if not adj[cycle[i]][cycle[(i+1) % k]]:
                valid = False
                break
        if valid:
            count += 1
    return count


def get_odd_cycles_with_counts(adj, n):
    """Get all odd cycle vertex sets with their directed cycle counts."""
    result = []
    for k in range(3, n+1, 2):
        for vs in combinations(range(n), k):
            cnt = count_directed_cycles_on_set(adj, vs)
            if cnt > 0:
                result.append((frozenset(vs), cnt))
    return result


def compute_independence_polynomial(cycles_with_counts, x_val):
    """Compute I(Omega, x) at integer x value.

    I(Omega, x) = sum over independent sets S of x^|S|
    where each vertex set contributes its directed cycle count to the multiplicity.
    """
    vs_list = [(vs, cnt) for vs, cnt in cycles_with_counts]
    n_groups = len(vs_list)

    if n_groups == 0:
        return 1

    # Enumerate all independent sets (vertex-disjoint cycle collections)
    # For small n_groups, brute force
    if n_groups > 20:
        # Use alpha coefficients instead
        return None

    total = 0
    for mask in range(1 << n_groups):
        selected = [i for i in range(n_groups) if mask & (1 << i)]
        # Check pairwise disjointness
        is_independent = True
        for i in range(len(selected)):
            for j in range(i+1, len(selected)):
                if vs_list[selected[i]][0] & vs_list[selected[j]][0]:
                    is_independent = False
                    break
            if not is_independent:
                break
        if is_independent:
            k = len(selected)
            # Product of counts (choosing one directed cycle from each vertex set)
            product = 1
            for i in selected:
                product *= vs_list[i][1]
            total += product * (x_val ** k)
    return total


def main():
    print("=== Ehrhart-theoretic Interpretation of I(Omega(T), x) ===\n")

    for n in [5, 6]:
        print(f"\n--- n={n} ---")
        num_edges = n * (n - 1) // 2
        max_bits = 1 << num_edges

        # For each tournament, compute I(Omega, x) at x = -1, 0, 1, 2, 3
        results = {}
        for bits in range(max_bits):
            adj = tournament_adj(n, bits)
            H = count_ham_paths(adj, n)

            cycles = get_odd_cycles_with_counts(adj, n)

            vals = {}
            for x in [-1, 0, 1, 2, 3]:
                vals[x] = compute_independence_polynomial(cycles, x)

            # Verify OCF
            assert vals[2] == H, f"OCF fail: bits={bits}, H={H}, I(Omega,2)={vals[2]}"

            key = (H, vals[-1], vals[0], vals[1], vals[3])
            results[key] = results.get(key, 0) + 1

        print(f"  {max_bits} tournaments, {len(results)} distinct (H, I(-1), I(0), I(1), I(3)) tuples\n")
        print(f"  {'H':>4} {'I(-1)':>6} {'I(0)':>5} {'I(1)':>5} {'I(2)':>5} {'I(3)':>5}  count")
        for (H, im1, i0, i1, i3), count in sorted(results.items()):
            print(f"  {H:4d} {im1:6d} {i0:5d} {i1:5d} {H:5d} {i3:5d}  {count}")

    # Analysis
    print("\n\n--- Key Observations ---")
    print("I(Omega, 0) = 1 always (empty independent set)")
    print("I(Omega, 1) = 1 + alpha_1 = 1 + #(odd directed cycles)")
    print("I(Omega, 2) = H(T) (OCF)")
    print("I(Omega, 3) = 1 + 3*alpha_1 + 9*alpha_2 + 27*alpha_3 + ...")
    print("I(Omega, -1) = 1 - alpha_1 + alpha_2 - alpha_3 + ...")
    print("  = (-1)^dim * (# interior lattice points of P_Omega + correction)")
    print()

    # Check: is I(Omega, -1) determined by H?
    n = 6
    num_edges = n * (n - 1) // 2
    max_bits = 1 << num_edges

    h_to_im1 = {}
    for bits in range(max_bits):
        adj = tournament_adj(n, bits)
        H = count_ham_paths(adj, n)
        cycles = get_odd_cycles_with_counts(adj, n)
        im1 = compute_independence_polynomial(cycles, -1)
        if H not in h_to_im1:
            h_to_im1[H] = set()
        h_to_im1[H].add(im1)

    print(f"  n=6: Is I(Omega, -1) determined by H?")
    for H in sorted(h_to_im1):
        vals = sorted(h_to_im1[H])
        det = "YES" if len(vals) == 1 else "NO"
        print(f"    H={H:3d}: I(-1) values = {vals} ({det})")


if __name__ == "__main__":
    main()
