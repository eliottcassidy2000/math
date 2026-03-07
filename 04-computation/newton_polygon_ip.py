#!/usr/bin/env python3
"""
Newton Polygon of I(Omega(T), x) and H-gap Analysis (kind-pasteur-S34).

FIXED: Count individual directed cycles, not just vertex sets.
A 5-vertex set can have 0, 1, or 2+ directed Hamiltonian 5-cycles.

Novel hypothesis: The Newton polygon of I(Omega(T), x) constrains
which H values are achievable.
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


def count_directed_cycles(adj, vertices):
    """Count the number of directed Hamiltonian cycles on a given vertex set.
    A directed cycle a->b->c->...->a is counted once per cyclic ordering.
    For k vertices, there are (k-1)! directed cycles total (fixing one vertex).
    We check each and count how many are valid in the tournament."""
    vlist = list(vertices)
    k = len(vlist)
    if k < 3:
        return 0
    # Fix first vertex, permute rest
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


def get_all_odd_cycles(adj, n):
    """Get all directed odd cycles as (vertex_set, count) pairs.
    Returns list of (frozenset, num_directed_cycles) for each vertex set."""
    cycles = []
    # 3-cycles: each cyclic 3-set has exactly 1 directed 3-cycle
    for vs in combinations(range(n), 3):
        a, b, c = vs
        if (adj[a][b] and adj[b][c] and adj[c][a]) or \
           (adj[a][c] and adj[c][b] and adj[b][a]):
            cycles.append((frozenset(vs), 1))
    # 5-cycles
    if n >= 5:
        for vs in combinations(range(n), 5):
            cnt = count_directed_cycles(adj, vs)
            if cnt > 0:
                cycles.append((frozenset(vs), cnt))
    # 7-cycles (expensive, skip for n > 7)
    if n == 7:
        for vs in combinations(range(n), 7):
            cnt = count_directed_cycles(adj, vs)
            if cnt > 0:
                cycles.append((frozenset(vs), cnt))
    return cycles


def compute_alpha_coefficients(cycles_with_counts):
    """Compute alpha_k coefficients for I(Omega, x).

    alpha_k = number of collections of k pairwise vertex-disjoint directed cycles.

    For the purpose of alpha_k, each vertex set with c directed cycles
    contributes c choices to alpha_1 (as individual cycles),
    but for alpha_2+, pairs must be on disjoint vertex sets, and we get
    c1*c2 pairs from two sets with c1 and c2 cycles respectively.
    """
    # Group by vertex set
    vs_list = [(vs, cnt) for vs, cnt in cycles_with_counts]
    n_groups = len(vs_list)

    # alpha_0 = 1 always
    # alpha_1 = sum of counts
    alpha_1 = sum(cnt for _, cnt in vs_list)

    # alpha_2 = sum over disjoint pairs of cnt_i * cnt_j
    alpha_2 = 0
    for i in range(n_groups):
        for j in range(i+1, n_groups):
            if not (vs_list[i][0] & vs_list[j][0]):  # disjoint
                alpha_2 += vs_list[i][1] * vs_list[j][1]

    # alpha_3 = sum over disjoint triples
    alpha_3 = 0
    for i in range(n_groups):
        for j in range(i+1, n_groups):
            if vs_list[i][0] & vs_list[j][0]:
                continue
            for k in range(j+1, n_groups):
                if not (vs_list[k][0] & vs_list[i][0]) and \
                   not (vs_list[k][0] & vs_list[j][0]):
                    alpha_3 += vs_list[i][1] * vs_list[j][1] * vs_list[k][1]

    return [1, alpha_1, alpha_2, alpha_3]


def v2(n):
    if n == 0:
        return float('inf')
    v = 0
    n = abs(n)
    while n % 2 == 0:
        n //= 2
        v += 1
    return v


def main():
    print("=== Newton Polygon of I(Omega(T), x) — FIXED ===\n")

    for n in [5, 6]:
        print(f"\n--- n={n} ---")
        num_edges = n * (n - 1) // 2
        max_bits = 1 << num_edges

        h_to_alphas = {}
        violations = 0

        for bits in range(max_bits):
            adj = tournament_adj(n, bits)
            H = count_ham_paths(adj, n)

            cycles = get_all_odd_cycles(adj, n)
            alpha = compute_alpha_coefficients(cycles)

            # Compute I(Omega, 2)
            ip_at_2 = sum(a * (2**k) for k, a in enumerate(alpha))

            if ip_at_2 != H:
                violations += 1
                if violations <= 5:
                    print(f"  OCF check: bits={bits}, H={H}, I(Omega,2)={ip_at_2}, "
                          f"alpha={alpha}")

            alpha_tuple = tuple(alpha[:4])
            if H not in h_to_alphas:
                h_to_alphas[H] = set()
            h_to_alphas[H].add(alpha_tuple)

        print(f"  {max_bits} tournaments, {violations} OCF mismatches")
        print(f"  {len(h_to_alphas)} distinct H values\n")

        for H in sorted(h_to_alphas):
            alphas_set = h_to_alphas[H]
            if len(alphas_set) <= 5:
                for alpha in sorted(alphas_set):
                    terms = [f"{a}*2^{k}" for k, a in enumerate(alpha) if a > 0]
                    print(f"    H={H:4d} (v2={v2(H)}): alpha={list(alpha)}, "
                          f"I = {' + '.join(terms)}")
            else:
                a1_vals = sorted(set(a[1] for a in alphas_set))
                print(f"    H={H:4d}: {len(alphas_set)} alpha tuples, "
                      f"alpha1 in {a1_vals}")

    # Part 2: What alpha tuples would give H=7 or H=21?
    print("\n\n--- What would H=7 or H=21 require? ---")
    for target_H in [7, 21, 23]:
        print(f"\n  H={target_H}: need 1 + 2*a1 + 4*a2 + 8*a3 = {target_H}")
        remaining = target_H - 1
        solutions = []
        for a1 in range(remaining // 2 + 1):
            r1 = remaining - 2 * a1
            for a2 in range(r1 // 4 + 1):
                r2 = r1 - 4 * a2
                for a3 in range(r2 // 8 + 1):
                    if r2 - 8 * a3 == 0:
                        solutions.append((1, a1, a2, a3))

        print(f"  {len(solutions)} alpha decompositions:")
        for sol in solutions:
            print(f"    alpha = {list(sol)}: "
                  f"a1={sol[1]} cycles, a2={sol[2]} disjoint pairs, "
                  f"a3={sol[3]} disjoint triples")

    # Part 3: n=7 sampling
    print("\n\n--- n=7 (sampling) ---")
    import random
    random.seed(42)
    n = 7
    num_edges = n * (n - 1) // 2
    max_bits = 1 << num_edges

    h_to_alphas = {}
    for trial in range(2000):
        bits = random.randint(0, max_bits - 1)
        adj = tournament_adj(n, bits)
        H = count_ham_paths(adj, n)

        cycles = get_all_odd_cycles(adj, n)
        alpha = compute_alpha_coefficients(cycles)

        ip_at_2 = sum(a * (2**k) for k, a in enumerate(alpha))
        if ip_at_2 != H and trial < 5:
            print(f"  MISMATCH: trial={trial}, H={H}, I={ip_at_2}, alpha={alpha}")

        if H not in h_to_alphas:
            h_to_alphas[H] = set()
        h_to_alphas[H].add(tuple(alpha[:4]))

    print(f"  2000 samples, {len(h_to_alphas)} distinct H values")

    # Show alpha structure for small H values
    for H in sorted(h_to_alphas):
        if H <= 25 or H in [7, 21, 23]:
            alphas = h_to_alphas[H]
            print(f"    H={H:4d}: {len(alphas)} alpha tuples: {sorted(alphas)[:5]}")


if __name__ == "__main__":
    main()
