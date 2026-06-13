#!/usr/bin/env python3
"""
Fast Ehrhart/IP analysis of I(Omega(T), x) (kind-pasteur-S34).

Uses optimized cycle counting: for k-cycles, fix first vertex and
check all (k-1)! orderings of the rest. For k=3, only 2 orderings.
For k=5, 24 orderings. For k=7, 720 orderings.

Focus on n=5 (exact) and n=6 (fast enough with 3-cycles only for alpha).
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


def get_3cycle_sets(adj, n):
    """Fast: get all 3-vertex sets that form a directed 3-cycle."""
    result = []
    for a, b, c in combinations(range(n), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or \
           (adj[a][c] and adj[c][b] and adj[b][a]):
            result.append(frozenset([a, b, c]))
    return result


def count_5cycles(adj, vs):
    """Count directed 5-cycles on a 5-vertex set. Fast check."""
    vl = list(vs)
    first = vl[0]
    count = 0
    for p in permutations(vl[1:]):
        cycle = [first] + list(p)
        if adj[cycle[0]][cycle[1]] and adj[cycle[1]][cycle[2]] and \
           adj[cycle[2]][cycle[3]] and adj[cycle[3]][cycle[4]] and \
           adj[cycle[4]][cycle[0]]:
            count += 1
    return count


def get_5cycle_sets(adj, n):
    """Get all 5-vertex sets with their directed 5-cycle counts."""
    result = []
    for vs in combinations(range(n), 5):
        cnt = count_5cycles(adj, vs)
        if cnt > 0:
            result.append((frozenset(vs), cnt))
    return result


def compute_alpha_from_cycles(cycle3_sets, cycle5_sets_with_counts, cycle7_count=0):
    """Compute alpha coefficients from cycle data."""
    # All cycles (3-cycles have count 1 each)
    all_cycles = [(vs, 1) for vs in cycle3_sets] + list(cycle5_sets_with_counts)
    ng = len(all_cycles)

    alpha_0 = 1
    alpha_1 = sum(cnt for _, cnt in all_cycles)

    # alpha_2: disjoint pairs
    alpha_2 = 0
    for i in range(ng):
        for j in range(i+1, ng):
            if not (all_cycles[i][0] & all_cycles[j][0]):
                alpha_2 += all_cycles[i][1] * all_cycles[j][1]

    # alpha_3: disjoint triples
    alpha_3 = 0
    for i in range(ng):
        for j in range(i+1, ng):
            if all_cycles[i][0] & all_cycles[j][0]:
                continue
            for k in range(j+1, ng):
                if not (all_cycles[k][0] & all_cycles[i][0]) and \
                   not (all_cycles[k][0] & all_cycles[j][0]):
                    alpha_3 += all_cycles[i][1] * all_cycles[j][1] * all_cycles[k][1]

    return [alpha_0, alpha_1, alpha_2, alpha_3]


def eval_ip(alpha, x):
    """Evaluate I(Omega, x) = sum alpha_k * x^k."""
    return sum(a * (x**k) for k, a in enumerate(alpha))


def main():
    print("=== Fast Ehrhart / Independence Polynomial Analysis ===\n")

    # n=5: exhaustive, 3-cycles + 5-cycles
    print("--- n=5 (exhaustive, 3+5 cycles) ---")
    n = 5
    num_edges = n * (n - 1) // 2
    max_bits = 1 << num_edges

    results_5 = {}
    for bits in range(max_bits):
        adj = tournament_adj(n, bits)
        H = count_ham_paths(adj, n)
        c3 = get_3cycle_sets(adj, n)
        c5 = get_5cycle_sets(adj, n)
        alpha = compute_alpha_from_cycles(c3, c5)
        assert eval_ip(alpha, 2) == H, f"OCF fail at n=5 bits={bits}"

        im1 = eval_ip(alpha, -1)
        i1 = eval_ip(alpha, 1)
        i3 = eval_ip(alpha, 3)

        key = (H, im1, i1, i3, tuple(alpha))
        results_5[key] = results_5.get(key, 0) + 1

    print(f"  {max_bits} tournaments, {len(results_5)} distinct tuples\n")
    print(f"  {'H':>4} {'I(-1)':>6} {'I(1)':>5} {'I(3)':>5} alpha  count")
    for (H, im1, i1, i3, alpha), cnt in sorted(results_5.items()):
        print(f"  {H:4d} {im1:6d} {i1:5d} {i3:5d} {list(alpha)}  {cnt}")

    # n=6: exhaustive, 3-cycles + 5-cycles
    print("\n\n--- n=6 (exhaustive, 3+5 cycles) ---")
    n = 6
    num_edges = n * (n - 1) // 2
    max_bits = 1 << num_edges

    h_to_im1 = {}
    h_to_i3 = {}
    alpha_report = {}

    for bits in range(max_bits):
        adj = tournament_adj(n, bits)
        H = count_ham_paths(adj, n)
        c3 = get_3cycle_sets(adj, n)
        c5 = get_5cycle_sets(adj, n)
        alpha = compute_alpha_from_cycles(c3, c5)

        ip2 = eval_ip(alpha, 2)
        assert ip2 == H, f"OCF fail at n=6 bits={bits}: H={H}, I={ip2}, alpha={alpha}"

        im1 = eval_ip(alpha, -1)
        i3 = eval_ip(alpha, 3)

        if H not in h_to_im1:
            h_to_im1[H] = set()
        h_to_im1[H].add(im1)

        if H not in h_to_i3:
            h_to_i3[H] = set()
        h_to_i3[H].add(i3)

        at = tuple(alpha)
        if H not in alpha_report:
            alpha_report[H] = set()
        alpha_report[H].add(at)

    print(f"  {max_bits} tournaments")
    print(f"\n  H -> I(Omega, -1): is it determined by H?")
    for H in sorted(h_to_im1):
        vals = sorted(h_to_im1[H])
        det = "UNIQUE" if len(vals) == 1 else f"{len(vals)} values"
        print(f"    H={H:3d}: I(-1) = {vals} ({det})")

    print(f"\n  H -> I(Omega, 3):")
    for H in sorted(h_to_i3):
        vals = sorted(h_to_i3[H])
        det = "UNIQUE" if len(vals) == 1 else f"{len(vals)} values"
        print(f"    H={H:3d}: I(3) = {vals} ({det})")

    print(f"\n  Alpha decompositions per H:")
    for H in sorted(alpha_report):
        alphas = sorted(alpha_report[H])
        print(f"    H={H:3d}: {len(alphas)} decomp: {alphas}")


if __name__ == "__main__":
    main()
