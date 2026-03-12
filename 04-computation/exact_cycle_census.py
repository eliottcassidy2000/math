"""
exact_cycle_census.py — Exact odd cycle counts for circulant tournaments

CRITICAL: Tr(A^ℓ)/ℓ = C_ℓ only for ℓ = 3 and ℓ = 5 in tournaments.
For ℓ ≥ 7, non-simple closed walks contribute to the trace.

This script computes EXACT directed cycle counts by DFS for p=7, 11.

Author: opus-2026-03-12-S60
"""
import sys
import time
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)


def circulant_tournaments(n):
    k = (n - 1) // 2
    pairs = [(d, n - d) for d in range(1, k + 1)]
    for bits in range(1 << k):
        S = set()
        for i in range(k):
            if bits & (1 << i):
                S.add(pairs[i][0])
            else:
                S.add(pairs[i][1])
        yield frozenset(S)


def is_paley(n, S):
    qr = set(pow(x, 2, n) for x in range(1, n))
    return set(S) == qr or set(S) == set(range(1, n)) - qr


def exact_cycle_counts(n, S):
    """Count ALL directed cycles by length using DFS from each vertex.

    For circulant: each cycle of length ℓ visits ℓ distinct vertices.
    To avoid double-counting, only count cycles where the starting vertex
    is the smallest vertex in the cycle.
    """
    S_set = set(S)
    counts = Counter()

    def dfs(start, current, visited, length):
        for s in S_set:
            nxt = (current + s) % n
            if nxt == start and length >= 3:
                counts[length] += 1
            elif nxt > start and nxt not in visited:
                visited.add(nxt)
                dfs(start, nxt, visited, length + 1)
                visited.remove(nxt)

    for v in range(n):
        dfs(v, v, {v}, 1)

    return dict(sorted(counts.items()))


def compute_indep_poly(n, S, max_alpha=5):
    """Compute independence polynomial of odd-cycle complex.

    α_k = number of collections of k pairwise vertex-disjoint odd cycles.
    """
    S_set = set(S)

    # First find all odd directed cycles
    all_cycles = []

    def dfs(start, current, visited_list, length):
        for s in S_set:
            nxt = (current + s) % n
            if nxt == start and length >= 3 and length % 2 == 1:
                all_cycles.append(frozenset(visited_list))
            elif nxt > start and nxt not in set(visited_list):
                dfs(start, nxt, visited_list + [nxt], length + 1)

    for v in range(n):
        dfs(v, v, [v], 1)

    m = len(all_cycles)
    vertex_sets = all_cycles

    # For large m, use greedy independence polynomial computation
    alpha = [0] * (max_alpha + 1)
    alpha[0] = 1

    if m <= 500:
        # Exact computation using inclusion
        def count_indep(idx, size, used):
            if size <= max_alpha:
                alpha[size] += 1
            if size >= max_alpha:
                return
            for i in range(idx, m):
                if vertex_sets[i].isdisjoint(used):
                    count_indep(i + 1, size + 1, used | vertex_sets[i])

        count_indep(0, 0, frozenset())
        # Fix: alpha[0] was set both initially and in recursion
        alpha[0] = 1  # Remove the duplicate

    else:
        # Too many cycles; approximate
        alpha[1] = m
        for i in range(m):
            for j in range(i+1, m):
                if vertex_sets[i].isdisjoint(vertex_sets[j]):
                    alpha[2] += 1

    return alpha


def circulant_eigenvalues(n, S):
    omega = np.exp(2j * np.pi / n)
    return [sum(omega ** (k * s) for s in S) for k in range(n)]


def main():
    print("EXACT CYCLE CENSUS FOR CIRCULANT TOURNAMENTS")
    print("=" * 75)

    for p in [7, 11]:
        print(f"\n{'='*75}")
        print(f"p = {p} (p mod 4 = {p % 4})")
        print(f"{'='*75}")

        tournaments = list(circulant_tournaments(p))
        results = []

        t0 = time.time()
        for S in tournaments:
            cycles = exact_cycle_counts(p, S)
            odd_cycles = {k: v for k, v in cycles.items() if k % 2 == 1}
            total_odd = sum(odd_cycles.values())

            evals = circulant_eigenvalues(p, S)
            y_vals = [evals[k].imag for k in range(1, p)]
            sum_y4 = sum(y**4 for y in y_vals)

            paley = is_paley(p, S)
            results.append({
                'S': S, 'paley': paley, 'odd_cycles': odd_cycles,
                'total_odd': total_odd, 'sum_y4': sum_y4
            })

        elapsed = time.time() - t0
        print(f"  Computed in {elapsed:.2f}s")

        # Sort by total_odd descending
        results.sort(key=lambda r: -r['total_odd'])

        # Print all
        print(f"\n  {'S':>25} {'Label':>8} {'C3':>5} {'C5':>5} {'C7':>6} {'C9':>7} {'C11':>7} {'Total':>8} {'Σy⁴':>10}")
        for r in results:
            label = "PALEY" if r['paley'] else "other"
            oc = r['odd_cycles']
            print(f"  {str(sorted(r['S'])):>25} {label:>8} "
                  f"{oc.get(3,0):>5} {oc.get(5,0):>5} {oc.get(7,0):>6} "
                  f"{oc.get(9,0):>7} {oc.get(11,0):>7} "
                  f"{r['total_odd']:>8} {r['sum_y4']:>10.2f}")

        # Key question: does Paley maximize total odd cycles?
        paley_total = max(r['total_odd'] for r in results if r['paley']) if any(r['paley'] for r in results) else 0
        max_total = results[0]['total_odd']
        print(f"\n  Paley total odd cycles: {paley_total}")
        print(f"  Max total odd cycles: {max_total}")
        print(f"  Paley is max: {paley_total == max_total}")

        # Independence polynomial for top and bottom
        print(f"\n  Independence polynomial (top 2 and bottom 2):")
        for r in [results[0], results[1], results[-2], results[-1]]:
            label = "PALEY" if r['paley'] else "other"
            alpha = compute_indep_poly(p, r['S'], max_alpha=4)
            H = sum(alpha[k] * 2**k for k in range(len(alpha)))
            print(f"    S={sorted(r['S'])} [{label}]: α=[{', '.join(str(a) for a in alpha)}] → H = {H}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
