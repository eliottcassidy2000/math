#!/usr/bin/env python3
"""
Test Open Problem 9 conjecture: alpha_1(T) = |C_3(T)| (mod 2) for all tournaments T.

alpha_1(T) = number of independent pairs (vertex-disjoint pairs) of odd cycles
              in the conflict graph Omega(T).
|C_3(T)|   = number of directed 3-cycles in tournament T.

Also computes H(T) = number of Hamiltonian paths via bitmask DP,
and examines H(T) mod 4, mod 8, and score-sequence patterns.

Exhaustive for n=3,4,5,6; sampled for n=7.
"""

import sys
from itertools import combinations
from collections import defaultdict
import random

def generate_all_tournaments(n):
    """Generate all tournaments on n vertices as adjacency sets."""
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    num_edges = len(edges)
    for mask in range(1 << num_edges):
        adj = [set() for _ in range(n)]
        for k, (i, j) in enumerate(edges):
            if mask & (1 << k):
                adj[i].add(j)
            else:
                adj[j].add(i)
        yield adj

def count_ham_paths_bitmask(adj, n):
    """Count Hamiltonian paths using bitmask DP. dp[mask][v] = # paths ending at v visiting mask."""
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    full = (1 << n) - 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[mask][v] == 0:
                continue
            for u in adj[v]:
                if not (mask & (1 << u)):
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[full])

def find_all_odd_cycles(adj, n):
    """Find all directed odd cycles as frozensets of vertices.
    Enumerate all subsets of size 3,5,7,... and check if they form a directed cycle."""
    cycles = set()
    for size in range(3, n+1, 2):  # odd sizes only
        for subset in combinations(range(n), size):
            # Try all cyclic orderings of this subset
            # A directed cycle on subset S exists if there's a permutation
            # (v0, v1, ..., v_{k-1}) with edges v0->v1, v1->v2, ..., v_{k-1}->v0
            # We use DFS to find all Hamiltonian cycles in the induced subgraph
            s_list = list(subset)
            s_set = frozenset(subset)
            # Find directed Hamiltonian cycles in induced tournament on subset
            # Fix first vertex to avoid counting rotations multiple times
            first = s_list[0]
            rest = s_list[1:]
            # DFS
            found = _has_directed_ham_cycle(adj, first, rest, s_set)
            if found:
                cycles.add(s_set)
    return cycles

def _has_directed_ham_cycle(adj, first, rest, s_set):
    """Check if the induced tournament on s_set has a directed Hamiltonian cycle.
    Fix first vertex, try all permutations of rest via DFS."""
    # DFS: build path starting from first, must use all vertices in s_set, last must connect to first
    stack = [(first, frozenset([first]))]
    target_len = len(s_set)
    while stack:
        v, visited = stack.pop()
        if len(visited) == target_len:
            if first in adj[v]:
                return True
            continue
        for u in adj[v]:
            if u in s_set and u not in visited:
                stack.append((u, visited | {u}))
    return False

def count_3cycles(adj, n):
    """Count directed 3-cycles in tournament."""
    count = 0
    for triple in combinations(range(n), 3):
        i, j, k = triple
        # Check both orientations of 3-cycle
        if j in adj[i] and k in adj[j] and i in adj[k]:
            count += 1
        elif k in adj[i] and j in adj[k] and i in adj[j]:
            count += 1
    return count

def build_conflict_graph(cycles):
    """Build conflict graph: two odd cycles are adjacent iff they share a vertex."""
    cycle_list = list(cycles)
    m = len(cycle_list)
    conflict_adj = [set() for _ in range(m)]
    for i in range(m):
        for j in range(i+1, m):
            if cycle_list[i] & cycle_list[j]:  # share vertex
                conflict_adj[i].add(j)
                conflict_adj[j].add(i)
    return cycle_list, conflict_adj

def count_independent_pairs(conflict_adj, m):
    """Count independent sets of size 2 in the conflict graph = pairs of non-adjacent vertices."""
    count = 0
    for i in range(m):
        for j in range(i+1, m):
            if j not in conflict_adj[i]:
                count += 1
    return count

def independence_polynomial_at_2(cycle_list, conflict_adj, m):
    """Compute I(Omega, 2) = sum over all independent sets S of 2^|S|.
    This equals H(T) the Hamiltonian path count (by the OCF).
    We compute it to cross-check."""
    # Enumerate all independent sets
    # For small m this is feasible
    if m > 25:
        return None  # too large
    total = 0
    # Use inclusion via subset enumeration
    # Independent set = subset with no two adjacent
    for mask in range(1 << m):
        subset = [i for i in range(m) if mask & (1 << i)]
        # Check independence
        indep = True
        for idx, i in enumerate(subset):
            for j in subset[idx+1:]:
                if j in conflict_adj[i]:
                    indep = False
                    break
            if not indep:
                break
        if indep:
            total += (1 << len(subset))  # 2^|S|
    return total

def score_sequence(adj, n):
    """Return sorted score sequence."""
    scores = sorted([len(adj[v]) for v in range(n)])
    return tuple(scores)

def random_tournament(n):
    """Generate a random tournament on n vertices."""
    adj = [set() for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                adj[i].add(j)
            else:
                adj[j].add(i)
    return adj

def analyze_tournament(adj, n, verbose=False):
    """Analyze a single tournament. Returns dict of results."""
    H = count_ham_paths_bitmask(adj, n)
    c3 = count_3cycles(adj, n)
    odd_cycles = find_all_odd_cycles(adj, n)
    alpha_0 = len(odd_cycles)
    cycle_list, conflict_adj = build_conflict_graph(odd_cycles)
    m = len(cycle_list)
    alpha_1 = count_independent_pairs(conflict_adj, m)

    # Cross-check: H(T) should equal I(Omega, 2) if OCF holds
    if m <= 20:
        I2 = independence_polynomial_at_2(cycle_list, conflict_adj, m)
    else:
        I2 = None

    scores = score_sequence(adj, n)

    conjecture_holds = (alpha_1 % 2) == (c3 % 2)

    result = {
        'H': H,
        'c3': c3,
        'alpha_0': alpha_0,
        'alpha_1': alpha_1,
        'I2': I2,
        'scores': scores,
        'conjecture_holds': conjecture_holds,
        'H_mod4': H % 4,
        'H_mod8': H % 8,
    }

    if verbose and not conjecture_holds:
        print(f"  COUNTEREXAMPLE: scores={scores}, H={H}, c3={c3}, alpha_0={alpha_0}, alpha_1={alpha_1}")

    return result

def main():
    random.seed(42)

    # Exhaustive tests
    for n in range(3, 7):
        print(f"\n{'='*60}")
        print(f"n = {n}: Exhaustive test over all tournaments")
        print(f"{'='*60}")

        total = 0
        counterexamples = 0
        ocf_mismatches = 0

        # Track distributions
        h_mod4_dist = defaultdict(int)
        h_mod8_dist = defaultdict(int)
        score_stats = defaultdict(lambda: {'count': 0, 'conjecture_holds': 0,
                                            'h_mod4': defaultdict(int),
                                            'h_mod8': defaultdict(int),
                                            'alpha1_parity': defaultdict(int),
                                            'c3_parity': defaultdict(int)})

        for adj in generate_all_tournaments(n):
            total += 1
            r = analyze_tournament(adj, n, verbose=True)

            if not r['conjecture_holds']:
                counterexamples += 1

            if r['I2'] is not None and r['I2'] != r['H']:
                ocf_mismatches += 1

            h_mod4_dist[r['H_mod4']] += 1
            h_mod8_dist[r['H_mod8']] += 1

            ss = r['scores']
            score_stats[ss]['count'] += 1
            if r['conjecture_holds']:
                score_stats[ss]['conjecture_holds'] += 1
            score_stats[ss]['h_mod4'][r['H_mod4']] += 1
            score_stats[ss]['h_mod8'][r['H_mod8']] += 1
            score_stats[ss]['alpha1_parity'][r['alpha_1'] % 2] += 1
            score_stats[ss]['c3_parity'][r['c3'] % 2] += 1

        print(f"\nTotal tournaments: {total}")
        print(f"Counterexamples to alpha_1 = c3 (mod 2): {counterexamples}")
        print(f"OCF mismatches (I(Omega,2) != H): {ocf_mismatches}")
        print(f"\nH mod 4 distribution: {dict(sorted(h_mod4_dist.items()))}")
        print(f"H mod 8 distribution: {dict(sorted(h_mod8_dist.items()))}")

        print(f"\nBy score sequence:")
        for ss in sorted(score_stats.keys()):
            s = score_stats[ss]
            print(f"  {ss}: count={s['count']}, conjecture_holds={s['conjecture_holds']}/{s['count']}, "
                  f"H mod 4={dict(s['h_mod4'])}, H mod 8={dict(s['h_mod8'])}, "
                  f"alpha1 par={dict(s['alpha1_parity'])}, c3 par={dict(s['c3_parity'])}")

    # Sampled test for n=7
    print(f"\n{'='*60}")
    print(f"n = 7: Sampled test (1000 random tournaments)")
    print(f"{'='*60}")

    n = 7
    total = 0
    counterexamples = 0
    ocf_mismatches = 0
    h_mod4_dist = defaultdict(int)
    h_mod8_dist = defaultdict(int)

    num_samples = 1000
    for _ in range(num_samples):
        adj = random_tournament(n)
        total += 1
        r = analyze_tournament(adj, n, verbose=True)

        if not r['conjecture_holds']:
            counterexamples += 1

        if r['I2'] is not None and r['I2'] != r['H']:
            ocf_mismatches += 1

        h_mod4_dist[r['H_mod4']] += 1
        h_mod8_dist[r['H_mod8']] += 1

    print(f"\nTotal sampled: {total}")
    print(f"Counterexamples to alpha_1 = c3 (mod 2): {counterexamples}")
    print(f"OCF mismatches: {ocf_mismatches}")
    print(f"H mod 4 distribution: {dict(sorted(h_mod4_dist.items()))}")
    print(f"H mod 8 distribution: {dict(sorted(h_mod8_dist.items()))}")

    # Additional analysis: look at H mod 8 vs alpha_1 relationship
    print(f"\n{'='*60}")
    print(f"Detailed mod-8 analysis for n=5 (exhaustive)")
    print(f"{'='*60}")

    n = 5
    for adj in generate_all_tournaments(n):
        r = analyze_tournament(adj, n)
        # H = 1 + 2*alpha_0 + 4*alpha_1 + ... (mod 8)
        # So (H - 1) / 2 mod 2 = alpha_0 mod 2
        # And (H - 1 - 2*alpha_0) / 4 mod 2 = alpha_1 mod 2
        H = r['H']
        a0 = r['alpha_0']
        a1 = r['alpha_1']
        c3 = r['c3']

        # Check: H mod 2 should always be 1 (Redei's theorem: H is odd)
        assert H % 2 == 1, f"H={H} is even!"

        # Check: (H-1)/2 mod 2 should equal alpha_0 mod 2
        h_alpha0_check = ((H - 1) // 2) % 2 == a0 % 2

        # Check: extract alpha_1 from H
        # H = I(Omega, 2) = sum_{k>=0} alpha_k * 2^k  (wait, not exactly)
        # Actually I(Omega, x) = sum_{k>=0} alpha_k * x^k where alpha_k = # indep sets of size k
        # So I(Omega, 2) = alpha_0_is + alpha_1_is*2 + alpha_2_is*4 + ...
        # where alpha_k_is is the count of independent sets of size k
        # Note: alpha_0_is = 1 (empty set), alpha_1_is = number of odd cycles, etc.

        # So H = 1 + (# odd cycles)*2 + (# indep pairs)*4 + ...
        # H mod 2 = 1 (always)
        # (H-1)/2 mod 2 = (# odd cycles) mod 2
        # ((H-1)/2 - # odd cycles) / 2 mod 2 = ... wait
        # H = 1 + 2*a0 + 4*a1 + 8*a2 + ...
        # H mod 4 = 1 + 2*(a0 mod 2)  -- so H mod 4 is 1 or 3
        # H mod 8 = 1 + 2*a0 + 4*(a1 mod 2)  -- but a0 can be > 1
        # Actually H mod 8 = (1 + 2*a0 + 4*a1) mod 8

        predicted_H_mod4 = (1 + 2 * a0) % 4
        predicted_H_mod8 = (1 + 2 * a0 + 4 * a1) % 8

        scores = score_sequence(adj, n)
        if H % 4 != predicted_H_mod4 or H % 8 != predicted_H_mod8:
            print(f"  MISMATCH: scores={scores}, H={H}, a0={a0}, a1={a1}, "
                  f"H%4={H%4} vs pred={predicted_H_mod4}, H%8={H%8} vs pred={predicted_H_mod8}")

    print("H mod 8 = 1 + 2*alpha_0 + 4*alpha_1 (mod 8) check complete for n=5")

    # Final summary: for n=5, print each tournament's key data
    print(f"\n{'='*60}")
    print(f"Full data table for n=4")
    print(f"{'='*60}")
    print(f"{'Scores':<20} {'H':>4} {'c3':>4} {'a0':>4} {'a1':>4} {'H%4':>4} {'H%8':>4} {'a1%2':>5} {'c3%2':>5} {'match':>6}")

    n = 4
    for adj in generate_all_tournaments(n):
        r = analyze_tournament(adj, n)
        scores = score_sequence(adj, n)
        print(f"{str(scores):<20} {r['H']:>4} {r['c3']:>4} {r['alpha_0']:>4} {r['alpha_1']:>4} "
              f"{r['H_mod4']:>4} {r['H_mod8']:>4} {r['alpha_1']%2:>5} {r['c3']%2:>5} "
              f"{'YES' if r['conjecture_holds'] else 'NO':>6}")

    print(f"\n{'='*60}")
    print(f"Full data table for n=5")
    print(f"{'='*60}")
    print(f"{'Scores':<20} {'H':>4} {'c3':>4} {'a0':>4} {'a1':>4} {'H%4':>4} {'H%8':>4} {'a1%2':>5} {'c3%2':>5} {'match':>6}")

    n = 5
    for adj in generate_all_tournaments(n):
        r = analyze_tournament(adj, n)
        scores = score_sequence(adj, n)
        print(f"{str(scores):<20} {r['H']:>4} {r['c3']:>4} {r['alpha_0']:>4} {r['alpha_1']:>4} "
              f"{r['H_mod4']:>4} {r['H_mod8']:>4} {r['alpha_1']%2:>5} {r['c3']%2:>5} "
              f"{'YES' if r['conjecture_holds'] else 'NO':>6}")

if __name__ == '__main__':
    main()
