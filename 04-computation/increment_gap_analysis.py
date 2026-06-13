#!/usr/bin/env python3
"""increment_gap_analysis.py -- Why are H=7 and H=21 the only forbidden values?

Session: kind-pasteur-2026-03-20-S2

The deletion increment Delta(T, v) = H(T) - H(T\v) is always even (Claim A).
For H=k to be forbidden, EVERY tournament T' with H(T') < k must have
the increment Delta = k - H(T') NOT achievable by any vertex addition.

KEY INSIGHT from Part B: H=7 is blocked because:
  - From H'=1: Delta=6 not in achievable deltas {0,2,4,8,10,14,16,18,22}
  - From H'=3: Delta=4 not in {0,2,6,10,12,14,20,22,24,26,30}
  - From H'=5: Delta=2 not in {0,4,6,8,10,12,14,18,20,24,26,28,32}

WHY are these specific deltas missing? The OCF gives:
  Delta(T,v) = 2 * sum_{C ni v} mu(C)
where C ranges over directed odd cycles through v, and mu(C) is a
positive integer (independence polynomial weight).

The minimum nonzero Delta is 2 (from one 3-cycle with mu=1).
The missing deltas correspond to cycle-count combinations that are
structurally impossible.

This script investigates the delta gap structure at n=5,6,7.
"""

from itertools import combinations, permutations
from collections import Counter, defaultdict
from math import comb, factorial

def adj_from_bits(bits, n):
    adj = [[0]*n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> k) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            k += 1
    return adj

def held_karp_H(adj, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))


def analyze_increment_gaps():
    print("=" * 70)
    print("INCREMENT GAP ANALYSIS — WHY 7 AND 21 ARE FORBIDDEN")
    print("=" * 70)

    for n in [5, 6, 7]:
        print(f"\n{'='*60}")
        print(f"  n = {n}")
        print(f"{'='*60}")

        m = n * (n - 1) // 2

        # For n=7, we can't enumerate all 2^21 = 2M tournaments easily
        # but we CAN enumerate all ways to extend an (n-1)-tournament
        # by adding vertex n-1 with n-1 new arcs (2^{n-1} extensions per T')

        n_prev = n - 1
        m_prev = n_prev * (n_prev - 1) // 2

        # Compute H for all (n-1)-tournaments
        H_prev_map = {}
        for bits in range(1 << m_prev):
            adj = adj_from_bits(bits, n_prev)
            H_prev_map[bits] = held_karp_H(adj, n_prev)

        # For each (n-1)-tournament, compute all possible deltas
        # when adding vertex n-1 with each of 2^{n-1} arc orientations
        delta_by_Hprev = defaultdict(set)
        H_achievable = set()

        for bits_prev in range(1 << m_prev):
            H_prev = H_prev_map[bits_prev]
            adj_prev = adj_from_bits(bits_prev, n_prev)

            # Try all 2^{n-1} orientations for the new vertex
            for new_bits in range(1 << (n - 1)):
                # Build full tournament
                adj = [row[:] + [0] for row in adj_prev]
                adj.append([0] * n)
                for i in range(n - 1):
                    if (new_bits >> i) & 1:
                        adj[n-1][i] = 1  # new vertex -> i
                    else:
                        adj[i][n-1] = 1  # i -> new vertex

                H = held_karp_H(adj, n)
                delta = H - H_prev
                delta_by_Hprev[H_prev].add(delta)
                H_achievable.add(H)

        # Report
        max_H = max(H_achievable)
        odd_achievable = sorted(h for h in H_achievable if h % 2 == 1)
        odd_expected = set(range(1, max_H + 1, 2))
        odd_missing = sorted(odd_expected - set(odd_achievable))

        print(f"\n  H values at n={n}: {len(odd_achievable)} odd values in [1, {max_H}]")
        print(f"  Missing: {odd_missing}")

        # For each missing H value, show WHY it's missing
        for target in odd_missing:
            print(f"\n  --- WHY H={target} IS FORBIDDEN ---")
            for H_prev in sorted(delta_by_Hprev.keys()):
                if H_prev >= target:
                    continue  # can only reach target from below
                needed_delta = target - H_prev
                if needed_delta % 2 != 0:
                    continue  # delta must be even
                available = sorted(delta_by_Hprev[H_prev])
                has_needed = needed_delta in available

                if not has_needed:
                    # Find the closest available deltas
                    below = [d for d in available if d < needed_delta]
                    above = [d for d in available if d > needed_delta]
                    closest_below = max(below) if below else None
                    closest_above = min(above) if above else None
                    print(f"    From H'={H_prev}: need delta={needed_delta}, "
                          f"MISSING (closest: {closest_below}<gap>{closest_above})")

        # Interesting: what is the SET of achievable deltas for each H'?
        print(f"\n  Delta gap structure:")
        for H_prev in sorted(delta_by_Hprev.keys()):
            deltas = sorted(delta_by_Hprev[H_prev])
            gaps = []
            for i in range(len(deltas) - 1):
                gap = deltas[i+1] - deltas[i]
                if gap > 2:
                    gaps.append((deltas[i], deltas[i+1], gap))
            if gaps:
                print(f"    H'={H_prev}: gaps in deltas: {gaps}")


def analyze_delta_via_ocf():
    """Connect deltas to the OCF formula."""
    print(f"\n\n{'='*70}")
    print(f"DELTA VIA OCF: Delta(T,v) = 2 * sum_C mu(C)")
    print(f"{'='*70}")

    # At n=5 (adding to n=4):
    # Delta = 2 * (number of 3-cycles through v, each with mu=1)
    # At n=4, all mu=1 (THM-008 triviality bound)
    # So Delta = 2 * c3(v) where c3(v) = number of 3-cycles containing v

    n = 5
    n_prev = 4
    m_prev = n_prev * (n_prev - 1) // 2

    print(f"\n  n={n}: Delta = 2 * c3(v) (all mu=1 at n<=5)")

    delta_by_c3v = defaultdict(set)

    for bits_prev in range(1 << m_prev):
        adj_prev = adj_from_bits(bits_prev, n_prev)
        H_prev = held_karp_H(adj_prev, n_prev)

        for new_bits in range(1 << (n - 1)):
            adj = [row[:] + [0] for row in adj_prev]
            adj.append([0] * n)
            for i in range(n - 1):
                if (new_bits >> i) & 1:
                    adj[n-1][i] = 1
                else:
                    adj[i][n-1] = 1

            H = held_karp_H(adj, n)
            delta = H - H_prev

            # Count 3-cycles through vertex n-1
            v = n - 1
            c3v = 0
            for i in range(v):
                for j in range(i+1, v):
                    if (adj[v][i] and adj[i][j] and adj[j][v]) or \
                       (adj[v][j] and adj[j][i] and adj[i][v]):
                        c3v += 1

            delta_by_c3v[c3v].add(delta)

    print(f"  c3(v) -> delta mapping:")
    for c3v in sorted(delta_by_c3v.keys()):
        deltas = sorted(delta_by_c3v[c3v])
        print(f"    c3(v)={c3v}: deltas = {deltas}, expected 2*c3(v) = {2*c3v}")
        # Check if delta = 2*c3(v) always
        # At n=5: delta = 2 * sum_{C ni v} mu(C) = 2 * c3v (since all mu=1 AND no 5-cycles through v... wait, there CAN be 5-cycles)
        # Actually at n=5, a 5-cycle through v uses all vertices. mu(5-cycle) = 1.
        # So delta = 2 * (c3v + c5v) where c5v = number of 5-cycles through v

    # Now check with 5-cycles
    print(f"\n  Including 5-cycles through v:")
    delta_check = defaultdict(list)
    for bits_prev in range(1 << m_prev):
        adj_prev = adj_from_bits(bits_prev, n_prev)
        H_prev = held_karp_H(adj_prev, n_prev)

        for new_bits in range(1 << (n - 1)):
            adj = [row[:] + [0] for row in adj_prev]
            adj.append([0] * n)
            for i in range(n - 1):
                if (new_bits >> i) & 1:
                    adj[n-1][i] = 1
                else:
                    adj[i][n-1] = 1

            H = held_karp_H(adj, n)
            delta = H - H_prev

            v = n - 1
            c3v = 0
            for i in range(v):
                for j in range(i+1, v):
                    if (adj[v][i] and adj[i][j] and adj[j][v]) or \
                       (adj[v][j] and adj[j][i] and adj[i][v]):
                        c3v += 1

            # Count 5-cycles through v (= Hamiltonian cycles through v)
            # At n=5, this is just the number of directed Hamiltonian cycles containing v
            c5v = 0
            others = [i for i in range(n) if i != v]
            for perm in permutations(others):
                path = [v] + list(perm)
                ok = True
                for idx in range(n):
                    if not adj[path[idx]][path[(idx+1) % n]]:
                        ok = False
                        break
                if ok:
                    c5v += 1

            predicted = 2 * (c3v + c5v)
            delta_check[(c3v, c5v)].append((delta, predicted, delta == predicted))

    # Verify
    mismatches = 0
    for (c3v, c5v), entries in sorted(delta_check.items()):
        all_match = all(e[2] for e in entries)
        delta_vals = set(e[0] for e in entries)
        pred_vals = set(e[1] for e in entries)
        if not all_match:
            mismatches += 1
            print(f"    c3v={c3v}, c5v={c5v}: MISMATCH! deltas={delta_vals}, predicted={pred_vals}")
        else:
            print(f"    c3v={c3v}, c5v={c5v}: delta = {list(delta_vals)[0]} = 2*({c3v}+{c5v}) EXACT")

    print(f"\n  Total mismatches: {mismatches}")
    if mismatches == 0:
        print(f"  CONFIRMED: Delta(T,v) = 2 * (c3(v) + c5(v)) at n=5")
        print(f"  This is OCF: Delta = 2 * sum_{'{C ni v}'} mu(C) with all mu=1")

    # KEY INSIGHT: H=7 requires Delta=6 from H'=1 (transitive)
    # Delta=6 means c3(v) + c5(v) = 3
    # From transitive T_4 (H=1): when we add v, c3(v) depends on how v connects
    # The maximum c3(v) when adding to transitive T_4: v must form cycles with T_4 arcs
    print(f"\n  From transitive T_4 (H'=1):")
    print(f"  Available deltas = {sorted(delta_by_c3v.keys())}")
    print(f"  For Delta=6: need c3(v)+c5(v)=3")
    print(f"  But from transitive T_4, the achievable (c3v, c5v) pairs determine")
    print(f"  which deltas are possible. If (c3v+c5v)=3 is not achievable,")
    print(f"  then H=7 is blocked from this direction.")


if __name__ == "__main__":
    analyze_increment_gaps()
    analyze_delta_via_ocf()
