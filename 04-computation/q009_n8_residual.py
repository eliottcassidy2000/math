"""
Q-009: Analyze the n=8 residual.

The formula DH = -2*sum(s_x*H(B_x)) + 2*sum(DL-CL) fails at n=8.
Residuals are multiples of 4, suggesting missing VD pair terms
involving 5-cycles (a 3-cycle and a 5-cycle can be VD at n=8).

Let's compute what the missing term is.

Author: opus-2026-03-05-S2
"""

import random
from itertools import combinations, permutations


def hamiltonian_path_count_dp(T, verts):
    n = len(verts)
    if n <= 1:
        return 1
    dp = [[0] * n for _ in range(1 << n)]
    for k in range(n):
        dp[1 << k][k] = 1
    full = (1 << n) - 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if T[verts[v]][verts[u]]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[full][v] for v in range(n))


def count_cycles_using_arc_by_length(T, n, i, j, min_length=5):
    others = [u for u in range(n) if u != i and u != j]
    counts = {}
    for L in range(min_length, n + 1, 2):
        extra_needed = L - 2
        if extra_needed > len(others):
            break
        count = 0
        for extra in combinations(others, extra_needed):
            verts = list(extra)
            nv = len(verts)
            dp = [[0] * nv for _ in range(1 << nv)]
            for k in range(nv):
                if T[j][verts[k]]:
                    dp[1 << k][k] = 1
            full = (1 << nv) - 1
            for mask in range(1, 1 << nv):
                for v in range(nv):
                    if not (mask & (1 << v)) or dp[mask][v] == 0:
                        continue
                    for u in range(nv):
                        if mask & (1 << u):
                            continue
                        if T[verts[v]][verts[u]]:
                            dp[mask | (1 << u)][u] += dp[mask][v]
            for v in range(nv):
                if dp[full][v] > 0 and T[verts[v]][i]:
                    count += dp[full][v]
        counts[L] = count
    return counts


def is_directed_cycle(T, cycle):
    """Check if cycle is a directed cycle in T."""
    for k in range(len(cycle)):
        if T[cycle[k]][cycle[(k+1) % len(cycle)]] != 1:
            return False
    return True


def count_directed_odd_cycles(T, n):
    """Return list of all directed odd cycles (as tuples of vertices)."""
    cycles = []
    for L in range(3, n + 1, 2):
        for verts in combinations(range(n), L):
            # Check all cyclic orderings (fix first vertex to avoid counting rotations)
            v0 = verts[0]
            rest = list(verts[1:])
            for perm in permutations(rest):
                cycle = (v0,) + perm
                if is_directed_cycle(T, cycle):
                    cycles.append(cycle)
                    break  # Only count one rotation per cycle... wait, this doesn't work
    # Actually, for directed cycles we need to be more careful.
    # A directed cycle on vertex set S is a cyclic permutation.
    # There are (|S|-1)! orderings, and |S| rotations of each cycle.
    # So #directed cycles on S = (#{valid cyclic perms}) / |S|.
    # Let me redo this properly.
    cycles = []
    for L in range(3, n + 1, 2):
        for verts in combinations(range(n), L):
            v0 = verts[0]
            rest = list(verts[1:])
            cycle_count = 0
            for perm in permutations(rest):
                cycle = (v0,) + perm
                if is_directed_cycle(T, cycle):
                    cycle_count += 1
            # Each distinct directed cycle appears (L-1)!/L * L = (L-1)! times
            # with v0 fixed... no. With v0 fixed, each directed cycle appears exactly once
            # (since fixing one vertex determines the rotation).
            # But there might be two directions (CW/CCW). For odd cycles in tournaments,
            # at most one direction works.
            # Actually, fixing v0 and permuting the rest: each directed cycle on these
            # vertices with v0 in it appears exactly once (the unique rotation starting at v0).
            for perm in permutations(rest):
                cycle = (v0,) + perm
                if is_directed_cycle(T, cycle):
                    cycles.append(frozenset(verts))  # Store vertex set
                    break  # We just need to know it exists, not count orderings
    # Hmm, but I want to count directed cycles properly.
    # Let me just count: for each vertex set, how many directed cycles exist?
    # A directed cycle on L vertices = a cyclic ordering where each arc goes forward.
    # Fixing vertex 0, there are (L-1)! linear orderings of the rest.
    # Each directed cycle contributes exactly 1 such ordering (the one starting at v0).
    result = []
    for L in range(3, n + 1, 2):
        for verts in combinations(range(n), L):
            v0 = verts[0]
            rest = list(verts[1:])
            for perm in permutations(rest):
                cycle = [v0] + list(perm)
                if is_directed_cycle(T, cycle):
                    result.append(tuple(cycle))
    return result


def count_vd_pairs_by_type(cycles):
    """Count VD pairs of odd cycles, grouped by (len1, len2)."""
    from collections import Counter
    type_counts = Counter()
    for a in range(len(cycles)):
        for b in range(a + 1, len(cycles)):
            sa = frozenset(cycles[a])
            sb = frozenset(cycles[b])
            if not (sa & sb):
                la, lb = len(cycles[a]), len(cycles[b])
                key = tuple(sorted([la, lb]))
                type_counts[key] += 1
    return type_counts


def flip_arc(T, i, j):
    T2 = [row[:] for row in T]
    T2[i][j] = 0
    T2[j][i] = 1
    return T2


def random_tournament(n):
    T = [[0]*n for _ in range(n)]
    for a in range(n):
        for b in range(a+1, n):
            if random.random() < 0.5:
                T[a][b] = 1
            else:
                T[b][a] = 1
    return T


if __name__ == "__main__":
    random.seed(42)
    n = 8

    for trial in range(5):
        T = random_tournament(n)
        i, j = random.sample(range(n), 2)
        if T[i][j] == 0:
            i, j = j, i

        T2 = flip_arc(T, i, j)
        others = [x for x in range(n) if x != i and x != j]
        all_verts = list(range(n))

        # Formula
        sx_hbx = 0
        for x in others:
            sx = 1 - T[x][i] - T[j][x]
            bx = [u for u in others if u != x]
            hbx = hamiltonian_path_count_dp(T, bx)
            sx_hbx += sx * hbx

        dc_T = count_cycles_using_arc_by_length(T, n, i, j, min_length=5)
        dc_T2 = count_cycles_using_arc_by_length(T2, n, j, i, min_length=5)
        cycle_term = sum(dc_T.get(L, 0) - dc_T2.get(L, 0) for L in range(5, n+1, 2))

        formula = -2 * sx_hbx + 2 * cycle_term

        # Ground truth
        ht = hamiltonian_path_count_dp(T, all_verts)
        ht2 = hamiltonian_path_count_dp(T2, all_verts)
        delta_H = ht - ht2

        residual = delta_H - formula

        # Count VD pairs
        cycles_T = count_directed_odd_cycles(T, n)
        cycles_T2 = count_directed_odd_cycles(T2, n)
        vd_T = count_vd_pairs_by_type(cycles_T)
        vd_T2 = count_vd_pairs_by_type(cycles_T2)

        print(f"\nTrial {trial}: flip {i}->{j}")
        print(f"  DH={delta_H}, formula={formula}, residual={residual}")
        print(f"  dc_T={dc_T}, dc_T2={dc_T2}")
        print(f"  VD pairs T:  {dict(vd_T)}")
        print(f"  VD pairs T': {dict(vd_T2)}")

        # Delta VD by type
        all_types = set(list(vd_T.keys()) + list(vd_T2.keys()))
        for tp in sorted(all_types):
            delta = vd_T[tp] - vd_T2[tp]
            if delta != 0:
                print(f"  Delta VD {tp}: {delta}")

        # Check: residual = 4 * delta(VD 3-5 pairs)?
        delta_35 = vd_T.get((3,5), 0) - vd_T2.get((3,5), 0)
        delta_55 = vd_T.get((5,5), 0) - vd_T2.get((5,5), 0)
        delta_37 = vd_T.get((3,7), 0) - vd_T2.get((3,7), 0)
        print(f"  residual={residual}, 4*d35={4*delta_35}, 4*d55={4*delta_55}, 4*d37={4*delta_37}")
        print(f"  4*(d35+d55+d37)={4*(delta_35+delta_55+delta_37)}")
