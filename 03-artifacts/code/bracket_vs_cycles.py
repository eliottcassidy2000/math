#!/usr/bin/env python3
"""
Connect the bracket structure B(u,w) to the cycle structure of Omega(T).

Key observation: B(u,w) encodes the "cross-interaction" of vertices u,w
with the flipped arc (i,j). The nonzero bracket pairs are:
  (M-, M-): +1   → both u,w form 3-cycles with i,j
  (M-, Z0): +1   → u forms a 3-cycle, w is beaten by both
  (M+, M+): -1   → both u,w would form 3-cycles in T'
  (M+, Z0): -1   → w forms a 3-cycle in T', u is beaten by both
  (Z1, M-): +1   → u beats both, w forms a 3-cycle
  (Z1, M+): -1   → u beats both, w forms a 3-cycle in T'

Key test: does the SUM Σ_S Σ_{u,w} B(u,w)*h_end*h_start decompose
into cycle-counting terms when we group by the type structure?

The (M-,M-) contribution involves pairs that BOTH form 3-cycles with {i,j}.
Vertex-disjoint 3-cycle PAIRS contribute to alpha_2 in Omega.
So the (M-, M-) bracket sum should relate to Delta(alpha_2).

Similarly, (M-, Z0) involves one cycle vertex and one non-cycle vertex.
This should relate to Delta(alpha_1) or the H(B_x) terms.

Instance: opus-2026-03-05-S3
"""

import sys
sys.path.insert(0, '.')
from tournament_lib import (
    hamiltonian_path_count, find_odd_cycles, independence_poly_at_fast
)
from itertools import permutations, combinations
import random


def h_end_sub(T, verts, target):
    """#Ham paths through verts ending at target (DP)."""
    vlist = sorted(verts)
    m = len(vlist)
    idx = {v: k for k, v in enumerate(vlist)}
    if target not in idx:
        return 0
    if m == 1:
        return 1
    dp = {}
    for v in range(m):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << m):
        for v in range(m):
            if not ((mask >> v) & 1) or (mask, v) not in dp:
                continue
            for u in range(m):
                if (mask >> u) & 1:
                    continue
                if T[vlist[v]][vlist[u]]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + dp[(mask, v)]
    full = (1 << m) - 1
    return dp.get((full, idx[target]), 0)


def h_start_sub(T, verts, source):
    """#Ham paths through verts starting at source (DP)."""
    vlist = sorted(verts)
    m = len(vlist)
    idx = {v: k for k, v in enumerate(vlist)}
    if source not in idx:
        return 0
    if m == 1:
        return 1
    dp = {(1 << idx[source], idx[source]): 1}
    for mask in range(1, 1 << m):
        for v in range(m):
            if not ((mask >> v) & 1) or (mask, v) not in dp:
                continue
            for u in range(m):
                if (mask >> u) & 1:
                    continue
                if T[vlist[v]][vlist[u]]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + dp[(mask, v)]
    full = (1 << m) - 1
    return sum(dp.get((full, v), 0) for v in range(m))


def decompose_delta_H_by_bracket_type(T, i, j):
    """Decompose Delta_H into contributions from each bracket type pair."""
    n = len(T)
    V0 = [v for v in range(n) if v != i and v != j]
    m = len(V0)
    v0_idx = {v: k for k, v in enumerate(V0)}

    # Classify
    s = {}
    typ = {}
    for x in V0:
        s[x] = 1 - T[x][i] - T[j][x]
        if s[x] == -1:
            typ[x] = 'M-'
        elif s[x] == 1:
            typ[x] = 'M+'
        elif T[x][i] == 1:
            typ[x] = 'Z1'
        else:
            typ[x] = 'Z0'

    # Bracket table (nonzero entries only)
    bracket_val = {
        ('M-', 'M-'): 1, ('M-', 'Z0'): 1,
        ('M+', 'M+'): -1, ('M+', 'Z0'): -1,
        ('Z1', 'M-'): 1, ('Z1', 'M+'): -1,
    }

    # Compute the convolution sum, decomposed by bracket type
    contributions = {}  # (ut, wt) -> total contribution
    total = 0

    for mask in range(1 << m):
        S = [V0[k] for k in range(m) if (mask >> k) & 1]
        R = [V0[k] for k in range(m) if not ((mask >> k) & 1)]

        for u in S:
            for w in R:
                key = (typ[u], typ[w])
                bval = bracket_val.get(key, 0)
                if bval == 0:
                    continue

                he = h_end_sub(T, S, u) if len(S) > 0 else 0
                hs = h_start_sub(T, R, w) if len(R) > 0 else 0
                contrib = bval * he * hs
                contributions[key] = contributions.get(key, 0) + contrib
                total += contrib

    # Also handle boundary: S=empty (f_i(empty)=1) and R=empty (g_j(empty)=1)
    # When S=empty: u doesn't exist, only g_j(R=V0) - g_i(R=V0) matters
    # When R=empty: w doesn't exist, only f_i(S=V0) - f_j(S=V0) matters
    # These contribute to total but don't have (u,w) pairs.

    # f_i(empty)*g_j(V0) - f_j(empty)*g_i(V0) = g_j(V0) - g_i(V0)
    gj_V0 = sum(h_start_sub(T, V0, w) for w in V0 if T[j][w])
    gi_V0 = sum(h_start_sub(T, V0, w) for w in V0 if T[i][w])
    boundary_empty_S = gj_V0 - gi_V0

    # f_i(V0)*g_j(empty) - f_j(V0)*g_i(empty) = f_i(V0) - f_j(V0)
    fi_V0 = sum(h_end_sub(T, V0, u) for u in V0 if T[u][i])
    fj_V0 = sum(h_end_sub(T, V0, u) for u in V0 if T[u][j])
    boundary_empty_R = fi_V0 - fj_V0

    total += boundary_empty_S + boundary_empty_R

    # Actual Delta_H
    Tp = [row[:] for row in T]
    Tp[i][j] = 0
    Tp[j][i] = 1
    H_T = hamiltonian_path_count(T)
    H_Tp = hamiltonian_path_count(Tp)
    delta_H_actual = H_T - H_Tp

    # Delta_I
    cycles = find_odd_cycles(T)
    I_T = independence_poly_at_fast(cycles, 2)
    cycles_p = find_odd_cycles(Tp)
    I_Tp = independence_poly_at_fast(cycles_p, 2)
    delta_I = I_T - I_Tp

    # H(B_x) for s_x != 0
    H_Bx = {}
    for x in V0:
        if s[x] != 0:
            B_x = [v for v in V0 if v != x]
            H_Bx[x] = hamiltonian_path_count(
                [[T[a][b] for b in B_x] for a in B_x]
            )  # Wrong: need proper sub-matrix
            # Fix: build sub-tournament
            idx_map = {v: k for k, v in enumerate(B_x)}
            mm = len(B_x)
            sub = [[0]*mm for _ in range(mm)]
            for a in B_x:
                for b in B_x:
                    if a != b:
                        sub[idx_map[a]][idx_map[b]] = T[a][b]
            H_Bx[x] = hamiltonian_path_count(sub)

    print(f"\nn={n}, arc ({i},{j})")
    print(f"Types: {typ}")
    print(f"Delta_H = {delta_H_actual}, Delta_I = {delta_I}")
    print(f"\nBracket-type contributions:")
    for key in sorted(contributions.keys()):
        print(f"  {key[0]:>2} x {key[1]:>2}: {contributions[key]:>6}")
    print(f"  Boundary (S=empty): {boundary_empty_S}")
    print(f"  Boundary (R=empty): {boundary_empty_R}")
    print(f"  TOTAL: {total} (should be {delta_H_actual})")

    # THM-013 formula terms
    sum_s = sum(s[x] for x in V0)
    sum_sH = sum(s[x] * H_Bx.get(x, 0) for x in V0 if x in H_Bx)
    print(f"\nTHM-013 terms:")
    print(f"  -2*sum(s_x) = {-2*sum_s} (3-cycle contribution)")
    print(f"  -2*sum(s_x*H(B_x)) = {-2*sum_sH}")
    print(f"  delta_I - (-2*sum_sH) = {delta_I - (-2*sum_sH)} (higher cycle corrections)")

    # Key test: do (M-,M-) and (M+,M+) contributions relate to VD 3-cycle pairs?
    # M- vertices form 3-cycles with {i,j} in T.
    # M+ vertices form 3-cycles with {i,j} in T'.
    # VD pairs among M- contribute to alpha_2 of Omega(T).
    M_minus = [x for x in V0 if typ[x] == 'M-']
    M_plus = [x for x in V0 if typ[x] == 'M+']

    # Count VD pairs among M- (all pairs are VD since cycles are {i,j,x})
    # Wait: two 3-cycles {i,j,x1} and {i,j,x2} SHARE vertices i,j,
    # so they are NOT vertex-disjoint! They're NOT independent in Omega.
    # Two 3-cycles {i,j,x1} and {i,j,x2} both contain i and j.
    # So they always share vertices i and j.
    # Therefore, no two affected 3-cycles are vertex-disjoint!
    # This means Delta(alpha_2) from 3-cycle changes is ZERO for n<=7.
    # Wait, that can't be right...

    # Actually, alpha_2 counts pairs of VD cycles in Omega, not just affected ones.
    # An affected cycle changes (appears/disappears), and it pairs with an
    # UNAFFECTED cycle in the complement.
    # The A-clique property says: at most one affected cycle in an independent set.
    # So Delta(alpha_2) = sum_{affected C} [+/- alpha_1(comp(C))]
    #                   = sum_{affected C} [+/- #{unaffected cycles in V\V(C)}]

    # For 3-cycle C = {i,j,x}: comp(C) = V\{i,j,x} = B_x.
    # alpha_1(comp(C)) = #{odd cycles in T[B_x]}.
    # By OCF for B_x: H(B_x) = 1 + 2*alpha_1(B_x) + 4*alpha_2(B_x) + ...
    # So alpha_1(B_x) = (H(B_x) - 1 - 4*alpha_2(B_x) - ...)/2

    # For small B_x (|B_x| <= 5), alpha_2 = 0, so alpha_1(B_x) = (H(B_x) - 1)/2.

    for x in M_minus:
        B_x = [v for v in V0 if v != x]
        cycles_Bx = find_odd_cycles([[T[a][b] if a != b else 0 for b in B_x] for a in B_x])
        # Fix: proper sub-tournament
        idx_m = {v: k for k, v in enumerate(B_x)}
        mm = len(B_x)
        sub = [[0]*mm for _ in range(mm)]
        for a in B_x:
            for b in B_x:
                if a != b:
                    sub[idx_m[a]][idx_m[b]] = T[a][b]
        cycles_Bx = find_odd_cycles(sub)
        print(f"\n  M- vertex {x}: H(B_x)={H_Bx[x]}, #cycles(B_x)={len(cycles_Bx)}")

    for x in M_plus:
        B_x = [v for v in V0 if v != x]
        idx_m = {v: k for k, v in enumerate(B_x)}
        mm = len(B_x)
        sub = [[0]*mm for _ in range(mm)]
        for a in B_x:
            for b in B_x:
                if a != b:
                    sub[idx_m[a]][idx_m[b]] = T[a][b]
        cycles_Bx = find_odd_cycles(sub)
        print(f"  M+ vertex {x}: H(B_x)={H_Bx[x]}, #cycles(B_x)={len(cycles_Bx)}")


def main():
    rng = random.Random(42)
    for n in [5, 6, 7]:
        T = [[0]*n for _ in range(n)]
        for a in range(n):
            for b in range(a+1, n):
                if rng.random() < 0.5:
                    T[a][b] = 1
                else:
                    T[b][a] = 1
        ii, jj = 0, 1
        if T[ii][jj] == 0:
            ii, jj = jj, ii
        decompose_delta_H_by_bracket_type(T, ii, jj)


if __name__ == "__main__":
    main()
