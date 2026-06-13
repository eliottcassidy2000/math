#!/usr/bin/env python3
"""
Verify the Euler zigzag identity E_T(-1) at n=9:

  E_T(-1) = A_n(-1) + sum_I 2^{parts(I)} * I(T) * A_{f+1}(-1) * (-2)^{d-f}

where d = n-1 = 8, and A_k(-1) are alternating (tangent/secant) numbers evaluated
at x = -1 via sum_j A(k,j)*(-1)^j.

Predicted coefficients for n=9:
  Constant:  A_9(-1)                                   =  7936
  t3:        2 * A_7(-1) * (-2)^2                      = -2176
  t5:        2 * A_5(-1) * (-2)^4                      =   512
  t7:        2 * A_3(-1) * (-2)^6                      =  -256
  t9:        2 * A_1(-1) * (-2)^8                      =   512
  bc33:      4 * A_5(-1) * (-2)^4                      =  1024
  bc35:      4 * A_3(-1) * (-2)^6                      =  -512
  bc37:      4 * A_1(-1) * (-2)^8                      =  1024
  a3:        8 * A_3(-1) * (-2)^6                      = -1024
"""

import random
from itertools import combinations
from collections import defaultdict
from math import comb


# ====================================================================
# Eulerian numbers and E(-1)
# ====================================================================

def eulerian_number(n, k):
    """A(n,k) = Eulerian number (0-indexed: k forward edges)."""
    return sum((-1)**j * comb(n+1, j) * (k+1-j)**n for j in range(k+2))


def euler_poly_at_minus1(n):
    """E_n(-1) = sum_{k=0}^{n-1} A(n,k) * (-1)^k."""
    return sum(eulerian_number(n, k) * (-1)**k for k in range(n))


# ====================================================================
# Tournament generation
# ====================================================================

def random_tournament(n, seed=None):
    rng = random.Random(seed)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


# ====================================================================
# Forward edge distribution via Hamiltonian path DP
# ====================================================================

def forward_edge_dist_dp(A, n):
    """Count permutations by number of forward edges using DP.
    Returns dict: {num_forward_edges: count}."""
    dp = {}
    for v in range(n):
        dp[(1 << v, v, 0)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            for fwd in range(n):
                c = dp.get((mask, v, fwd), 0)
                if c == 0:
                    continue
                for u in range(n):
                    if mask & (1 << u):
                        continue
                    new_fwd = fwd + A[v][u]
                    key = (mask | (1 << u), u, new_fwd)
                    dp[key] = dp.get(key, 0) + c
    full = (1 << n) - 1
    dist = defaultdict(int)
    for v in range(n):
        for fwd in range(n):
            dist[fwd] += dp.get((full, v, fwd), 0)
    return dict(dist)


def E_at_minus1(dist, n):
    """Evaluate E_T(-1) = sum_k a_k * (-1)^k from the distribution."""
    return sum(count * (-1)**k for k, count in dist.items())


# ====================================================================
# Cycle counting invariants
# ====================================================================

def count_t3(A, n):
    """Count directed 3-cycles."""
    total = 0
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]:
            total += 1
        if A[i][k] and A[k][j] and A[j][i]:
            total += 1
    return total


def count_directed_cycles(A, n, cl):
    """Count directed cycles of length cl using Hamiltonian path DP on subsets."""
    if n < cl:
        return 0
    total = 0
    for verts in combinations(range(n), cl):
        sub = [[A[verts[i]][verts[j]] for j in range(cl)] for i in range(cl)]
        dp = [[0]*cl for _ in range(1 << cl)]
        dp[1][0] = 1  # start at vertex 0 of the subset
        for m in range(1, 1 << cl):
            for v in range(cl):
                if not (m & (1 << v)) or dp[m][v] == 0:
                    continue
                for u in range(cl):
                    if m & (1 << u):
                        continue
                    if sub[v][u]:
                        dp[m | (1 << u)][u] += dp[m][v]
        full = (1 << cl) - 1
        total += sum(dp[full][v] for v in range(1, cl) if sub[v][0])
    return total


def get_3cycles(A, n):
    """Return list of (frozenset_of_vertices, direction_count) for all 3-cycle vertex sets."""
    cycles = []
    for t in combinations(range(n), 3):
        ct = 0
        if A[t[0]][t[1]] and A[t[1]][t[2]] and A[t[2]][t[0]]:
            ct += 1
        if A[t[0]][t[2]] and A[t[2]][t[1]] and A[t[1]][t[0]]:
            ct += 1
        if ct > 0:
            cycles.append((frozenset(t), ct))
    return cycles


def get_kcycles(A, n, cl):
    """Return list of (frozenset_of_vertices, cycle_count) for all cl-cycle vertex sets."""
    cycles = []
    for verts in combinations(range(n), cl):
        sub = [[A[verts[i]][verts[j]] for j in range(cl)] for i in range(cl)]
        dp = [[0]*cl for _ in range(1 << cl)]
        dp[1][0] = 1
        for m in range(1, 1 << cl):
            for v in range(cl):
                if not (m & (1 << v)) or dp[m][v] == 0:
                    continue
                for u in range(cl):
                    if m & (1 << u):
                        continue
                    if sub[v][u]:
                        dp[m | (1 << u)][u] += dp[m][v]
        full = (1 << cl) - 1
        ct = sum(dp[full][v] for v in range(1, cl) if sub[v][0])
        if ct > 0:
            cycles.append((frozenset(verts), ct))
    return cycles


def count_bc33(A, n):
    """Count vertex-disjoint pairs of 3-cycles (weighted by product of cycle counts)."""
    cyc3 = get_3cycles(A, n)
    total = 0
    for i in range(len(cyc3)):
        for j in range(i+1, len(cyc3)):
            if cyc3[i][0].isdisjoint(cyc3[j][0]):
                total += cyc3[i][1] * cyc3[j][1]
    return total


def count_bc35(A, n):
    """Count vertex-disjoint (3-cycle, 5-cycle) pairs, weighted by cycle counts."""
    cyc3 = get_3cycles(A, n)
    cyc5 = get_kcycles(A, n, 5)
    total = 0
    for c3_verts, c3_ct in cyc3:
        for c5_verts, c5_ct in cyc5:
            if c3_verts.isdisjoint(c5_verts):
                total += c3_ct * c5_ct
    return total


def count_bc37(A, n):
    """Count vertex-disjoint (3-cycle, 7-cycle) pairs, weighted by cycle counts."""
    cyc3 = get_3cycles(A, n)
    cyc7 = get_kcycles(A, n, 7)
    total = 0
    for c3_verts, c3_ct in cyc3:
        for c7_verts, c7_ct in cyc7:
            if c3_verts.isdisjoint(c7_verts):
                total += c3_ct * c7_ct
    return total


def count_a3(A, n):
    """Count vertex-disjoint triples of 3-cycles, weighted by product of cycle counts."""
    cyc3 = get_3cycles(A, n)
    total = 0
    nc = len(cyc3)
    for i in range(nc):
        for j in range(i+1, nc):
            if not cyc3[i][0].isdisjoint(cyc3[j][0]):
                continue
            for k in range(j+1, nc):
                if cyc3[k][0].isdisjoint(cyc3[i][0]) and cyc3[k][0].isdisjoint(cyc3[j][0]):
                    total += cyc3[i][1] * cyc3[j][1] * cyc3[k][1]
    return total


# ====================================================================
# Main verification
# ====================================================================

if __name__ == "__main__":
    n = 9
    d = n - 1  # = 8

    # Precompute A_k(-1) values
    A1 = euler_poly_at_minus1(1)
    A3 = euler_poly_at_minus1(3)
    A5 = euler_poly_at_minus1(5)
    A7 = euler_poly_at_minus1(7)
    A9 = euler_poly_at_minus1(9)

    print(f"Euler zigzag values A_k(-1):")
    print(f"  A_1(-1) = {A1}")
    print(f"  A_3(-1) = {A3}")
    print(f"  A_5(-1) = {A5}")
    print(f"  A_7(-1) = {A7}")
    print(f"  A_9(-1) = {A9}")

    # Predicted formula coefficients
    coeff_const = A9
    coeff_t3  = 2 * A7 * (-2)**2
    coeff_t5  = 2 * A5 * (-2)**4
    coeff_t7  = 2 * A3 * (-2)**6
    coeff_t9  = 2 * A1 * (-2)**8
    coeff_bc33 = 4 * A5 * (-2)**4
    coeff_bc35 = 4 * A3 * (-2)**6
    coeff_bc37 = 4 * A1 * (-2)**8
    coeff_a3  = 8 * A3 * (-2)**6

    print(f"\nPredicted coefficients:")
    print(f"  constant = {coeff_const}")
    print(f"  t3       = {coeff_t3}")
    print(f"  t5       = {coeff_t5}")
    print(f"  t7       = {coeff_t7}")
    print(f"  t9       = {coeff_t9}")
    print(f"  bc33     = {coeff_bc33}")
    print(f"  bc35     = {coeff_bc35}")
    print(f"  bc37     = {coeff_bc37}")
    print(f"  a3       = {coeff_a3}")

    print(f"\n{'='*70}")
    print(f"VERIFICATION: E_T(-1) formula on random n={n} tournaments")
    print(f"{'='*70}")

    all_pass = True
    num_trials = 3

    for trial in range(num_trials):
        seed = 1000 + trial * 137
        print(f"\nTrial {trial} (seed={seed}):")
        A = random_tournament(n, seed)

        # Step 1: Compute E_T(-1) directly via DP
        print("  Computing forward edge distribution...", flush=True)
        dist = forward_edge_dist_dp(A, n)
        actual = E_at_minus1(dist, n)
        print(f"  E_T(-1) [actual] = {actual}")

        # Step 2: Count all invariants
        print("  Counting invariants...", flush=True)
        t3  = count_t3(A, n)
        t5  = count_directed_cycles(A, n, 5)
        t7  = count_directed_cycles(A, n, 7)
        t9  = count_directed_cycles(A, n, 9)
        bc33 = count_bc33(A, n)
        bc35 = count_bc35(A, n)
        bc37 = count_bc37(A, n)
        a3  = count_a3(A, n)

        print(f"  t3={t3}, t5={t5}, t7={t7}, t9={t9}")
        print(f"  bc33={bc33}, bc35={bc35}, bc37={bc37}, a3={a3}")

        # Step 3: Evaluate formula
        predicted = (coeff_const
                     + coeff_t3 * t3
                     + coeff_t5 * t5
                     + coeff_t7 * t7
                     + coeff_t9 * t9
                     + coeff_bc33 * bc33
                     + coeff_bc35 * bc35
                     + coeff_bc37 * bc37
                     + coeff_a3 * a3)

        print(f"  E_T(-1) [formula] = {predicted}")

        if predicted == actual:
            print(f"  ==> MATCH")
        else:
            print(f"  ==> MISMATCH (diff = {actual - predicted})")
            all_pass = False

    print(f"\n{'='*70}")
    if all_pass:
        print(f"RESULT: ALL {num_trials} TRIALS PASS -- formula verified at n=9")
    else:
        print(f"RESULT: SOME TRIALS FAILED")
    print(f"{'='*70}")
