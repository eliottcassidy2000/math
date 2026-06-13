#!/usr/bin/env python3
"""
Analyze the bracket B(u,w) = T[u][i]*T[j][w] - T[u][j]*T[i][w] structure.

Key discovery: B(u,w) depends only on the s-types and t-types of u and w:
  s_x = 1 - T[x][i] - T[j][x]
  For s=0: t_x = T[x][i] (which of the two s=0 sub-types)

The bracket table is:
  B(u,w) by type:
       | M-  M+  Z1  Z0
  -----|------------------
  M-   |  1   0   0   1
  M+   |  0  -1   0  -1
  Z1   |  1  -1   0   0
  Z0   |  0   0   0   0

where M- = {s=-1}, M+ = {s=+1}, Z1 = {s=0, t=1}, Z0 = {s=0, t=0}.

CRITICAL OBSERVATION: Z0 rows and Z1 columns are all zero!
This means:
- Z0-type vertices in u-position contribute NOTHING to Delta_H
- Z1-type vertices in w-position contribute NOTHING
- The "effective" vertices are: u from M- ∪ M+ ∪ Z1, w from M- ∪ M+ ∪ Z0

Instance: opus-2026-03-05-S3
"""

import sys
sys.path.insert(0, '.')
from tournament_lib import hamiltonian_path_count, find_odd_cycles, independence_poly_at_fast
from itertools import permutations, combinations
import random


def classify_vertex(T, x, i, j):
    """Classify x into M-, M+, Z1, Z0."""
    s = 1 - T[x][i] - T[j][x]
    if s == -1:
        return 'M-'
    elif s == 1:
        return 'M+'
    else:  # s == 0
        return 'Z1' if T[x][i] == 1 else 'Z0'


def bracket(T, u, w, i, j):
    """Compute B(u,w) = T[u][i]*T[j][w] - T[u][j]*T[i][w]."""
    return T[u][i] * T[j][w] - T[u][j] * T[i][w]


def verify_bracket_table():
    """Verify the bracket table by enumeration."""
    print("Verifying bracket table B(u,w) by type:")
    print("="*50)

    # For each type pair, compute B
    # Types defined by (p_u = T[u][i], q_u = T[j][u]):
    # M-: p=1, q=1  (s = 1-1-1 = -1)
    # M+: p=0, q=0  (s = 1-0-0 = +1)
    # Z1: p=1, q=0  (s = 1-1-0 = 0, t=1)
    # Z0: p=0, q=1  (s = 1-0-1 = 0, t=0)

    types = {'M-': (1,1), 'M+': (0,0), 'Z1': (1,0), 'Z0': (0,1)}

    print(f"{'':>4} | ", end='')
    for wt in types:
        print(f"{wt:>4}", end=' ')
    print()
    print("-"*30)

    for ut, (pu, qu) in types.items():
        print(f"{ut:>4} | ", end='')
        for wt, (pw, qw) in types.items():
            # B = pu*qw - (1-qu)*(1-pw)
            B = pu * qw - (1 - qu) * (1 - pw)
            print(f"{B:>4}", end=' ')
        print()

    print("\nKey observations:")
    print("- Z0 rows are all zero (vertices beaten by both i,j are invisible in u-position)")
    print("- Z1 columns are all zero (vertices beating both i,j are invisible in w-position)")
    print("- M- row: (1, 0, 0, 1) -> contributes from M- and Z0 neighbors")
    print("- M+ row: (0, -1, 0, -1) -> contributes from M+ and Z0 neighbors (negatively)")
    print("- Z1 row: (1, -1, 0, 0) -> contributes from M- (positive) and M+ (negative)")


def test_effective_decomposition(T, i, j):
    """Test that the bracket structure correctly predicts Delta_H."""
    n = len(T)
    V0 = [v for v in range(n) if v != i and v != j]

    # Classify vertices
    types = {}
    for x in V0:
        types[x] = classify_vertex(T, x, i, j)

    # Compute actual Delta_H
    Tp = [row[:] for row in T]
    Tp[i][j] = 0
    Tp[j][i] = 1

    H_T = hamiltonian_path_count(T)
    H_Tp = hamiltonian_path_count(Tp)
    delta_H = H_T - H_Tp

    # Compute Delta_I
    cycles = find_odd_cycles(T)
    cycles_p = find_odd_cycles(Tp)
    I_T = independence_poly_at_fast(cycles, 2)
    I_Tp = independence_poly_at_fast(cycles_p, 2)
    delta_I = I_T - I_Tp

    print(f"\nn={n}, arc ({i},{j})")
    print(f"Types: {types}")
    type_counts = {}
    for t in types.values():
        type_counts[t] = type_counts.get(t, 0) + 1
    print(f"Counts: {type_counts}")
    print(f"Delta_H = {delta_H}, Delta_I = {delta_I}")

    # Verify bracket values for all (u,w) pairs
    for u in V0:
        for w in V0:
            if u == w:
                continue
            B = bracket(T, u, w, i, j)
            ut = types[u]
            wt = types[w]
            # Expected from table
            expected = {
                ('M-', 'M-'): 1, ('M-', 'M+'): 0, ('M-', 'Z1'): 0, ('M-', 'Z0'): 1,
                ('M+', 'M-'): 0, ('M+', 'M+'): -1, ('M+', 'Z1'): 0, ('M+', 'Z0'): -1,
                ('Z1', 'M-'): 1, ('Z1', 'M+'): -1, ('Z1', 'Z1'): 0, ('Z1', 'Z0'): 0,
                ('Z0', 'M-'): 0, ('Z0', 'M+'): 0, ('Z0', 'Z1'): 0, ('Z0', 'Z0'): 0,
            }
            exp = expected[(ut, wt)]
            if B != exp:
                print(f"  BRACKET MISMATCH: u={u}({ut}), w={w}({wt}), B={B}, expected={exp}")

    # Now test: does the decomposition
    # Delta_H = sum_S sum_{u in S, w in R} B(u,w) * h_end(S,u) * h_start(R,w) + boundaries
    # correctly give Delta_H?

    # Count "effective" vertex pair interactions
    eff_positive = 0  # (M-, M-), (M-, Z0), (Z1, M-)
    eff_negative = 0  # (M+, M+), (M+, Z0), (Z1, M+)
    for u in V0:
        for w in V0:
            if u == w:
                continue
            B = bracket(T, u, w, i, j)
            if B > 0:
                eff_positive += B
            elif B < 0:
                eff_negative += B

    print(f"Total positive brackets: {eff_positive}, negative: {eff_negative}")
    print(f"Effective pairs: positive from (M-,M-),(M-,Z0),(Z1,M-); negative from (M+,M+),(M+,Z0),(Z1,M+)")

    # The s-value decomposition
    s = {x: 1 - T[x][i] - T[j][x] for x in V0}
    s_sum = sum(s[x] for x in V0)
    print(f"sum(s_x) = {s_sum}, -2*sum(s_x) = {-2*s_sum}")

    # |M-| - |M+| relationship:
    nM = type_counts.get('M-', 0)
    nP = type_counts.get('M+', 0)
    nZ1 = type_counts.get('Z1', 0)
    nZ0 = type_counts.get('Z0', 0)
    print(f"|M-|={nM}, |M+|={nP}, |Z1|={nZ1}, |Z0|={nZ0}")
    print(f"|M-| - |M+| = {nM - nP} = -sum(s_x) = {-s_sum}")
    # Check: s_x for M- = -1, M+ = +1, Z = 0. So sum(s) = -|M-| + |M+|.
    # So |M-| - |M+| = -sum(s).


def main():
    verify_bracket_table()

    rng = random.Random(42)
    for n in [5, 6, 7]:
        for trial in range(2):
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
            test_effective_decomposition(T, ii, jj)


if __name__ == "__main__":
    main()
