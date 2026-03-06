#!/usr/bin/env python3
"""
Analyze the inclusion-exclusion structure of U_T - U_T' to understand
what a general proof requires.

The key identity (THM-014 + THM-015):
  delta_H = H(T') - H(T) = U_T' - U_T

And (THM-013):
  delta_I = I(Omega(T'), 2) - I(Omega(T), 2) = sum 2^k Delta(alpha_k)

By inclusion-exclusion over blocking vertices:
  U_T = sum_{x: s=-1} P_T(x) - sum_{x<y: s_x=s_y=-1} P_T(x,y) + ...
  U_T' = sum_{x: s=+1} P_T'(x) - sum_{x<y: s_x=s_y=+1} P_T'(x,y) + ...

where P_T(x) = nadj(x,i,j) + nadj(i,j,x) = #(T-paths using i->j with x as neighbor)
      P_T(x,y) = nadj(x,i,j,y) = #(T-paths with x before i and y after j)

Key insight to test: can P_T(x) be expressed in terms of H(B_x)?

nadj(x,i,j) = sum_{S subset B_x} h_end(T[S+{x}], x) * h_start(T[B_x\S + {j}], j)
            but j is not in B_x... Let me be more careful.

Actually: nadj_T(x,i,j) counts Ham paths of T of the form (prefix, x, i, j, suffix)
where prefix = Ham path through some S subset B_x starting somewhere ending at v -> x,
and suffix = Ham path through B_x\S starting at j -> w.

So nadj_T(x,i,j) = sum_{S subset B_x} [#paths through S ending adj to x] * [#paths through B_x\S starting from j adj to w]

Hmm, this is a convolution. Let me compute it directly and compare.

Instance: opus-2026-03-05-S3
"""

import sys
sys.path.insert(0, '.')
from tournament_lib import hamiltonian_path_count, find_odd_cycles, all_tournaments
from itertools import permutations, combinations
import random


def nadj(T, seq):
    """Count Ham paths containing consecutive subsequence `seq`."""
    n = len(T)
    count = 0
    for perm in permutations(range(n)):
        if not all(T[perm[k]][perm[k+1]] for k in range(n-1)):
            continue
        for start in range(n - len(seq) + 1):
            if all(perm[start+k] == seq[k] for k in range(len(seq))):
                count += 1
                break
    return count


def h_sub(T, vertices):
    """Ham path count of sub-tournament on given vertices."""
    if len(vertices) <= 1:
        return 1
    vlist = sorted(vertices)
    m = len(vlist)
    idx = {v: k for k, v in enumerate(vlist)}
    sub = [[0]*m for _ in range(m)]
    for a in vlist:
        for b in vlist:
            if a != b:
                sub[idx[a]][idx[b]] = T[a][b]
    return hamiltonian_path_count(sub)


def h_end(T, vertices, target):
    """#Ham paths through `vertices` ending at `target`."""
    if len(vertices) == 1:
        return 1 if target in vertices else 0
    vlist = sorted(vertices)
    m = len(vlist)
    idx = {v: k for k, v in enumerate(vlist)}
    sub = [[0]*m for _ in range(m)]
    for a in vlist:
        for b in vlist:
            if a != b:
                sub[idx[a]][idx[b]] = T[a][b]
    # DP
    dp = {}
    for v in range(m):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << m):
        for v in range(m):
            if not ((mask >> v) & 1):
                continue
            if (mask, v) not in dp:
                continue
            for u in range(m):
                if (mask >> u) & 1:
                    continue
                if sub[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + dp[(mask, v)]
    full = (1 << m) - 1
    return dp.get((full, idx[target]), 0)


def h_start(T, vertices, source):
    """#Ham paths through `vertices` starting at `source`."""
    if len(vertices) == 1:
        return 1 if source in vertices else 0
    vlist = sorted(vertices)
    m = len(vlist)
    idx = {v: k for k, v in enumerate(vlist)}
    sub = [[0]*m for _ in range(m)]
    for a in vlist:
        for b in vlist:
            if a != b:
                sub[idx[a]][idx[b]] = T[a][b]
    # DP starting from source
    dp = {(1 << idx[source], idx[source]): 1}
    for mask in range(1, 1 << m):
        for v in range(m):
            if not ((mask >> v) & 1):
                continue
            if (mask, v) not in dp:
                continue
            for u in range(m):
                if (mask >> u) & 1:
                    continue
                if sub[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + dp[(mask, v)]
    full = (1 << m) - 1
    return sum(dp.get((full, v), 0) for v in range(m))


def analyze_nadj_factorization(T, i, j):
    """Test: does nadj(x,i,j) factorize through B_x sub-tournament?"""
    n = len(T)
    others = [v for v in range(n) if v != i and v != j]
    s = {x: 1 - T[x][i] - T[j][x] for x in others}

    print(f"\nn={n}, arc ({i},{j}), s-values: {s}")

    for x in others:
        B_x = [v for v in others if v != x]

        # nadj(x,i,j): paths (..., x, i, j, ...)
        # Requires T[x][i]=1. If T[x][i]=0, this is 0.
        n_xij = nadj(T, (x, i, j))

        # Factorization: nadj(x,i,j) = sum_{S subset B_x}
        #   h_end(T[S union {x}], x) * h_start(T[B_x\S], first_after_j)
        # where first_after_j must satisfy T[j][first]=1.
        # More precisely: the suffix is a Ham path through B_x\S starting at some w
        # with T[j][w]=1.

        # Compute factorized version
        n_xij_fact = 0
        for mask in range(1 << len(B_x)):
            S = [B_x[k] for k in range(len(B_x)) if (mask >> k) & 1]
            R = [B_x[k] for k in range(len(B_x)) if not ((mask >> k) & 1)]
            S_plus_x = S + [x]

            he = h_end(T, S_plus_x, x) if T[x][i] else 0

            # h_start from j: paths through R starting at w with T[j][w]=1
            hs = 0
            for w in R:
                if T[j][w]:
                    hs += h_start(T, R, w)

            n_xij_fact += he * hs

        # Similarly for nadj(i,j,x): paths (..., i, j, x, ...)
        # Requires T[j][x]=1.
        n_ijx = nadj(T, (i, j, x))

        n_ijx_fact = 0
        for mask in range(1 << len(B_x)):
            S = [B_x[k] for k in range(len(B_x)) if (mask >> k) & 1]
            R = [B_x[k] for k in range(len(B_x)) if not ((mask >> k) & 1)]

            # Prefix ends at some v with T[v][i]=1, path through S
            hp = 0
            for v in S:
                if T[v][i]:
                    hp += h_end(T, S, v)

            # Suffix starts at x, path through R+{x}
            R_plus_x = R + [x]
            hs = h_start(T, R_plus_x, x) if T[j][x] else 0

            n_ijx_fact += hp * hs

        print(f"  x={x}: s={s[x]}, nadj(x,i,j)={n_xij} (fact={n_xij_fact}), "
              f"nadj(i,j,x)={n_ijx} (fact={n_ijx_fact}), "
              f"H(B_x)={h_sub(T, B_x)}")

        if n_xij != n_xij_fact or n_ijx != n_ijx_fact:
            print(f"    FACTORIZATION MISMATCH!")

        # Key test: for s_x=-1, is nadj(x,i,j) + nadj(i,j,x) related to H(B_x)?
        if s[x] == -1:
            total = n_xij + n_ijx
            H_Bx = h_sub(T, B_x)
            print(f"    s=-1: nadj_sum={total}, H(B_x)={H_Bx}, ratio={total/H_Bx if H_Bx else 'inf'}")

        if s[x] == 1:
            # For T'-paths
            Tp = [row[:] for row in T]
            Tp[i][j] = 0
            Tp[j][i] = 1
            n_xji_p = nadj(Tp, (x, j, i))
            n_jix_p = nadj(Tp, (j, i, x))
            total_p = n_xji_p + n_jix_p
            H_Bx = h_sub(T, B_x)
            print(f"    s=+1: nadj'_sum={total_p}, H(B_x)={H_Bx}, ratio={total_p/H_Bx if H_Bx else 'inf'}")


def main():
    rng = random.Random(42)
    for n in [5, 6, 7]:
        print(f"\n{'='*70}")
        print(f"n = {n}")
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
            analyze_nadj_factorization(T, ii, jj)


if __name__ == "__main__":
    main()
