#!/usr/bin/env python3
"""
Test inductive structure of the polynomial identity U_T' - U_T = delta_I.

Key idea: Can we express the n-vertex identity in terms of (n-1)-vertex identities?

For each vertex x in V\{i,j}, consider:
- Paths where x is at position 0 (start): prefix is just x, suffix through B_x
- Paths where x is at position n-1 (end): prefix through B_x, suffix is just x
- Paths where x is in the middle: x is between two B_x vertices

The unmatched count U_T decomposes by the position of x relative to (i,j).

Instance: opus-2026-03-05-S3
"""

import sys
sys.path.insert(0, '.')
from tournament_lib import hamiltonian_path_count, find_odd_cycles, independence_poly_at_fast
from itertools import permutations, combinations
import random


def nadj_count(T, seq):
    """Count Ham paths of T containing the consecutive subsequence `seq`.
    Uses brute force enumeration."""
    n = len(T)
    count = 0
    for perm in permutations(range(n)):
        # Check if perm is a valid Ham path
        valid = all(T[perm[k]][perm[k+1]] for k in range(n-1))
        if not valid:
            continue
        # Check if seq appears consecutively
        for start in range(n - len(seq) + 1):
            if all(perm[start + k] == seq[k] for k in range(len(seq))):
                count += 1
                break
    return count


def test_nadj_decomposition(T, i, j):
    """Test if U_T = sum_{x: s=-1} [nadj(x,i,j) + nadj(i,j,x)] - double_count
    and U_T' = sum_{x: s=+1} [nadj'(x,j,i) + nadj'(j,i,x)] - double_count."""
    n = len(T)
    others = [v for v in range(n) if v != i and v != j]

    Tp = [row[:] for row in T]
    Tp[i][j] = 0
    Tp[j][i] = 1

    s = {x: 1 - T[x][i] - T[j][x] for x in others}

    # Compute nadj terms for T
    for x in others:
        if s[x] == -1:
            n_xij = nadj_count(T, (x, i, j))
            n_ijx = nadj_count(T, (i, j, x))
            n_xijx_impossible = 0  # x can't appear twice

            # What about paths blocked at BOTH pred and succ?
            # (..., x, i, j, y, ...) where both x blocks AND y blocks
            n_both = 0
            for y in others:
                if y != x and s[y] == -1:
                    n_both += nadj_count(T, (x, i, j, y))

            print(f"  T-side: x={x}, s={s[x]}, nadj(x,i,j)={n_xij}, "
                  f"nadj(i,j,x)={n_ijx}, both={n_both}")

    for x in others:
        if s[x] == 1:
            n_xji = nadj_count(Tp, (x, j, i))
            n_jix = nadj_count(Tp, (j, i, x))
            n_both = 0
            for y in others:
                if y != x and s[y] == 1:
                    n_both += nadj_count(Tp, (x, j, i, y))
            print(f"  T'-side: x={x}, s={s[x]}, nadj'(x,j,i)={n_xji}, "
                  f"nadj'(j,i,x)={n_jix}, both={n_both}")

    # Key relationship: nadj_T(x, i, j) = #Ham paths of T with x immediately before i
    # and i immediately before j. This is a "3-consecutive" count.
    # nadj_T(x, i, j) = sum over subsets S of B_x:
    #   [#paths through {x}∪S ending at x in T] * [#paths through {j}∪(B_x\S) starting at j in T]
    # This is the factorization through the cut point.

    # The crucial observation: nadj_T(x, i, j) involves paths in B_x that END at some
    # vertex adjacent to x, then x->i->j, then paths in B_x\S starting from j.
    # This is exactly a convolution over B_x partitions.

    # And H(B_x) = sum over all partitions S of B_x:
    #   [#paths through S ending at some v] * [#paths through B_x\S starting from some w]
    # with v->w being an arc.

    # So nadj_T(x,i,j) and nadj_T(i,j,x) are "projections" of the B_x path structure
    # through the arcs x->i and j->x.


def test_key_relationship(T, i, j):
    """Test: nadj(x,i,j) + nadj(i,j,x) = H(B_x) when s_x=-1? NO, that's wrong.

    But maybe: sum_{x: s=-1} [nadj(x,i,j) + nadj(i,j,x)]
             - sum_{x,y: s=-1} nadj(x,i,j,y) = U_T

    And the T' version:
    sum_{x: s=+1} [nadj'(x,j,i) + nadj'(j,i,x)]
             - sum_{x,y: s=+1} nadj'(x,j,i,y) = U_T'

    So U_T'-U_T = [sum_{s=+1} stuff] - [sum_{s=-1} stuff].

    Can we relate nadj to H(B_x)?
    """
    n = len(T)
    others = [v for v in range(n) if v != i and v != j]
    Tp = [row[:] for row in T]
    Tp[i][j] = 0
    Tp[j][i] = 1
    s = {x: 1 - T[x][i] - T[j][x] for x in others}

    print(f"\nn={n}, arc ({i},{j}), s-values: {s}")

    # For each x, compute nadj_T(x,i,j) and nadj_T(i,j,x) and H(B_x)
    for x in others:
        B_x = [v for v in others if v != x]
        # H(B_x)
        m = len(B_x)
        sub = [[0]*m for _ in range(m)]
        idx_map = {v: k for k, v in enumerate(B_x)}
        for a in B_x:
            for b in B_x:
                if a != b:
                    sub[idx_map[a]][idx_map[b]] = T[a][b]
        H_Bx = hamiltonian_path_count(sub)

        n_xij = nadj_count(T, (x, i, j))
        n_ijx = nadj_count(T, (i, j, x))
        n_xji_p = nadj_count(Tp, (x, j, i))
        n_jix_p = nadj_count(Tp, (j, i, x))

        # nadj_T(x,i,j): paths (..., x, i, j, ...) in T
        # x->i requires T[x][i]=1, i->j requires T[i][j]=1 (=1 by construction)
        # So nadj_T(x,i,j) = 0 if T[x][i]=0.
        # Similarly nadj_T(i,j,x) = 0 if T[j][x]=0.

        # For T'-paths: nadj_{T'}(x,j,i) requires T'[x][j]=T[x][j] and T'[j][i]=1
        # So nadj_{T'}(x,j,i) = 0 if T[x][j]=0.
        # nadj_{T'}(j,i,x) requires T'[i][x]=T[i][x].

        print(f"  x={x}: s={s[x]}, H(B_x)={H_Bx}, "
              f"T: nadj(x,i,j)={n_xij}, nadj(i,j,x)={n_ijx}, sum={n_xij+n_ijx} | "
              f"T': nadj'(x,j,i)={n_xji_p}, nadj'(j,i,x)={n_jix_p}, sum={n_xji_p+n_jix_p}")

        # Key test: is there a relationship between nadj sums and H(B_x)?
        # For x with s_x = -1: T[x][i]=1, T[j][x]=1
        #   nadj(x,i,j) counts paths with x->i->j (x beats i, arc i->j exists)
        #   nadj(i,j,x) counts paths with i->j->x (arc i->j, j beats x) -- but T[j][x]=1
        #   Both are nonzero in general.
        # For x with s_x = +1: T[x][i]=0, T[j][x]=0
        #   nadj(x,i,j) = 0 (need T[x][i]=1)
        #   nadj(i,j,x) = 0 (need T[j][x]=1)
        # For x with s_x = 0: exactly one of T[x][i], T[j][x] is 1.


def main():
    rng = random.Random(42)
    for n in [5, 6]:
        print(f"\n{'='*70}")
        print(f"Testing n={n}")

        for trial in range(3):
            T = [[0]*n for _ in range(n)]
            for a in range(n):
                for b in range(a+1, n):
                    if rng.random() < 0.5:
                        T[a][b] = 1
                    else:
                        T[b][a] = 1
            i, j = 0, 1
            if T[i][j] == 0:
                i, j = j, i
            test_key_relationship(T, i, j)


if __name__ == "__main__":
    main()
