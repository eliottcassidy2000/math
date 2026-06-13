#!/usr/bin/env python3
"""
Analyze the full tiling class size = 2^(n-2)+1 conjecture.

Key question: WHY is the full class (all arcs flipped from transitive)
exactly 2^(n-2)+1 for n >= 4?

The full tiling corresponds to the "anti-transitive" tournament: the unique
tournament T on {1,...,n} where vertex k beats vertex l iff k < l (the
REVERSE of transitive where k beats l iff k > l).

Class size = H(T) / |Aut(T)|. If |Aut| = 1, then class_size = H(T).
So we need: H(anti-transitive_n) = 2^(n-2) + 1 for n >= 4.

The anti-transitive tournament is the transitive tournament with ALL arcs
reversed. But actually, the "full tiling" tournament depends on the choice
of fixed Hamiltonian path and grid encoding.

Let's verify what tournament the full tiling gives, and compute H(T) directly.

Instance: kind-pasteur-2026-03-05-S10
"""

from itertools import permutations

def full_tiling_tournament(n):
    """Build the tournament from the all-ones tiling."""
    verts = list(range(n, 0, -1))  # [n, n-1, ..., 1]

    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            tiles.append((x, y))

    # All bits = 1
    A = [[0]*n for _ in range(n)]
    # Fixed Ham path: indices 0->1->...->n-1 (labels n->n-1->...->1)
    for k in range(n-1):
        A[k][k+1] = 1

    # Each tile (x,y) with bit=1: arc goes yi -> xi (reversed from default)
    for (xL, yL) in tiles:
        xi = verts.index(xL)
        yi = verts.index(yL)
        A[yi][xi] = 1

    return A


def transitive_tournament(n):
    """Build the transitive tournament (all tiles off)."""
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            A[i][j] = 1  # vertex i (label n-i) beats vertex j (label n-j)
    return A


def count_ham_paths(A, n):
    """Count Hamiltonian paths."""
    count = 0
    for p in permutations(range(n)):
        if all(A[p[k]][p[k+1]] for k in range(n-1)):
            count += 1
    return count


def count_automorphisms(A, n):
    """Count automorphisms."""
    count = 0
    for p in permutations(range(n)):
        if all(A[p[i]][p[j]] == A[i][j] for i in range(n) for j in range(n)):
            count += 1
    return count


def score_sequence(A, n):
    """Out-degree sequence (sorted descending)."""
    return tuple(sorted([sum(A[i]) for i in range(n)], reverse=True))


def analyze_full_tiling(n):
    """Analyze the full tiling tournament at size n."""
    A_full = full_tiling_tournament(n)
    A_trans = transitive_tournament(n)

    h_full = count_ham_paths(A_full, n)
    h_trans = count_ham_paths(A_trans, n)
    aut_full = count_automorphisms(A_full, n)
    scores_full = score_sequence(A_full, n)

    print(f"\nn={n}:")
    print(f"  Full tiling tournament scores: {scores_full}")
    print(f"  Transitive tournament scores: {score_sequence(A_trans, n)}")
    print(f"  H(full) = {h_full}")
    print(f"  H(trans) = {h_trans}")
    print(f"  |Aut(full)| = {aut_full}")
    print(f"  class_size = H/|Aut| = {h_full // aut_full}")
    print(f"  2^(n-2)+1 = {2**(n-2)+1}")
    print(f"  Match: {h_full // aut_full == 2**(n-2)+1}")

    # Check: is the full tiling the same as the reversed transitive?
    A_rev = [[A_trans[j][i] for j in range(n)] for i in range(n)]
    # Check if full == reversed transitive
    is_rev = all(A_full[i][j] == A_rev[i][j] for i in range(n) for j in range(n))
    print(f"  Full = reversed transitive? {is_rev}")

    # Display adjacency matrix
    print(f"  Full adj matrix:")
    for row in A_full:
        print(f"    {row}")

    return h_full


def count_ham_paths_starting_at(A, n, start):
    """Count Ham paths starting at vertex 'start'."""
    count = 0
    for p in permutations(range(n)):
        if p[0] != start:
            continue
        if all(A[p[k]][p[k+1]] for k in range(n-1)):
            count += 1
    return count


def analyze_ham_path_structure(n):
    """Analyze the structure of Ham paths in the full tiling tournament."""
    A = full_tiling_tournament(n)
    h = count_ham_paths(A, n)

    print(f"\nn={n}: H(full) = {h}")
    print(f"  Ham paths by starting vertex:")
    for v in range(n):
        h_v = count_ham_paths_starting_at(A, n, v)
        print(f"    start={v} (label {n-v}): {h_v} paths")

    print(f"  Ham paths by ending vertex:")
    for v in range(n):
        count = 0
        for p in permutations(range(n)):
            if p[-1] != v:
                continue
            if all(A[p[k]][p[k+1]] for k in range(n-1)):
                count += 1
        print(f"    end={v} (label {n-v}): {count} paths")


if __name__ == "__main__":
    for n in range(3, 9):
        analyze_full_tiling(n)

    for n in range(3, 8):
        analyze_ham_path_structure(n)
