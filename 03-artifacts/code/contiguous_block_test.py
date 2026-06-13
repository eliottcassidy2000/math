#!/usr/bin/env python3
"""
Test contiguous block decomposition of paths using an arc.

For a cycle C using arc i->j, define:
  f_C(T, i->j) = #{Ham paths where C appears as contiguous block AND i is before j}

Conjecture: f_C(T, i->j) = 2 * H(T[V\\V(C)])

Also test: f_0(T, i->j) = #{paths using i->j with NO cycle as contiguous block}
And whether f_0(T, i->j) = f_0(T', j->i) under arc flip.

Instance: kind-pasteur-2026-03-05-S5
"""

import sys
sys.path.insert(0, r"C:\Users\Eliott\Documents\GitHub\math\03-artifacts\code")
from tournament_lib import *
import random


def flip_arc(T, i, j):
    Tp = [row[:] for row in T]
    Tp[i][j] = 1 - T[i][j]
    Tp[j][i] = 1 - T[j][i]
    return Tp


def enum_ham_paths(T):
    """Enumerate all Hamiltonian paths of T."""
    n = len(T)
    paths = []
    for perm in permutations(range(n)):
        valid = True
        for k in range(n-1):
            if not T[perm[k]][perm[k+1]]:
                valid = False
                break
        if valid:
            paths.append(perm)
    return paths


def paths_using_arc_enum(T, i, j):
    """Enumerate Ham paths using arc i->j (i immediately before j)."""
    return [p for p in enum_ham_paths(T)
            if any(p[k]==i and p[k+1]==j for k in range(len(p)-1))]


def cycle_is_contiguous_in_path(path, cycle):
    """Check if cycle C appears as contiguous block in path P, following cycle direction."""
    n = len(path)
    L = len(cycle)
    if L > n:
        return False

    # All L rotations of the cycle as linearizations
    for rot in range(L):
        block = tuple(cycle[(rot + m) % L] for m in range(L))
        # Check if block appears as contiguous substring of path
        for start in range(n - L + 1):
            if path[start:start+L] == block:
                return True
    return False


def cycle_contiguous_using_arc(path, cycle, i, j):
    """Check if cycle is contiguous in path AND the path has i immediately before j."""
    n = len(path)
    L = len(cycle)
    # First check i immediately before j
    has_arc = any(path[k]==i and path[k+1]==j for k in range(n-1))
    if not has_arc:
        return False
    # Check cycle is contiguous
    return cycle_is_contiguous_in_path(path, cycle)


def H_complement(T, verts):
    """H of complement subtournament. H(empty) = 1."""
    n = len(T)
    keep = [v for v in range(n) if v not in verts]
    if not keep:
        return 1
    k = len(keep)
    sub = [[0]*k for _ in range(k)]
    for a in range(k):
        for b in range(k):
            sub[a][b] = T[keep[a]][keep[b]]
    return hamiltonian_path_count(sub)


def uses_arc(cyc, a, b):
    L = len(cyc)
    return any(cyc[k] == a and cyc[(k+1) % L] == b for k in range(L))


def main():
    print("=== Contiguous Block Decomposition Test ===\n")

    rng = random.Random(42)

    for n in [4, 5]:
        print(f"--- n={n} ---")
        fc_match = fc_total = 0
        f0_match = f0_total = 0
        total = 0

        tournaments = list(all_tournaments(n))
        for T in tournaments:
            for i in range(n):
                for j in range(n):
                    if i == j or not T[i][j]:
                        continue

                    paths_ij = paths_using_arc_enum(T, i, j)
                    f_ij = len(paths_ij)

                    # Find cycles using i->j
                    all_cyc = find_odd_cycles(T)
                    cyc_ij = [c for c in all_cyc if uses_arc(c, i, j)]

                    # For each cycle, count paths where cycle is contiguous
                    f_cycle_total = 0
                    for c in cyc_ij:
                        fc = sum(1 for p in paths_ij
                                 if cycle_is_contiguous_in_path(p, c))
                        h_comp = H_complement(T, set(c))
                        predicted_fc = 2 * h_comp

                        fc_total += 1
                        if fc == predicted_fc:
                            fc_match += 1

                        f_cycle_total += fc

                    # f_0 = paths using i->j with no cycle contiguous
                    f0 = sum(1 for p in paths_ij
                             if not any(cycle_is_contiguous_in_path(p, c) for c in cyc_ij))

                    # Check f_0 flip invariance
                    Tp = flip_arc(T, i, j)
                    paths_ji = paths_using_arc_enum(Tp, j, i)
                    cyc_ji = [c for c in find_odd_cycles(Tp) if uses_arc(c, j, i)]
                    f0_prime = sum(1 for p in paths_ji
                                   if not any(cycle_is_contiguous_in_path(p, c) for c in cyc_ji))

                    f0_total += 1
                    if f0 == f0_prime:
                        f0_match += 1

                    total += 1

        print(f"  Total arc tests: {total}")
        print(f"  f_C = 2*H(comp) per cycle: {fc_match}/{fc_total}")
        print(f"  f_0 flip-invariant: {f0_match}/{f0_total}")

    # n=6 random
    print(f"\n--- n=6 (random 50) ---")
    fc_match = fc_total = 0
    f0_match = f0_total = 0
    total = 0

    for trial in range(50):
        T = random_tournament(6, rng)
        i, j = rng.sample(range(6), 2)
        if not T[i][j]:
            i, j = j, i

        paths_ij = paths_using_arc_enum(T, i, j)
        all_cyc = find_odd_cycles(T)
        cyc_ij = [c for c in all_cyc if uses_arc(c, i, j)]

        for c in cyc_ij:
            fc = sum(1 for p in paths_ij if cycle_is_contiguous_in_path(p, c))
            h_comp = H_complement(T, set(c))
            predicted_fc = 2 * h_comp
            fc_total += 1
            if fc == predicted_fc:
                fc_match += 1

        f0 = sum(1 for p in paths_ij
                 if not any(cycle_is_contiguous_in_path(p, c) for c in cyc_ij))

        Tp = flip_arc(T, i, j)
        paths_ji = paths_using_arc_enum(Tp, j, i)
        cyc_ji = [c for c in find_odd_cycles(Tp) if uses_arc(c, j, i)]
        f0_prime = sum(1 for p in paths_ji
                       if not any(cycle_is_contiguous_in_path(p, c) for c in cyc_ji))

        f0_total += 1
        if f0 == f0_prime:
            f0_match += 1

        total += 1

    print(f"  Total arc tests: {total}")
    print(f"  f_C = 2*H(comp) per cycle: {fc_match}/{fc_total}")
    print(f"  f_0 flip-invariant: {f0_match}/{f0_total}")


if __name__ == "__main__":
    main()
