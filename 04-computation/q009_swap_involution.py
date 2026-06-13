"""
Q-009: Swap involution approach to proving adj(i,j) - adj'(j,i) = DeltaI.

For a T-path pi = (..., x, i, j, y, ...), define swap(pi) = (..., x, j, i, y, ...).
swap(pi) is a valid T'-path iff:
  - T[x][j] = 1 (predecessor of i must also beat j)
  - T[i][y] = 1 (i must beat successor of j)
  (boundary: if i is first, no x constraint; if j is last, no y constraint)

Matched T-paths pair perfectly with matched T'-paths.
So: adj(i,j) - adj'(j,i) = #unmatched_T - #unmatched_T'

An unmatched T-path (..., x, i, j, y, ...) has:
  T[x][j]=0 (j->x, so s_x=-1) OR T[i][y]=0 (y->i, so s_y=-1)

An unmatched T'-path (..., x, j, i, y, ...) has:
  T[x][i]=0 (i->x, so s_x=+1) OR T[j][y]=0 (y->j... wait)

Let's verify and analyze computationally.

Author: opus-2026-03-05-S2
"""

import random
from itertools import permutations


def random_tournament(n):
    T = [[0]*n for _ in range(n)]
    for a in range(n):
        for b in range(a+1, n):
            if random.random() < 0.5:
                T[a][b] = 1
            else:
                T[b][a] = 1
    return T


def flip_arc(T, i, j):
    T2 = [row[:] for row in T]
    T2[i][j] = 0
    T2[j][i] = 1
    return T2


def ham_paths(T, n):
    """Return all Ham paths of T."""
    paths = []
    for perm in permutations(range(n)):
        valid = True
        for k in range(n - 1):
            if T[perm[k]][perm[k+1]] != 1:
                valid = False
                break
        if valid:
            paths.append(perm)
    return paths


if __name__ == "__main__":
    random.seed(42)

    for n in [5, 6, 7]:
        print(f"\n{'='*60}")
        print(f"n = {n}")
        n_trials = {5: 20, 6: 10, 7: 5}[n]

        for trial in range(n_trials):
            T = random_tournament(n)
            i, j = random.sample(range(n), 2)
            if T[i][j] == 0:
                i, j = j, i

            T2 = flip_arc(T, i, j)
            others = [x for x in range(n) if x != i and x != j]

            # Find all Ham paths
            paths_T = ham_paths(T, n)
            paths_T2 = ham_paths(T2, n)

            # T-paths using i->j
            using_ij = []
            for p in paths_T:
                for k in range(n-1):
                    if p[k] == i and p[k+1] == j:
                        using_ij.append((p, k))

            # T'-paths using j->i
            using_ji = []
            for p in paths_T2:
                for k in range(n-1):
                    if p[k] == j and p[k+1] == i:
                        using_ji.append((p, k))

            adj_ij = len(using_ij)
            adj_ji = len(using_ji)

            # Classify: matched vs unmatched
            matched_T = 0
            unmatched_T = []
            for (p, k) in using_ij:
                # Predecessor and successor
                has_pred = k > 0
                has_succ = k + 1 < n - 1
                pred_ok = (not has_pred) or T[p[k-1]][j] == 1
                succ_ok = (not has_succ) or T[i][p[k+2]] == 1
                if pred_ok and succ_ok:
                    matched_T += 1
                else:
                    # Classify blocking
                    block_pred = has_pred and T[p[k-1]][j] == 0
                    block_succ = has_succ and T[i][p[k+2]] == 0
                    unmatched_T.append((p, k, block_pred, block_succ))

            matched_T2 = 0
            unmatched_T2 = []
            for (p, k) in using_ji:
                has_pred = k > 0
                has_succ = k + 1 < n - 1
                pred_ok = (not has_pred) or T[p[k-1]][i] == 1  # Note: T not T', since these arcs unchanged
                succ_ok = (not has_succ) or T[j][p[k+2]] == 1
                if pred_ok and succ_ok:
                    matched_T2 += 1
                else:
                    block_pred = has_pred and T[p[k-1]][i] == 0
                    block_succ = has_succ and T[j][p[k+2]] == 0
                    unmatched_T2.append((p, k, block_pred, block_succ))

            delta = adj_ij - adj_ji
            unmatched_delta = len(unmatched_T) - len(unmatched_T2)

            # Verify matched counts equal
            if matched_T != matched_T2:
                print(f"  *** MATCHED MISMATCH: T={matched_T}, T'={matched_T2}")

            if trial < 3:
                print(f"\nTrial {trial}: flip {i}->{j}, delta={delta}")
                print(f"  adj(i,j)={adj_ij}: {matched_T} matched + {len(unmatched_T)} unmatched")
                print(f"  adj'(j,i)={adj_ji}: {matched_T2} matched + {len(unmatched_T2)} unmatched")
                print(f"  delta = {len(unmatched_T)} - {len(unmatched_T2)} = {unmatched_delta}")
                print(f"  Match: {unmatched_delta == delta}")

                # Analyze unmatched T-paths by blocking type
                pred_only = sum(1 for (_, _, bp, bs) in unmatched_T if bp and not bs)
                succ_only = sum(1 for (_, _, bp, bs) in unmatched_T if bs and not bp)
                both = sum(1 for (_, _, bp, bs) in unmatched_T if bp and bs)
                print(f"  Unmatched T: pred_block={pred_only}, succ_block={succ_only}, both={both}")

                pred_only2 = sum(1 for (_, _, bp, bs) in unmatched_T2 if bp and not bs)
                succ_only2 = sum(1 for (_, _, bp, bs) in unmatched_T2 if bs and not bp)
                both2 = sum(1 for (_, _, bp, bs) in unmatched_T2 if bp and bs)
                print(f"  Unmatched T': pred_block={pred_only2}, succ_block={succ_only2}, both={both2}")

                # s_x values of blocking vertices
                sx = {x: 1 - T[x][i] - T[j][x] for x in others}
                print(f"  s_x: {sx}")

                # For unmatched T-paths, what are the blocking vertices?
                block_verts_T = {}
                for (p, k, bp, bs) in unmatched_T:
                    if bp:
                        x = p[k-1]
                        block_verts_T[x] = block_verts_T.get(x, 0) + 1
                    if bs:
                        y = p[k+2]
                        block_verts_T[y] = block_verts_T.get(y, 0) + 1

                block_verts_T2 = {}
                for (p, k, bp, bs) in unmatched_T2:
                    if bp:
                        x = p[k-1]
                        block_verts_T2[x] = block_verts_T2.get(x, 0) + 1
                    if bs:
                        y = p[k+2]
                        block_verts_T2[y] = block_verts_T2.get(y, 0) + 1

                print(f"  Blocking vertices T: {block_verts_T}")
                print(f"  Blocking vertices T': {block_verts_T2}")

            # Verify match for all trials
            if unmatched_delta != delta:
                print(f"  *** DELTA MISMATCH at trial {trial}!")
