"""
Algebraic analysis of WHY beta_2 = 0 for all tournaments.

Key structure:
- Omega_1 = A_1 = all directed edges (no constraint), dim = C(n,2)
- Omega_2: constrained by non-allowed 1-faces
  A path (a,b,c) has face (a,c) at position 1.
  (a,c) is NOT in A_1 iff c->a.
  So the constraint is: for each (a,c) with c->a,
  sum_{b: a->b, b->c} coeff(a,b,c) = 0
  (where the sum is over all intermediate vertices on 2-paths from a to c)

- The number of constraints = number of "backward pairs" (a,c) with c->a
  that have at least one intermediate b with a->b->c.
  For a tournament: this is the number of pairs (a,c) with c->a that are
  NOT in a "source-sink" relationship (no common successor-predecessor).
  Actually: the number of such pairs is exactly the number of 3-cycles
  (since (a,b,c) with a->b->c and c->a is a 3-cycle {a,b,c}).

Let me verify this and compute the constraint matrix structure.

Then: beta_2 = ker(d_2|Omega_2) - im(d_3|Omega_3).
We need to show this is always 0.

Since beta_2 = 0 while beta_4 can be nonzero, something special happens
at dimension 2 that doesn't happen at dimension 4.

Hypothesis: the constraint structure at Omega_2 is "rigid" —
the non-allowed faces create enough constraints to force exactness.
At Omega_4, the constraint structure has more slack.
"""

import numpy as np
from math import comb
from itertools import combinations

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def enumerate_allowed_paths(A, n, p):
    if p < 0: return []
    if p == 0: return [(v,) for v in range(n)]
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if A[i][j] == 1: adj[i].append(j)
    paths = []
    stack = []
    for start in range(n):
        stack.append(([start], 1 << start))
        while stack:
            path, visited = stack.pop()
            if len(path) == p + 1:
                paths.append(tuple(path))
                continue
            v = path[-1]
            for u in adj[v]:
                if not (visited & (1 << u)):
                    stack.append((path + [u], visited | (1 << u)))
    return paths

def boundary_coeffs(path):
    return [((-1)**i, path[:i] + path[i+1:]) for i in range(len(path))]

def main():
    print("=" * 70)
    print("ALGEBRAIC ANALYSIS OF BETA_2 = 0")
    print("=" * 70)

    # Part 1: Constraint structure for Omega_2
    print("\n--- Part 1: Omega_2 constraint analysis ---")

    for n in [4, 5, 6]:
        print(f"\n  n={n}:")
        N = 2**(n*(n-1)//2)

        for trial in range(min(N, 5)):
            A = bits_to_adj(trial, n)

            a1 = enumerate_allowed_paths(A, n, 1)
            a2 = enumerate_allowed_paths(A, n, 2)
            a1_set = set(a1)

            # Find non-allowed faces
            na_faces = {}  # face -> list of contributing paths
            for path in a2:
                a, b, c = path
                for i, face_entry in enumerate(boundary_coeffs(path)):
                    sign, face = face_entry
                    if len(face) == 2 and face not in a1_set:
                        if face not in na_faces:
                            na_faces[face] = []
                        na_faces[face].append((sign, path))

            # Count 3-cycles
            c3 = 0
            for i, j, k in combinations(range(n), 3):
                if A[i][j] and A[j][k] and A[k][i]: c3 += 1
                if A[i][k] and A[k][j] and A[j][i]: c3 += 1

            # Number of backward pairs with intermediate
            backward_pairs = 0
            for a_v in range(n):
                for c_v in range(n):
                    if a_v != c_v and A[c_v][a_v]:  # c->a
                        has_intermediate = False
                        for b_v in range(n):
                            if b_v != a_v and b_v != c_v and A[a_v][b_v] and A[b_v][c_v]:
                                has_intermediate = True
                                break
                        if has_intermediate:
                            backward_pairs += 1

            print(f"    trial {trial}: |A_2|={len(a2)}, non-allowed faces={len(na_faces)}, "
                  f"c3={c3}, backward_pairs_with_intermediate={backward_pairs}")

            # Each non-allowed face (a,c) with c->a comes from deleting
            # the middle vertex from a 2-path (a,b,c).
            # The sign is (-1)^1 = -1.
            # So the constraint at face (a,c) is:
            # -sum_{b: a->b, b->c} coeff(a,b,c) = 0
            # i.e., sum_{b} coeff(a,b,c) = 0

    # Part 2: The constraint matrix for Omega_2 at n=5
    print("\n--- Part 2: Detailed constraint matrix at n=5 ---")

    for trial in range(min(2**(5*4//2), 10)):
        A = bits_to_adj(trial, 5)
        n = 5

        a1 = enumerate_allowed_paths(A, n, 1)
        a2 = enumerate_allowed_paths(A, n, 2)
        a1_set = set(a1)

        na_faces = {}
        for j, path in enumerate(a2):
            for sign, face in boundary_coeffs(path):
                if len(face) == 2 and tuple(face) not in a1_set:
                    face_t = tuple(face)
                    if face_t not in na_faces:
                        na_faces[face_t] = []
                    na_faces[face_t].append((sign, j))

        if len(na_faces) == 0:
            # Transitive tournament — no constraints, Omega_2 = A_2
            pass
        else:
            # Build constraint matrix
            na_list = list(na_faces.keys())
            P = np.zeros((len(na_list), len(a2)))
            for i, face in enumerate(na_list):
                for sign, j in na_faces[face]:
                    P[i, j] += sign

            sv = np.linalg.svd(P, compute_uv=False)
            rank = sum(s > 1e-10 for s in sv)
            dim_omega2 = len(a2) - rank

            c3 = sum(1 for i,j,k in combinations(range(n), 3)
                     if (A[i][j] and A[j][k] and A[k][i]) or
                        (A[i][k] and A[k][j] and A[j][i]))

            if trial < 5 or dim_omega2 != len(a2) - len(na_list):
                print(f"  trial {trial}: |A_2|={len(a2)}, constraints={len(na_list)}, "
                      f"rank(P)={rank}, dim(Omega_2)={dim_omega2}, c3={c3}")

    # Part 3: Key insight — the constraint rows
    # For each backward pair (a,c) with c->a, the constraint row has
    # a +1 for each path (a,b,c) with a->b->c. The number of such b
    # equals the number of common out-neighbors of a that are in-neighbors of c.
    # In other words: the constraint is sum_{b in N_out(a) cap N_in(c)} coeff(a,b,c) = 0.

    # For Omega_2: how many paths per constraint?
    print("\n--- Part 3: Paths per constraint at n=6 ---")
    for trial in [0, 1, 100, 1000, 5000]:
        A = bits_to_adj(trial, 6)
        n = 6

        a1 = enumerate_allowed_paths(A, n, 1)
        a2 = enumerate_allowed_paths(A, n, 2)
        a1_set = set(a1)

        na_faces = {}
        for j, path in enumerate(a2):
            for sign, face in boundary_coeffs(path):
                if len(face) == 2 and tuple(face) not in a1_set:
                    face_t = tuple(face)
                    if face_t not in na_faces:
                        na_faces[face_t] = []
                    na_faces[face_t].append((sign, j))

        paths_per = [len(v) for v in na_faces.values()]
        c3 = sum(1 for i,j,k in combinations(range(n), 3)
                 if (A[i][j] and A[j][k] and A[k][i]) or
                    (A[i][k] and A[k][j] and A[j][i]))

        print(f"  trial {trial}: c3={c3}, constraints={len(na_faces)}, "
              f"paths_per_constraint={sorted(paths_per) if len(paths_per) <= 20 else f'min={min(paths_per)}, max={max(paths_per)}, avg={np.mean(paths_per):.1f}'}")

    # Part 4: Compare constraints at Omega_2 vs Omega_4
    print("\n--- Part 4: Constraint ratio at Omega_p for various p ---")
    for n in [6, 7, 8]:
        print(f"\n  n={n}:")
        import random
        rng = random.Random(42)
        for _ in range(3):
            A = np.zeros((n, n), dtype=int)
            for i in range(n):
                for j in range(i+1, n):
                    if rng.random() < 0.5:
                        A[i][j] = 1
                    else:
                        A[j][i] = 1

            for p in range(1, min(n, 6)):
                ap = enumerate_allowed_paths(A, n, p)
                apm1 = enumerate_allowed_paths(A, n, p-1)
                apm1_set = set(apm1)

                na_count = 0
                for path in ap:
                    for sign, face in boundary_coeffs(path):
                        if len(set(face)) == len(face) and tuple(face) not in apm1_set:
                            na_count += 1
                            break  # Just counting paths with at least 1 NA face

                # Quick: # distinct NA faces
                na_faces = set()
                for path in ap:
                    for sign, face in boundary_coeffs(path):
                        if len(set(face)) == len(face) and tuple(face) not in apm1_set:
                            na_faces.add(tuple(face))

                if len(ap) > 0:
                    ratio = len(na_faces) / len(ap) if len(ap) > 0 else 0
                    print(f"    p={p}: |A_p|={len(ap)}, NA_faces={len(na_faces)}, "
                          f"ratio={ratio:.3f}")

if __name__ == '__main__':
    main()
