"""
Omega_2 formula: Can dim(Omega_2) be expressed in terms of c3?

For p=2, a 2-path is (a,b,c) with a->b and b->c (distinct vertices).
|A_2| = sum_{|S|=3} H(T[S]) = C(n,3) + 2*c3.

The Omega_2 constraint: the boundary of a path in A_2 must have all its
faces in A_1 (= set of tournament edges). But face 1 = (a,c) may have
the REVERSE edge c->a, making it non-allowed.

For each "non-allowed face" (u,v) with v->u:
  The constraint: sum of (-1)^1 * [all 2-paths (u,w,v) with u->w, w->v] = 0
  i.e., sum of paths through (u,?,v) where u->?, ?->v must vanish in Omega_2.

The number of non-allowed constraints equals the number of ordered pairs
(u,v) where v->u AND there exists w with u->w->v.

Let's compute dim(Omega_2) and check if it equals |A_2| - rank(constraints).
"""
import numpy as np
from math import comb
from collections import defaultdict
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

def random_tournament(n, rng):
    A = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
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

def compute_omega_dim(A, n, p, allowed_p, allowed_pm1):
    dim_Ap = len(allowed_p)
    if dim_Ap == 0: return 0
    if p == 0: return dim_Ap
    allowed_pm1_set = set(allowed_pm1)
    non_allowed = {}
    na_count = 0
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if len(set(face)) == len(face) and face not in allowed_pm1_set:
                if face not in non_allowed:
                    non_allowed[face] = na_count
                    na_count += 1
    if na_count == 0: return dim_Ap
    P = np.zeros((na_count, dim_Ap))
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if face in non_allowed:
                P[non_allowed[face], j] += sign
    sv = np.linalg.svd(P, compute_uv=False)
    rank = sum(s > 1e-10 for s in sv)
    return dim_Ap - rank

def main():
    print("=" * 70)
    print("OMEGA_2 FORMULA: dim(Omega_2) in terms of tournament invariants")
    print("=" * 70)

    # Part 1: Exhaustive at n=3,4,5,6; sampled at n=7,8
    for n in range(3, 9):
        rng = np.random.RandomState(42 + n)
        if n <= 5:
            N = 2**(n*(n-1)//2)
        elif n == 6:
            N = 2**(n*(n-1)//2)  # exhaustive at n=6 too
        else:
            N = 500

        data = []  # (c3, omega_2, A2, n_na_faces, rank_constraints)

        for trial in range(N):
            if n <= 6:
                A = bits_to_adj(trial, n)
            else:
                A = random_tournament(n, rng)

            # Count c3
            c3 = 0
            for a, b, c in combinations(range(n), 3):
                if A[a][b] and A[b][c] and A[c][a]: c3 += 1
                if A[a][c] and A[c][b] and A[b][a]: c3 += 1

            # Compute |A_2| and Omega_2
            a2 = enumerate_allowed_paths(A, n, 2)
            a1 = enumerate_allowed_paths(A, n, 1)
            od = compute_omega_dim(A, n, 2, a2, a1)

            # Count non-allowed faces
            a1_set = set(a1)
            na_faces = set()
            for path in a2:
                for sign, face in boundary_coeffs(path):
                    if len(set(face)) == len(face) and face not in a1_set:
                        na_faces.add(face)
            n_na = len(na_faces)

            # Count rank of constraint matrix
            rank = len(a2) - od

            data.append((c3, od, len(a2), n_na, rank))

        print(f"\n  n={n} (d={n-1}), {N} samples:")
        print(f"    |A_2| = C(n,3) + 2*c3 = {comb(n,3)} + 2*c3")

        # Check if Omega_2 is a function of c3 alone
        by_c3 = defaultdict(list)
        for c3, od, a2, n_na, rank in data:
            by_c3[c3].append((od, a2, n_na, rank))

        deterministic = True
        for c3_val in sorted(by_c3.keys()):
            ods = [d[0] for d in by_c3[c3_val]]
            ranks = [d[3] for d in by_c3[c3_val]]
            n_nas = [d[2] for d in by_c3[c3_val]]

            od_set = sorted(set(ods))
            rank_set = sorted(set(ranks))
            c_val = comb(n, 3) + 2 * c3_val

            if len(od_set) > 1:
                deterministic = False

            if len(by_c3[c3_val]) >= 3:  # Only report if enough samples
                print(f"    c3={c3_val}: Omega_2 in {od_set}, |A_2|={c_val}, "
                      f"rank in {rank_set}, count={len(by_c3[c3_val])}")

        print(f"    Omega_2 is {'a function of c3 alone' if deterministic else 'NOT a function of c3 alone'}")

        # Is there a formula? Try: Omega_2 = a*c3 + b*n + c
        if n >= 4:
            c3_arr = np.array([d[0] for d in data], dtype=float)
            od_arr = np.array([d[1] for d in data], dtype=float)

            # Linear fit: Omega_2 = alpha * c3 + beta
            if len(set(c3_arr)) > 1:
                coeffs = np.polyfit(c3_arr, od_arr, 1)
                pred = np.polyval(coeffs, c3_arr)
                residual = np.max(np.abs(od_arr - pred))
                print(f"    Linear fit: Omega_2 ~ {coeffs[0]:.4f}*c3 + {coeffs[1]:.4f}, max residual={residual:.4f}")

            # Quadratic fit
            if len(set(c3_arr)) > 2:
                coeffs2 = np.polyfit(c3_arr, od_arr, 2)
                pred2 = np.polyval(coeffs2, c3_arr)
                residual2 = np.max(np.abs(od_arr - pred2))
                print(f"    Quadratic fit: max residual={residual2:.4f}")

    # Part 2: What ADDITIONAL invariant is needed beyond c3?
    # At n >= 5, Omega_2 is not a function of c3 alone.
    # Try: score sequence, or specific graph properties
    print("\n" + "=" * 70)
    print("Part 2: What determines Omega_2 beyond c3?")
    print("=" * 70)

    for n in [5, 6]:
        total = 2**(n*(n-1)//2)
        # Find c3 values where Omega_2 varies
        by_c3 = defaultdict(list)
        for bits in range(total):
            A = bits_to_adj(bits, n)
            c3 = 0
            for a, b, c in combinations(range(n), 3):
                if A[a][b] and A[b][c] and A[c][a]: c3 += 1
                if A[a][c] and A[c][b] and A[b][a]: c3 += 1
            a2 = enumerate_allowed_paths(A, n, 2)
            a1 = enumerate_allowed_paths(A, n, 1)
            od = compute_omega_dim(A, n, 2, a2, a1)

            # Score sequence
            scores = tuple(sorted([sum(A[i]) for i in range(n)]))

            by_c3[c3].append((od, scores, bits))

        # Report the ambiguous c3 values
        for c3_val in sorted(by_c3.keys()):
            entries = by_c3[c3_val]
            od_vals = sorted(set(e[0] for e in entries))
            if len(od_vals) > 1:
                print(f"\n  n={n}, c3={c3_val}: Omega_2 in {od_vals}")
                for od_val in od_vals:
                    subset = [e for e in entries if e[0] == od_val]
                    score_dist = defaultdict(int)
                    for _, score, _ in subset:
                        score_dist[score] += 1
                    top_scores = sorted(score_dist.items(), key=lambda x: -x[1])[:3]
                    print(f"    Omega_2={od_val}: {len(subset)} tours, "
                          f"top scores: {[(list(s), c) for s, c in top_scores]}")

    # Part 3: Is Omega_2 determined by (c3, score_sequence)?
    print("\n" + "=" * 70)
    print("Part 3: Is Omega_2 determined by (c3, score)?")
    print("=" * 70)

    for n in [5, 6]:
        total = 2**(n*(n-1)//2)
        by_key = defaultdict(list)
        for bits in range(total):
            A = bits_to_adj(bits, n)
            c3 = 0
            for a, b, c in combinations(range(n), 3):
                if A[a][b] and A[b][c] and A[c][a]: c3 += 1
                if A[a][c] and A[c][b] and A[b][a]: c3 += 1
            a2 = enumerate_allowed_paths(A, n, 2)
            a1 = enumerate_allowed_paths(A, n, 1)
            od = compute_omega_dim(A, n, 2, a2, a1)
            scores = tuple(sorted([sum(A[i]) for i in range(n)]))
            key = (c3, scores)
            by_key[key].append(od)

        deterministic = True
        for key in sorted(by_key.keys()):
            od_vals = sorted(set(by_key[key]))
            if len(od_vals) > 1:
                deterministic = False
                print(f"  n={n}, (c3={key[0]}, score={list(key[1])}): Omega_2 in {od_vals}")

        if deterministic:
            print(f"  n={n}: Omega_2 IS determined by (c3, score)")
        else:
            print(f"  n={n}: Omega_2 is NOT determined by (c3, score)")

if __name__ == '__main__':
    main()
