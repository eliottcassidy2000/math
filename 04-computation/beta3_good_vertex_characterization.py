"""
beta3_good_vertex_characterization.py — What makes a vertex "good" for beta_3?

For the beta_3 <= 1 proof, we need: for every T with beta_3(T) > 0,
exists v with beta_3(T\v) = 0 ("good vertex for beta_3").

At n=6: ALL vertices are good (completely fragile). Why?
At n=7: most but not all vertices are good. What distinguishes good from bad?

Key quantities to check:
- Score (out-degree) of v
- Number of 3-cycles through v (c3(v))
- Position in topological/score ordering
- Whether v is in a 3-cycle (cycle vertex vs non-cycle vertex)
- Local structure: in-neighbors, out-neighbors connectivity

Author: kind-pasteur-S46 (2026-03-09)
"""
import sys
import numpy as np
from math import comb
from itertools import combinations
from collections import Counter, defaultdict
sys.stdout.reconfigure(line_buffering=True)

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

def compute_omega_basis(ap, p):
    if not ap.get(p, []):
        return np.zeros((0, 0)), 0
    if p == 0:
        return np.eye(len(ap[p])), len(ap[p])

    apm1_set = set(ap.get(p-1, []))
    non_allowed = {}
    na_count = 0
    for j, path in enumerate(ap[p]):
        for sign, face in boundary_coeffs(path):
            if len(set(face)) == len(face) and face not in apm1_set:
                if face not in non_allowed:
                    non_allowed[face] = na_count
                    na_count += 1
    if na_count == 0:
        return np.eye(len(ap[p])), len(ap[p])

    P = np.zeros((na_count, len(ap[p])))
    for j, path in enumerate(ap[p]):
        for sign, face in boundary_coeffs(path):
            if face in non_allowed:
                P[non_allowed[face], j] += sign
    U, S, Vt = np.linalg.svd(P, full_matrices=True)
    rank = int(sum(s > 1e-10 for s in S))
    ns = Vt[rank:].T
    if ns.shape[1] > 0:
        return ns, ns.shape[1]
    else:
        return np.zeros((len(ap[p]), 0)), 0


def compute_beta3(A, n):
    ap = {}
    for p in range(min(6, n)):
        ap[p] = enumerate_allowed_paths(A, n, p)

    omega_bases = {}
    omega_dims = {}
    for p in range(min(6, n)):
        if not ap.get(p, []):
            omega_bases[p] = np.zeros((0, 0))
            omega_dims[p] = 0
        else:
            omega_bases[p], omega_dims[p] = compute_omega_basis(ap, p)

    dim3 = omega_dims.get(3, 0)
    if dim3 == 0:
        return 0

    bd3 = np.zeros((len(ap.get(2, [])), len(ap[3])))
    idx2 = {path: i for i, path in enumerate(ap.get(2, []))}
    for j, path in enumerate(ap[3]):
        for sign, face in boundary_coeffs(path):
            if face in idx2:
                bd3[idx2[face], j] += sign

    O3 = omega_bases[3]
    d3_om = bd3 @ O3
    sv3 = np.linalg.svd(d3_om, compute_uv=False)
    rank_d3 = int(sum(s > 1e-8 for s in sv3))
    ker_d3 = dim3 - rank_d3

    if ker_d3 == 0:
        return 0

    dim4 = omega_dims.get(4, 0)
    if dim4 == 0:
        return ker_d3

    bd4 = np.zeros((len(ap[3]), len(ap.get(4, []))))
    idx3 = {path: i for i, path in enumerate(ap[3])}
    for j, path in enumerate(ap.get(4, [])):
        for sign, face in boundary_coeffs(path):
            if face in idx3:
                bd4[idx3[face], j] += sign

    O4 = omega_bases[4]
    d4_om = bd4 @ O4
    O3pinv = np.linalg.pinv(O3)
    d4_omega3 = O3pinv @ d4_om
    sv4 = np.linalg.svd(d4_omega3, compute_uv=False)
    rank_d4 = int(sum(s > 1e-8 for s in sv4))

    return ker_d3 - rank_d4


def deletion_adj(A, n, v):
    vertices = [i for i in range(n) if i != v]
    n1 = len(vertices)
    A1 = np.zeros((n1, n1), dtype=int)
    for i, vi in enumerate(vertices):
        for j, vj in enumerate(vertices):
            A1[i][j] = A[vi][vj]
    return A1, n1


def count_c3_through_v(A, n, v):
    """Count 3-cycles through vertex v."""
    c = 0
    for i in range(n):
        if i == v: continue
        for j in range(n):
            if j == v or j == i: continue
            if (A[v][i] == 1 and A[i][j] == 1 and A[j][v] == 1) or \
               (A[v][j] == 1 and A[j][i] == 1 and A[i][v] == 1):
                c += 1
    return c // 2  # Each 3-cycle counted twice


def is_strongly_connected(A, n):
    if n <= 1: return True
    visited = set([0])
    queue = [0]
    while queue:
        v = queue.pop(0)
        for u in range(n):
            if A[v][u] == 1 and u not in visited:
                visited.add(u)
                queue.append(u)
    if len(visited) < n: return False
    visited = set([0])
    queue = [0]
    while queue:
        v = queue.pop(0)
        for u in range(n):
            if A[u][v] == 1 and u not in visited:
                visited.add(u)
                queue.append(u)
    return len(visited) == n


def main():
    print("=" * 70)
    print("BETA_3 GOOD VERTEX CHARACTERIZATION")
    print("=" * 70)

    # Part 1: At n=7, characterize good vs bad vertices for beta_3
    print("\n--- Part 1: Good vs bad vertex properties at n=7 ---")
    n = 7
    rng = np.random.RandomState(42)
    N = 500
    good_scores = Counter()
    bad_scores = Counter()
    good_c3 = Counter()
    bad_c3 = Counter()
    good_sc = Counter()
    bad_sc = Counter()
    b3_count = 0

    for trial in range(N):
        A = random_tournament(n, rng)
        b3 = compute_beta3(A, n)
        if b3 == 0:
            continue
        b3_count += 1

        scores = [int(sum(A[i])) for i in range(n)]
        c3_per_v = [count_c3_through_v(A, n, v) for v in range(n)]

        for v in range(n):
            A1, n1 = deletion_adj(A, n, v)
            b3v = compute_beta3(A1, n1)

            sc = is_strongly_connected(A1, n1)

            if b3v == 0:  # good
                good_scores[scores[v]] += 1
                good_c3[c3_per_v[v]] += 1
                good_sc[sc] += 1
            else:  # bad
                bad_scores[scores[v]] += 1
                bad_c3[c3_per_v[v]] += 1
                bad_sc[sc] += 1

        if (trial + 1) % 100 == 0:
            print(f"  n=7: {trial+1}/{N} done, {b3_count} with beta_3>0", flush=True)

    print(f"\n  beta_3>0: {b3_count}")
    print(f"\n  Good vertex (beta_3(T\\v)=0) score distribution: {dict(sorted(good_scores.items()))}")
    print(f"  Bad vertex (beta_3(T\\v)>0) score distribution: {dict(sorted(bad_scores.items()))}")
    print(f"\n  Good vertex c3(v) distribution: {dict(sorted(good_c3.items()))}")
    print(f"  Bad vertex c3(v) distribution: {dict(sorted(bad_c3.items()))}")
    print(f"\n  Good vertex: T\\v SC? {dict(good_sc)}")
    print(f"  Bad vertex: T\\v SC? {dict(bad_sc)}")

    # Part 2: Score-based vertex selection rules
    print("\n--- Part 2: Score-based vertex selection rules at n=7 ---")
    rng2 = np.random.RandomState(123)
    N2 = 500
    rules_success = defaultdict(int)
    rules_total = 0
    b3_count2 = 0

    for trial in range(N2):
        A = random_tournament(n, rng2)
        b3 = compute_beta3(A, n)
        if b3 == 0:
            continue
        b3_count2 += 1
        rules_total += 1

        scores = [int(sum(A[i])) for i in range(n)]
        c3_per_v = [count_c3_through_v(A, n, v) for v in range(n)]

        # Compute beta_3 for all deletions
        b3_del = []
        for v in range(n):
            A1, n1 = deletion_adj(A, n, v)
            b3_del.append(compute_beta3(A1, n1))

        # Rule 1: Delete vertex with max score
        v_max_score = max(range(n), key=lambda v: scores[v])
        if b3_del[v_max_score] == 0:
            rules_success['max_score'] += 1

        # Rule 2: Delete vertex with min score
        v_min_score = min(range(n), key=lambda v: scores[v])
        if b3_del[v_min_score] == 0:
            rules_success['min_score'] += 1

        # Rule 3: Delete vertex with max c3(v)
        v_max_c3 = max(range(n), key=lambda v: c3_per_v[v])
        if b3_del[v_max_c3] == 0:
            rules_success['max_c3'] += 1

        # Rule 4: Delete vertex with min c3(v)
        v_min_c3 = min(range(n), key=lambda v: c3_per_v[v])
        if b3_del[v_min_c3] == 0:
            rules_success['min_c3'] += 1

        # Rule 5: Delete vertex with max or min score (whichever extreme)
        v_extreme = v_max_score if (max(scores) - (n-1)/2) >= ((n-1)/2 - min(scores)) else v_min_score
        if b3_del[v_extreme] == 0:
            rules_success['extreme_score'] += 1

        # Rule 6: Delete median-score vertex
        sorted_v = sorted(range(n), key=lambda v: scores[v])
        v_median = sorted_v[n // 2]
        if b3_del[v_median] == 0:
            rules_success['median_score'] += 1

        # Rule 7: Delete the vertex that maximizes the MINIMUM deletion score
        # (Most "disruptive" vertex — removes the most score imbalance)
        if (trial + 1) % 100 == 0:
            print(f"  n=7: {trial+1}/{N2} done, {b3_count2} with beta_3>0", flush=True)

    print(f"\n  beta_3>0: {b3_count2}")
    print(f"  Vertex selection rule success rates:")
    for rule, count in sorted(rules_success.items(), key=lambda x: -x[1]):
        pct = 100 * count / rules_total if rules_total > 0 else 0
        print(f"    {rule}: {count}/{rules_total} = {pct:.1f}%")

    # Part 3: n=8 — same characterization
    print("\n--- Part 3: Good vs bad vertex properties at n=8 ---")
    n = 8
    rng3 = np.random.RandomState(77)
    N3 = 200
    good_scores8 = Counter()
    bad_scores8 = Counter()
    b3_count3 = 0

    for trial in range(N3):
        A = random_tournament(n, rng3)
        try:
            b3 = compute_beta3(A, n)
        except:
            continue
        if b3 == 0:
            continue
        b3_count3 += 1

        scores = [int(sum(A[i])) for i in range(n)]

        for v in range(n):
            A1, n1 = deletion_adj(A, n, v)
            b3v = compute_beta3(A1, n1)
            if b3v == 0:
                good_scores8[scores[v]] += 1
            else:
                bad_scores8[scores[v]] += 1

        if (trial + 1) % 50 == 0:
            print(f"  n=8: {trial+1}/{N3} done, {b3_count3} with beta_3>0", flush=True)

    print(f"\n  beta_3>0: {b3_count3}")
    print(f"  Good vertex score dist: {dict(sorted(good_scores8.items()))}")
    print(f"  Bad vertex score dist: {dict(sorted(bad_scores8.items()))}")

    # Part 4: Fraction of good vertices vs n
    print("\n--- Part 4: Good vertex fraction vs n ---")
    for n in [6, 7]:
        rng4 = np.random.RandomState(42)
        if n == 6:
            N4 = 2**(n*(n-1)//2)
            sample = False
        else:
            N4 = 500
            sample = True

        good_frac_data = []
        b3_c = 0

        for trial in range(N4):
            if sample:
                A = random_tournament(n, rng4)
            else:
                A = bits_to_adj(trial, n)

            b3 = compute_beta3(A, n)
            if b3 == 0:
                continue
            b3_c += 1

            num_good = sum(1 for v in range(n)
                          if compute_beta3(*deletion_adj(A, n, v)) == 0)
            good_frac_data.append(num_good / n)

            if not sample and (trial + 1) % 10000 == 0:
                print(f"  n={n}: {trial+1}/{N4}", flush=True)

        avg_frac = sum(good_frac_data) / len(good_frac_data) if good_frac_data else 0
        min_frac = min(good_frac_data) if good_frac_data else 0
        print(f"  n={n}: {b3_c} tours with beta_3>0, good fraction: avg={avg_frac:.3f}, min={min_frac:.3f}")

    print("\nDONE.")


if __name__ == '__main__':
    main()
