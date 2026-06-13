"""
beta34_coexistence.py — Investigate beta_3 + beta_4 coexistence at n=8.

Profile [1,0,0,1,1,0,0,0] is very rare at n=8.
When beta_3 and beta_4 are both nonzero, what structure does the tournament have?

Also check: does beta_3+beta_4 imply beta_5=0?
And: what's the MAXIMUM chi observed?

Author: kind-pasteur-S45 (2026-03-09)
"""
import sys
import numpy as np
from math import comb
from itertools import combinations
from collections import defaultdict, Counter
sys.stdout.reconfigure(line_buffering=True)

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

def compute_betti_vector(A, n, max_p=None):
    if max_p is None:
        max_p = n - 1

    ap = {}
    for p in range(max_p + 2):
        if p <= n - 1:
            ap[p] = enumerate_allowed_paths(A, n, p)
        else:
            ap[p] = []

    omega_bases = {}
    omega_dims = {}
    for p in range(max_p + 2):
        if p > n - 1 or len(ap[p]) == 0:
            omega_bases[p] = np.zeros((0, 0))
            omega_dims[p] = 0
            continue
        if p == 0:
            omega_bases[p] = np.eye(len(ap[p]))
            omega_dims[p] = len(ap[p])
            continue

        apm1_set = set(ap[p-1])
        non_allowed = {}
        na_count = 0
        for j, path in enumerate(ap[p]):
            for sign, face in boundary_coeffs(path):
                if len(set(face)) == len(face) and face not in apm1_set:
                    if face not in non_allowed:
                        non_allowed[face] = na_count
                        na_count += 1
        if na_count == 0:
            omega_bases[p] = np.eye(len(ap[p]))
            omega_dims[p] = len(ap[p])
        else:
            P = np.zeros((na_count, len(ap[p])))
            for j, path in enumerate(ap[p]):
                for sign, face in boundary_coeffs(path):
                    if face in non_allowed:
                        P[non_allowed[face], j] += sign
            U, S, Vt = np.linalg.svd(P, full_matrices=True)
            rank = int(sum(s > 1e-10 for s in S))
            ns = Vt[rank:].T
            if ns.shape[1] > 0:
                omega_bases[p] = ns
                omega_dims[p] = ns.shape[1]
            else:
                omega_bases[p] = np.zeros((len(ap[p]), 0))
                omega_dims[p] = 0

    betti = []
    for p in range(max_p + 1):
        dim_p = omega_dims.get(p, 0)
        # rank(d_p)
        if p == 0 or dim_p == 0:
            rank_dp = 0
        else:
            bd = np.zeros((len(ap[p-1]), len(ap[p])))
            idx_pm1 = {path: i for i, path in enumerate(ap[p-1])}
            for j, path in enumerate(ap[p]):
                for sign, face in boundary_coeffs(path):
                    if face in idx_pm1:
                        bd[idx_pm1[face], j] += sign
            dp_om = bd @ omega_bases[p]
            sv = np.linalg.svd(dp_om, compute_uv=False)
            rank_dp = int(sum(s > 1e-8 for s in sv))

        ker_dp = dim_p - rank_dp

        # im(d_{p+1})
        dim_p1 = omega_dims.get(p+1, 0)
        if p >= n - 1 or p >= max_p or dim_p1 == 0:
            im_dp1 = 0
        else:
            bd1 = np.zeros((len(ap[p]), len(ap[p+1])))
            idx_p = {path: i for i, path in enumerate(ap[p])}
            for j, path in enumerate(ap[p+1]):
                for sign, face in boundary_coeffs(path):
                    if face in idx_p:
                        bd1[idx_p[face], j] += sign
            dp1_om = bd1 @ omega_bases[p+1]
            sv1 = np.linalg.svd(dp1_om, compute_uv=False)
            im_dp1 = int(sum(s > 1e-8 for s in sv1))

        betti.append(ker_dp - im_dp1)

    return betti

def count_3cycles(A, n):
    c3 = 0
    for i, j, k in combinations(range(n), 3):
        if (A[i][j] and A[j][k] and A[k][i]) or (A[i][k] and A[k][j] and A[j][i]):
            c3 += 1
    return c3


def main():
    print("=" * 70)
    print("BETA_3 + BETA_4 COEXISTENCE INVESTIGATION")
    print("=" * 70)

    n = 8
    rng = np.random.RandomState(42)
    N = 2000

    profiles = Counter()
    interesting = []  # beta_3 AND beta_4 both nonzero

    for trial in range(N):
        A = random_tournament(n, rng)
        bv = compute_betti_vector(A, n)
        profiles[tuple(bv)] += 1

        if bv[3] > 0 and bv[4] > 0:
            scores = tuple(sorted([int(sum(A[i])) for i in range(n)]))
            c3 = count_3cycles(A, n)
            chi = sum((-1)**p * bv[p] for p in range(n))
            interesting.append({
                'trial': trial,
                'betti': bv,
                'scores': scores,
                'c3': c3,
                'chi': chi,
                'A': A.copy()
            })

        if (trial + 1) % 500 == 0:
            print(f"  {trial+1}/{N}: {len(interesting)} coexistence cases found", flush=True)

    print(f"\n  Total beta_3+beta_4 coexistence: {len(interesting)}/{N} "
          f"({100*len(interesting)/N:.2f}%)")

    print(f"\n  All Betti profiles ({len(profiles)} distinct):")
    for profile in sorted(profiles.keys(), key=lambda x: -profiles[x]):
        cnt = profiles[profile]
        chi = sum((-1)**p * profile[p] for p in range(n))
        print(f"    {list(profile)}: {cnt} ({100*cnt/N:.1f}%), chi={chi}")

    if interesting:
        print(f"\n  Coexistence cases detail:")
        for case in interesting[:10]:
            print(f"    trial {case['trial']}: betti={case['betti']}, "
                  f"scores={case['scores']}, c3={case['c3']}, chi={case['chi']}")

        # Check: is there a score pattern?
        score_counts = Counter(case['scores'] for case in interesting)
        print(f"\n  Score distributions for coexistence:")
        for sc, cnt in score_counts.most_common():
            print(f"    {sc}: {cnt}")

    # Part 2: Focus on the STRONGEST seesaw claim
    print("\n--- Part 2: Adjacent-odd seesaw test on FULL 2000-sample set ---")
    adj_odd_violations = {(1,3): 0, (3,5): 0, (5,7): 0}
    for profile, cnt in profiles.items():
        for i, j in adj_odd_violations:
            if profile[i] > 0 and profile[j] > 0:
                adj_odd_violations[(i,j)] += cnt
                print(f"  VIOLATION beta_{i}*beta_{j}: {list(profile)} x {cnt}")

    for pair, cnt in adj_odd_violations.items():
        print(f"  beta_{pair[0]}*beta_{pair[1]} = 0: {N-cnt}/{N} satisfied "
              f"({100*(N-cnt)/N:.1f}%)")

    # Part 3: Does beta_3+beta_4 imply beta_5=0?
    print("\n--- Part 3: beta_3 + beta_4 implies beta_5 = 0? ---")
    if interesting:
        b5_vals = [case['betti'][5] for case in interesting]
        print(f"  beta_5 values when beta_3+beta_4>0: {sorted(set(b5_vals))}")
        if all(v == 0 for v in b5_vals):
            print(f"  YES: beta_5 = 0 in all {len(interesting)} coexistence cases")
        else:
            print(f"  NO: some have beta_5 > 0")

    # Part 4: Non-adjacent odd seesaw
    print("\n--- Part 4: Non-adjacent pairs ---")
    for i, j in [(1,5), (1,7), (3,7)]:
        cnt = 0
        for profile, pcnt in profiles.items():
            if profile[i] > 0 and profile[j] > 0:
                cnt += pcnt
        print(f"  beta_{i}*beta_{j}: {cnt}/{N} coexist")

    # Part 5: Odd-even pairs
    print("\n--- Part 5: Odd-even coexistence ---")
    for i in [1, 3, 5, 7]:
        for j in [2, 4, 6]:
            cnt = 0
            for profile, pcnt in profiles.items():
                if profile[i] > 0 and profile[j] > 0:
                    cnt += pcnt
            if cnt > 0:
                print(f"  beta_{i}*beta_{j}: {cnt}/{N} coexist ({100*cnt/N:.2f}%)")

    print("\nDONE.")


if __name__ == '__main__':
    main()
