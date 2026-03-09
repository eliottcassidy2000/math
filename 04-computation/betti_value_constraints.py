"""
betti_value_constraints.py — What values can each beta_p take?

Known:
  beta_0 = 1 always (connectivity)
  beta_1 in {0, 1} always (rank(d_2) takes exactly 2 values)
  beta_2 = 0 always (proved)
  beta_3 = ? (gap always 1 in tested cases)
  beta_4 = ? (values 1,2,5,6 observed)
  beta_5 = ? (values 0,1 observed)

Questions:
  1. Is beta_3 in {0, 1} always?
  2. What values can beta_4 take?
  3. Is there a pattern: beta_{2k-1} in {0, 1} for all k?

Author: kind-pasteur-S45 (2026-03-09)
"""
import sys
import numpy as np
from math import comb
from itertools import combinations
from collections import defaultdict, Counter
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


def main():
    print("=" * 70)
    print("BETTI VALUE CONSTRAINTS")
    print("=" * 70)

    # Part 1: Exhaustive at n=5,6
    for n in [5, 6]:
        N = 2**(n*(n-1)//2)
        betti_vals = defaultdict(set)

        print(f"\n--- n={n}: Exhaustive ({N} tournaments) ---")

        for bits in range(N):
            A = bits_to_adj(bits, n)
            bv = compute_betti_vector(A, n)
            for p, bp in enumerate(bv):
                betti_vals[p].add(bp)

        print(f"  Value sets:")
        for p in range(n):
            vals = sorted(betti_vals[p])
            print(f"    beta_{p}: {vals}")

    # Part 2: Sampled at n=7,8
    for n in [7, 8]:
        rng = np.random.RandomState(42 + n)
        N = 500 if n == 7 else 200
        betti_vals = defaultdict(set)

        print(f"\n--- n={n}: Sampled ({N} tournaments) ---")

        for trial in range(N):
            A = random_tournament(n, rng)
            bv = compute_betti_vector(A, n)
            for p, bp in enumerate(bv):
                betti_vals[p].add(bp)

            if (trial + 1) % 100 == 0:
                print(f"  {trial+1}/{N} done", flush=True)

        print(f"  Value sets:")
        for p in range(n):
            vals = sorted(betti_vals[p])
            print(f"    beta_{p}: {vals}")

    # Part 3: Focus on beta_3 values at n=7,8
    # Is beta_3 always in {0, 1}?
    print("\n--- Part 3: beta_3 value distribution ---")
    for n in [7, 8]:
        rng = np.random.RandomState(99 + n)
        N = 1000 if n == 7 else 300
        b3_counter = Counter()

        for trial in range(N):
            A = random_tournament(n, rng)
            bv = compute_betti_vector(A, n)
            b3_counter[bv[3]] += 1

            if (trial + 1) % 200 == 0:
                print(f"  n={n}: {trial+1}/{N} done", flush=True)

        print(f"  n={n}: beta_3 distribution:")
        for val in sorted(b3_counter.keys()):
            cnt = b3_counter[val]
            print(f"    beta_3={val}: {cnt}/{N} ({100*cnt/N:.1f}%)")

    print("\nDONE.")


if __name__ == '__main__':
    main()
