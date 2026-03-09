"""
defect_rate_ushape.py — Investigate the "U-shape" in defect rates P(beta_p > 0).

S44 observed that defect rates form a U-shape across p: high at extremes, low in middle.
THM-095 proves beta_1 * beta_3 = 0 (seesaw), so there's structure in the defect pattern.

Plan:
  1. Compute full Betti vectors at n=6 (exhaustive) and n=7,8 (sampled)
  2. Track which (p,n) pairs have beta_p > 0
  3. Look for patterns: is the defect rate symmetric? monotone in p?
  4. Relate to filling ratio and cycle counts

Also: investigate WHY certain beta_p are always 0.
  - beta_0 = 1 always (connectivity)
  - beta_2 = 0 always (THM-108+109)
  - beta_4 = 0 at n<=7, first nonzero at n=8
  - beta_even = 0? (conjecture from S44)

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
    """Compute full Betti vector beta_0, ..., beta_{n-1}."""
    if max_p is None:
        max_p = n - 1

    # Enumerate all allowed paths
    ap = {}
    for p in range(max_p + 2):  # need p+1 for im(d_{p+1})
        if p <= n - 1:
            ap[p] = enumerate_allowed_paths(A, n, p)
        else:
            ap[p] = []

    # Compute Omega bases
    omega_bases = {}
    for p in range(max_p + 2):
        if p > n - 1 or len(ap[p]) == 0:
            omega_bases[p] = np.zeros((0, 0))
            continue
        if p == 0:
            omega_bases[p] = np.eye(len(ap[p]))
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
        else:
            P = np.zeros((na_count, len(ap[p])))
            for j, path in enumerate(ap[p]):
                for sign, face in boundary_coeffs(path):
                    if face in non_allowed:
                        P[non_allowed[face], j] += sign
            U, S, Vt = np.linalg.svd(P, full_matrices=True)
            rank = int(sum(s > 1e-10 for s in S))
            ns = Vt[rank:].T
            omega_bases[p] = ns if ns.shape[1] > 0 else np.zeros((len(ap[p]), 0))

    # Compute Betti numbers
    betti = []
    for p in range(max_p + 1):
        dim_p = omega_bases[p].shape[1] if omega_bases[p].ndim == 2 and omega_bases[p].shape[0] > 0 else 0

        # ker(d_p)
        if p == 0 or dim_p == 0:
            ker_dp = dim_p
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
        if p >= n - 1 or p >= max_p:
            im_dp1 = 0
        else:
            dim_p1 = omega_bases[p+1].shape[1] if omega_bases[p+1].ndim == 2 and omega_bases[p+1].shape[0] > 0 else 0
            if dim_p1 == 0:
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

        bp = ker_dp - im_dp1
        betti.append(bp)

    return betti


def main():
    print("=" * 70)
    print("DEFECT RATE U-SHAPE INVESTIGATION")
    print("=" * 70)

    # Part 1: Exhaustive at n=5,6
    for n in [5, 6]:
        N = 2**(n*(n-1)//2)
        betti_profiles = Counter()
        betti_nonzero = defaultdict(int)
        total = 0

        print(f"\n--- n={n}: Exhaustive ({N} tournaments) ---", flush=True)

        for bits in range(N):
            A = bits_to_adj(bits, n)
            bv = compute_betti_vector(A, n)
            betti_profiles[tuple(bv)] += 1
            for p, bp in enumerate(bv):
                if bp > 0:
                    betti_nonzero[p] += 1
            total += 1

        print(f"  Defect rates:")
        for p in range(n):
            rate = betti_nonzero[p] / total
            print(f"    beta_{p}: {rate:.4f} ({betti_nonzero[p]}/{total})")

        print(f"\n  Distinct Betti profiles: {len(betti_profiles)}")
        for profile in sorted(betti_profiles.keys(), key=lambda x: -betti_profiles[x]):
            cnt = betti_profiles[profile]
            if cnt >= total * 0.001 or len(betti_profiles) <= 20:
                print(f"    {list(profile)}: {cnt} ({100*cnt/total:.1f}%)")

    # Part 2: Sampled at n=7
    for n in [7]:
        N = 500
        rng = np.random.RandomState(42 + n)
        betti_profiles = Counter()
        betti_nonzero = defaultdict(int)

        print(f"\n--- n={n}: Sampled ({N} tournaments) ---", flush=True)

        for trial in range(N):
            A = random_tournament(n, rng)
            bv = compute_betti_vector(A, n)
            betti_profiles[tuple(bv)] += 1
            for p, bp in enumerate(bv):
                if bp > 0:
                    betti_nonzero[p] += 1

            if (trial + 1) % 100 == 0:
                print(f"  {trial+1}/{N} done", flush=True)

        print(f"\n  Defect rates:")
        for p in range(n):
            rate = betti_nonzero[p] / N
            print(f"    beta_{p}: {rate:.4f} ({betti_nonzero[p]}/{N})")

        print(f"\n  Distinct Betti profiles: {len(betti_profiles)}")
        for profile in sorted(betti_profiles.keys(), key=lambda x: -betti_profiles[x]):
            cnt = betti_profiles[profile]
            if cnt >= 2:
                print(f"    {list(profile)}: {cnt} ({100*cnt/N:.1f}%)")

    # Part 3: beta_even vanishing check at n=7
    print("\n--- Part 3: Even Betti vanishing at n=7 ---")
    print(f"  beta_0 nonzero: {betti_nonzero[0]} (always 1, expected)")
    print(f"  beta_2 nonzero: {betti_nonzero[2]} (THM-108+109 proves 0)")
    print(f"  beta_4 nonzero: {betti_nonzero[4]} (0 at n<=7?)")
    print(f"  beta_6 nonzero: {betti_nonzero[6]} (top dimension)")

    # Part 4: Seesaw check: beta_1 * beta_3 = 0 (THM-095)
    print("\n--- Part 4: Seesaw beta_1 * beta_3 = 0 ---")
    violations = 0
    for profile, cnt in betti_profiles.items():
        if profile[1] > 0 and profile[3] > 0:
            violations += cnt
            print(f"  VIOLATION: {list(profile)} x {cnt}")
    print(f"  Total seesaw violations: {violations}/{N}")

    # Part 5: n=8 partial (just check which beta_p can be nonzero)
    print(f"\n--- n=8: Sampled (100 tournaments, checking nonzero beta) ---", flush=True)
    n = 8
    rng = np.random.RandomState(42 + n)
    N = 100
    betti_nonzero_8 = defaultdict(int)
    betti_profiles_8 = Counter()

    for trial in range(N):
        A = random_tournament(n, rng)
        bv = compute_betti_vector(A, n)
        betti_profiles_8[tuple(bv)] += 1
        for p, bp in enumerate(bv):
            if bp > 0:
                betti_nonzero_8[p] += 1

        if (trial + 1) % 25 == 0:
            print(f"  {trial+1}/{N} done", flush=True)

    print(f"\n  n=8 Defect rates:")
    for p in range(n):
        rate = betti_nonzero_8[p] / N
        print(f"    beta_{p}: {rate:.4f} ({betti_nonzero_8[p]}/{N})")

    print(f"\n  Even Betti check:")
    print(f"    beta_2: {betti_nonzero_8[2]}/{N}")
    print(f"    beta_4: {betti_nonzero_8[4]}/{N}")
    print(f"    beta_6: {betti_nonzero_8[6]}/{N}")

    print(f"\n  Seesaw checks:")
    for profile, cnt in betti_profiles_8.items():
        if profile[1] > 0 and profile[3] > 0:
            print(f"    beta_1*beta_3 VIOLATION: {list(profile)} x {cnt}")
        if profile[3] > 0 and profile[5] > 0:
            print(f"    beta_3*beta_5 VIOLATION: {list(profile)} x {cnt}")
        if profile[5] > 0 and profile[7] > 0:
            print(f"    beta_5*beta_7 VIOLATION: {list(profile)} x {cnt}")

    print(f"\n  Distinct Betti profiles at n=8: {len(betti_profiles_8)}")
    for profile in sorted(betti_profiles_8.keys(), key=lambda x: -betti_profiles_8[x]):
        cnt = betti_profiles_8[profile]
        if cnt >= 1:
            print(f"    {list(profile)}: {cnt} ({100*cnt/N:.1f}%)")

    print("\nDONE.")


if __name__ == '__main__':
    main()
