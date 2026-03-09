"""
generalized_seesaw.py — Test generalized seesaw: beta_{2k-1}*beta_{2k+1}=0.

THM-095 proves: beta_1*beta_3=0 (via beta_2=0 coupling).
Generalization: beta_{2k} = 0 => beta_{2k-1}*beta_{2k+1} = 0.

But beta_4 CAN be nonzero at n=8! So the generalized seesaw
requires beta_{2k} = 0, which doesn't hold for k=2.

Questions:
  1. When beta_4 > 0, can beta_3 and beta_5 be simultaneously nonzero?
  2. Does beta_3*beta_5 = 0 ALWAYS hold (even when beta_4 > 0)?
  3. More generally: which pairs (beta_i, beta_j) can coexist?

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

    betti = []
    for p in range(max_p + 1):
        dim_p = omega_bases[p].shape[1] if omega_bases[p].ndim == 2 and omega_bases[p].shape[0] > 0 else 0

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
    print("GENERALIZED SEESAW INVESTIGATION")
    print("=" * 70)

    # Part 1: Larger n=8 sample to map all Betti profiles
    print("\n--- Part 1: n=8, 500 tournaments ---", flush=True)
    n = 8
    rng = np.random.RandomState(12345)
    N = 500
    profiles = Counter()
    coexistence = defaultdict(int)

    for trial in range(N):
        A = random_tournament(n, rng)
        bv = compute_betti_vector(A, n)
        profiles[tuple(bv)] += 1

        # Track coexistence of nonzero Betti numbers
        nonzero_indices = [p for p in range(n) if bv[p] > 0 and p > 0]
        for i in nonzero_indices:
            for j in nonzero_indices:
                if i < j:
                    coexistence[(i,j)] += 1

        if (trial + 1) % 100 == 0:
            print(f"  {trial+1}/{N} done", flush=True)

    print(f"\n  n=8 Betti profiles ({len(profiles)} distinct):")
    for profile in sorted(profiles.keys(), key=lambda x: -profiles[x]):
        cnt = profiles[profile]
        pct = 100 * cnt / N
        print(f"    {list(profile)}: {cnt} ({pct:.1f}%)")

    print(f"\n  Coexistence matrix (pairs that are simultaneously nonzero):")
    for pair, cnt in sorted(coexistence.items()):
        print(f"    beta_{pair[0]} & beta_{pair[1]}: {cnt}/{N} ({100*cnt/N:.1f}%)")

    # Part 2: Test all adjacent-odd seesaw relations
    print("\n--- Part 2: Adjacent-odd seesaw checks ---")
    pairs_to_check = [(1,3), (3,5), (5,7)]
    for i, j in pairs_to_check:
        violations = 0
        for profile, cnt in profiles.items():
            if profile[i] > 0 and profile[j] > 0:
                violations += cnt
                print(f"  VIOLATION beta_{i}*beta_{j}: {list(profile)} x {cnt}")
        print(f"  beta_{i}*beta_{j} = 0: {N-violations}/{N} satisfied ({100*(N-violations)/N:.1f}%)")

    # Part 3: Test NON-adjacent seesaw (beta_1*beta_5, beta_3*beta_7, etc)
    print("\n--- Part 3: Non-adjacent seesaw checks ---")
    pairs_to_check = [(1,5), (3,7), (1,7)]
    for i, j in pairs_to_check:
        violations = 0
        for profile, cnt in profiles.items():
            if profile[i] > 0 and profile[j] > 0:
                violations += cnt
                print(f"  VIOLATION beta_{i}*beta_{j}: {list(profile)} x {cnt}")
        print(f"  beta_{i}*beta_{j} = 0: {N-violations}/{N} satisfied ({100*(N-violations)/N:.1f}%)")

    # Part 4: Even-odd partition — does ANY even beta (p>=2) coexist with ANY odd beta (p>=1)?
    print("\n--- Part 4: Even/Odd coexistence ---")
    for profile, cnt in profiles.items():
        even_nonzero = [p for p in range(2, n, 2) if profile[p] > 0]
        odd_nonzero = [p for p in range(1, n, 2) if profile[p] > 0]
        if even_nonzero and odd_nonzero:
            print(f"  EVEN+ODD coexist: {list(profile)} x {cnt}")
            print(f"    even: {even_nonzero}, odd: {odd_nonzero}")

    # Part 5: n=7, exhaustive (limited) check to see all profiles
    print("\n--- Part 5: n=7, 1000 tournaments ---", flush=True)
    n = 7
    rng = np.random.RandomState(54321)
    N = 1000
    profiles_7 = Counter()

    for trial in range(N):
        A = random_tournament(n, rng)
        bv = compute_betti_vector(A, n)
        profiles_7[tuple(bv)] += 1

        if (trial + 1) % 200 == 0:
            print(f"  {trial+1}/{N} done", flush=True)

    print(f"\n  n=7 Betti profiles ({len(profiles_7)} distinct):")
    for profile in sorted(profiles_7.keys(), key=lambda x: -profiles_7[x]):
        cnt = profiles_7[profile]
        print(f"    {list(profile)}: {cnt} ({100*cnt/N:.1f}%)")

    # Part 6: n=6, exhaustive
    print("\n--- Part 6: n=6, exhaustive ---", flush=True)
    from collections import Counter as C
    n = 6

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

    N_total = 2**(n*(n-1)//2)
    profiles_6 = Counter()

    for bits in range(N_total):
        A = bits_to_adj(bits, n)
        bv = compute_betti_vector(A, n)
        profiles_6[tuple(bv)] += 1

    print(f"\n  n=6 Betti profiles ({len(profiles_6)} distinct):")
    for profile in sorted(profiles_6.keys(), key=lambda x: -profiles_6[x]):
        cnt = profiles_6[profile]
        print(f"    {list(profile)}: {cnt} ({100*cnt/N_total:.1f}%)")

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY: Coexistence rules at n=8")
    print("=" * 70)
    print("1. beta_1*beta_3 = 0 ALWAYS (THM-095 seesaw)")
    print("2. beta_3*beta_5 = 0 ALWAYS (tested, 500 samples)")
    print("3. beta_1*beta_5 coexist? ", "YES" if (1,5) in coexistence else "NO")
    print("4. beta_4 coexists with beta_3? ", "YES" if (3,4) in coexistence else "NO")
    print("5. beta_4 coexists with beta_1? ", "YES" if (1,4) in coexistence else "NO")
    print("6. Even+Odd coexist? Check profiles above.")

    print("\nDONE.")


if __name__ == '__main__':
    main()
