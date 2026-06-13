"""
Filling ratio formula: dim(Omega_p) / C(n, p+1) for random tournaments.

The filling ratio f_p = <dim(Omega_p)> / C(n, p+1) measures how cyclic
content inflates the chain complex beyond the simplicial baseline.

From the crossover data:
  f_p appears to depend on p/(n-1) but also grows with n at fixed t=p/(n-1).

Key observations from data:
  - f_0 = 1 always (vertices)
  - f_1 = 1 always (all edges are tournament edges)
  - f_p > 1 for p >= 3 at n >= 6
  - f_{n-1} grows roughly as n/2 (ratio at top dimension)

E[|A_p|] = C(n,p+1) * (p+1)!/2^p for random tournament.
This is because E[Ham paths on k vertices] = k!/2^{k-1}.

So E[|A_p|]/C(n,p+1) = (p+1)!/2^p — the PATH MULTIPLICITY.

But dim(Omega_p) is NOT |A_p|. The Omega constraint kills some of these.
The question is: what fraction of A_p survives into Omega_p?

Survival fraction s_p = dim(Omega_p) / |A_p|
f_p = s_p * (p+1)!/2^p

Let's compute s_p directly.
"""
import numpy as np
from math import comb, factorial
from collections import defaultdict

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
    print("FILLING RATIO DECOMPOSITION: f_p = s_p * (p+1)!/2^p")
    print("s_p = dim(Omega_p) / |A_p|  (survival fraction)")
    print("=" * 70)

    # Part 1: Compute s_p for each n
    for n in range(3, 9):
        rng = np.random.RandomState(42 + n)
        N = min(200, 2**(n*(n-1)//2))

        Ap_sum = defaultdict(float)
        Om_sum = defaultdict(float)
        cnt = 0

        for _ in range(N):
            A = random_tournament(n, rng)
            for p in range(n):
                ap = enumerate_allowed_paths(A, n, p)
                apm1 = enumerate_allowed_paths(A, n, p-1) if p > 0 else []
                od = compute_omega_dim(A, n, p, ap, apm1)
                Ap_sum[p] += len(ap)
                Om_sum[p] += od
            cnt += 1

        print(f"\n  n={n} (d={n-1}), {cnt} samples:")
        print(f"  {'p':>3} | {'|A_p|':>8} | {'Omega_p':>8} | {'C(n,p+1)':>8} | {'(p+1)!/2^p':>10} | {'s_p':>8} | {'f_p':>8} | {'s_p theory':>10}")
        print(f"  {'-'*3}-+-{'-'*8}-+-{'-'*8}-+-{'-'*8}-+-{'-'*10}-+-{'-'*8}-+-{'-'*8}-+-{'-'*10}")
        for p in range(n):
            c = comb(n, p+1)
            path_mult = factorial(p+1) / (2**p)
            avg_Ap = Ap_sum[p] / cnt
            avg_Om = Om_sum[p] / cnt
            s_p = avg_Om / avg_Ap if avg_Ap > 0 else 0
            f_p = avg_Om / c if c > 0 else 0
            # s_p_theory: what would s_p be if Omega_p = C(n,p+1)?
            s_theory = c / (c * path_mult) if path_mult > 0 else 0  # = 1/path_mult = 2^p/(p+1)!
            print(f"  {p:3d} | {avg_Ap:8.1f} | {avg_Om:8.1f} | {c:8d} | {path_mult:10.2f} | {s_p:8.4f} | {f_p:8.4f} | {s_theory:10.4f}")

    # Part 2: Is there a formula for s_p?
    # For transitive: s_p = 1 and |A_p| = C(n,p+1), so trivially s_p = 1.
    # For random: s_p = dim(Omega_p) / |A_p|
    #   At p=1: s_1 = C(n,2) / C(n,2) = 1 (all edges survive)
    #   At p=2: s_2 varies
    #   At high p: s_p * path_mult = f_p > 1

    print("\n" + "=" * 70)
    print("Part 2: Normalized filling — f_p * 2^p / (p+1)! = s_p * C(n,p+1) / C(n,p+1) = s_p")
    print("If s_p ~ 2^p / (p+1)! then f_p ~ 1 (simplex-like)")
    print("If s_p > 2^p / (p+1)! then f_p > 1 (inflated)")
    print("=" * 70)

    # Part 3: Track how NON-ALLOWED constraints scale
    # The number of non-allowed (p-1)-paths determines how many constraints
    # are imposed on Omega_p. More constraints => smaller Omega_p.
    print("\n" + "=" * 70)
    print("Part 3: Non-allowed face count (constraints on Omega_p)")
    print("=" * 70)

    for n in range(3, 9):
        rng = np.random.RandomState(42 + n)
        N = min(200, 2**(n*(n-1)//2))

        na_sum = defaultdict(float)
        cnt = 0

        for _ in range(N):
            A = random_tournament(n, rng)
            for p in range(1, n):
                ap = enumerate_allowed_paths(A, n, p)
                apm1 = enumerate_allowed_paths(A, n, p-1)
                apm1_set = set(apm1)
                na_faces = set()
                for path in ap:
                    for sign, face in boundary_coeffs(path):
                        if len(set(face)) == len(face) and face not in apm1_set:
                            na_faces.add(face)
                na_sum[p] += len(na_faces)
            cnt += 1

        print(f"\n  n={n}:")
        for p in range(1, n):
            avg_na = na_sum[p] / cnt
            # Compare to total possible (p-1)-paths on n vertices: n!/(n-p)!
            total_possible = factorial(n) // factorial(n - p) if p <= n else 0
            c_pm1 = comb(n, p)
            avg_Ap_m1 = c_pm1 * factorial(p) / (2**(p-1))  # expected |A_{p-1}|
            na_ratio = avg_na / total_possible if total_possible > 0 else 0
            # Non-allowed = those (p-1)-paths that are NOT allowed but appear as faces
            print(f"    p={p}: avg non-allowed faces = {avg_na:.1f}, "
                  f"total possible (p-1)-seqs = {total_possible}, "
                  f"na/total = {na_ratio:.4f}")

    # Part 4: The filling ratio should relate to the CYCLIC CONTENT
    # For each (p+1)-vertex set, the number of Ham paths is:
    # = (p+1)! / 2^{C(p+1,2)} * sum over all orientations
    # For non-cyclic subtournaments: = 1 path
    # For cyclic: >= 2 paths
    # The EXCESS paths come from the cyclic structure
    print("\n" + "=" * 70)
    print("Part 4: Cyclic content — paths per vertex set distribution")
    print("=" * 70)

    for n in [5, 7, 8]:
        rng = np.random.RandomState(42 + n)
        N = 50

        for p in [2, 3, min(4, n-1)]:
            ppvs_total = defaultdict(float)  # paths-per-vertex-set distribution
            cnt = 0
            for _ in range(N):
                A = random_tournament(n, rng)
                paths = enumerate_allowed_paths(A, n, p)
                by_vset = defaultdict(int)
                for path in paths:
                    by_vset[frozenset(path)] += 1
                for v_count in by_vset.values():
                    ppvs_total[v_count] += 1.0 / N
                cnt += 1

            print(f"\n  n={n}, p={p}: avg paths-per-vertex-set distribution:")
            for k in sorted(ppvs_total.keys()):
                print(f"    {k} paths: {ppvs_total[k]:.1f} vertex sets")

if __name__ == '__main__':
    main()
