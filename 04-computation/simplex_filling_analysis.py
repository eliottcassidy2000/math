"""
Simplex filling analysis: tournament orientation vs simplicial structure.

Key insight: Transitive tournament = contractible simplex.
- dim(Omega_p) = C(n, p+1) for transitive (one path per vertex set)
- For non-transitive: cyclic subtournaments add paths, changing Omega_p

The FILLING RATIO f_p(T) = dim(Omega_p) / C(n, p+1) measures how the
tournament orientation "inflates" or "deflates" the simplex at dimension p.

For transitive: f_p = 1 for all p (exact simplex)
For random: f_p depends on the density of cycles

Questions:
1. How does f_p scale with n for fixed p?
2. Is there a universal filling curve f(p/n)?
3. Does the filling ratio predict homological non-triviality?
4. What's the relationship between f_p and the number of directed cycles?

Also: Omega_p for transitive = Pascal's triangle row = C(n, p+1).
This means the path complex of the transitive tournament IS the simplicial
chain complex of Delta_{n-1}. The boundary map in path homology
COINCIDES with the simplicial boundary for the transitive case.
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

def ham_path_count_on_subset(A, n, subset):
    """Count directed Hamiltonian paths on subset of vertices."""
    sub = sorted(subset)
    k = len(sub)
    if k <= 1: return 1
    # Build sub-adjacency
    idx_map = {v: i for i, v in enumerate(sub)}
    sub_A = np.zeros((k, k), dtype=int)
    for i, u in enumerate(sub):
        for j, v in enumerate(sub):
            sub_A[i][j] = A[u][v]
    # Count Ham paths via bitmask DP
    dp = {}
    for v in range(k):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << k):
        bits = bin(mask).count('1')
        if bits == 1: continue
        for v in range(k):
            if not (mask & (1 << v)): continue
            prev_mask = mask ^ (1 << v)
            for u in range(k):
                if (prev_mask & (1 << u)) and sub_A[u][v]:
                    if (mask, v) not in dp:
                        dp[(mask, v)] = 0
                    dp[(mask, v)] += dp.get((prev_mask, u), 0)
    full = (1 << k) - 1
    return sum(dp.get((full, v), 0) for v in range(k))

def main():
    print("=" * 70)
    print("SIMPLEX FILLING ANALYSIS")
    print("Transitive = contractible simplex. Cycles = topological twists.")
    print("=" * 70)

    # Part 1: Average filling ratio for random tournaments
    print("\n--- Part 1: Average filling ratio f_p = <Omega_p> / C(n,p+1) ---")
    print("(Transitive: f_p = 1.0 for all p)")

    for n in range(3, 9):
        rng = np.random.RandomState(42 + n)
        N_samples = min(200, 2 ** (n*(n-1)//2))

        # Accumulate
        omega_sums = defaultdict(float)
        Ap_sums = defaultdict(float)
        count = 0

        for _ in range(N_samples):
            A = random_tournament(n, rng)
            for p in range(n):
                ap = enumerate_allowed_paths(A, n, p)
                apm1 = enumerate_allowed_paths(A, n, p-1) if p > 0 else []
                od = compute_omega_dim(A, n, p, ap, apm1)
                omega_sums[p] += od
                Ap_sums[p] += len(ap)
            count += 1

        print(f"\n  n={n} (d={n-1}), {count} samples:")
        for p in range(n):
            c = comb(n, p+1)
            f_omega = (omega_sums[p] / count) / c
            f_Ap = (Ap_sums[p] / count) / c
            # Expected |A_p| for random tournament: C(n,p+1) * (p+1)!/2^p
            # because each subset of size p+1 has (p+1)!/2^p expected Ham paths
            # (actually (p+1)!/2^{binom(p+1,2)} * 2^{binom(p+1,2)} = (p+1)!...no)
            # For random tournament: E[#Ham paths on k vertices] = k!/2^{C(k,2)-1}... no
            # Actually: each of the C(k,2) edge orientations is equally likely
            # Expected Ham paths = k!/2^{C(k,2)} * 2^{C(k,2)} ... = k! but divided by...
            # For random tournament on k vertices: E[H] = k!/2^{k-1}
            # So E[|A_p|] = C(n,p+1) * (p+1)!/2^p
            expected_mult = factorial(p+1) / (2**p)
            print(f"    p={p}: f_Ap={f_Ap:.4f} (expect {expected_mult:.2f}), "
                  f"f_Omega={f_omega:.4f}, "
                  f"Omega/Ap={omega_sums[p]/Ap_sums[p]:.4f}" if Ap_sums[p] > 0 else "")

    # Part 2: Paths per vertex set — what determines the "twist"?
    print("\n" + "=" * 70)
    print("Part 2: Paths per vertex set at each dimension")
    print("Transitive: exactly 1 path per set. Cyclic: > 1.")
    print("=" * 70)

    for n in [5, 7]:
        # One transitive, one Paley, one random
        configs = []
        # Transitive
        A_trans = np.array([[1 if j > i else 0 for j in range(n)] for i in range(n)])
        configs.append(("Transitive", A_trans))
        # Paley (if prime)
        if n in [3, 7]:
            A_paley = np.zeros((n, n), dtype=int)
            qr = {k*k % n for k in range(1, n)}
            for i in range(n):
                for j in range(n):
                    if i != j and ((j - i) % n) in qr:
                        A_paley[i][j] = 1
            configs.append(("Paley", A_paley))
        # Random
        rng = np.random.RandomState(42)
        A_rand = random_tournament(n, rng)
        configs.append(("Random", A_rand))

        for label, A in configs:
            print(f"\n  {label} T_{n}:")
            for p in range(n):
                paths = enumerate_allowed_paths(A, n, p)
                # Group by vertex set
                by_vset = defaultdict(int)
                for path in paths:
                    by_vset[frozenset(path)] += 1
                n_vsets = len(by_vset)
                c = comb(n, p+1)
                ppvs = sorted(by_vset.values())
                ppvs_dist = defaultdict(int)
                for v in ppvs:
                    ppvs_dist[v] += 1

                # Compute omega dim
                apm1 = enumerate_allowed_paths(A, n, p-1) if p > 0 else []
                od = compute_omega_dim(A, n, p, paths, apm1)
                filling = od / c if c > 0 else 0

                print(f"    p={p}: |A_p|={len(paths)}, vsets={n_vsets}/{c}, "
                      f"paths/vset={dict(sorted(ppvs_dist.items()))}, "
                      f"Omega={od}, f={filling:.3f}")

    # Part 3: Omega_p / C(n,p+1) as function of p/(n-1) — universal curve?
    print("\n" + "=" * 70)
    print("Part 3: Universal filling curve f(t) where t = p/(n-1)?")
    print("=" * 70)

    for n in range(4, 9):
        rng = np.random.RandomState(42 + n)
        N = min(100, 2**(n*(n-1)//2))

        omega_avg = defaultdict(float)
        cnt = 0
        for _ in range(N):
            A = random_tournament(n, rng)
            for p in range(n):
                ap = enumerate_allowed_paths(A, n, p)
                apm1 = enumerate_allowed_paths(A, n, p-1) if p > 0 else []
                omega_avg[p] += compute_omega_dim(A, n, p, ap, apm1)
            cnt += 1

        print(f"\n  n={n}: t = p/(n-1), f = <Omega_p> / C(n,p+1)")
        for p in range(n):
            c = comb(n, p+1)
            f = (omega_avg[p] / cnt) / c
            t = p / (n - 1)
            print(f"    t={t:.3f} (p={p}): f={f:.4f}")

if __name__ == '__main__':
    main()
