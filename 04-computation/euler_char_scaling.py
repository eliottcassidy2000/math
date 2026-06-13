"""
Euler characteristic scaling: chi(T) = sum (-1)^p beta_p across dimensions.

Known: chi(transitive) = 1 always. chi(random) is mostly 0 or 1.
But at n=8+, chi can be larger.

Key question: Is chi(T) determined by simpler invariants?

For the path complex:
chi = sum (-1)^p dim(Omega_p) - sum (-1)^p rank(d_p)
    = sum (-1)^p dim(Omega_p)  [since rank terms cancel in Euler char]

So chi depends ONLY on the Omega_p dimensions, not on the boundary maps!

By the rank-nullity identity for chain complexes:
chi = sum (-1)^p dim(Omega_p)

This is a PURELY COMBINATORIAL quantity computable from |A_p| and
the constraint matrices, without computing homology.

Let's verify this and see how chi relates to tournament structure.
"""
import numpy as np
from math import comb
from collections import defaultdict
from itertools import combinations

def random_tournament(n, rng):
    A = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

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

def is_strongly_connected(A, n):
    visited = set()
    stack = [0]
    while stack:
        v = stack.pop()
        if v in visited: continue
        visited.add(v)
        for u in range(n):
            if A[v][u] and u not in visited:
                stack.append(u)
    if len(visited) != n: return False
    visited = set()
    stack = [0]
    while stack:
        v = stack.pop()
        if v in visited: continue
        visited.add(v)
        for u in range(n):
            if A[u][v] and u not in visited:
                stack.append(u)
    return len(visited) == n

def main():
    print("=" * 70)
    print("EULER CHARACTERISTIC SCALING")
    print("chi = sum (-1)^p dim(Omega_p)")
    print("=" * 70)

    # Part 1: chi via Omega dims (no homology needed!)
    print("\n--- Part 1: chi = sum (-1)^p dim(Omega_p) verification ---")

    for n in range(3, 9):
        rng = np.random.RandomState(42 + n)
        N = min(300, 2**(n*(n-1)//2))

        chi_dist = defaultdict(int)
        chi_by_c3 = defaultdict(list)
        chi_by_sc = defaultdict(list)
        H_by_chi = defaultdict(list)
        cnt = 0

        for trial in range(N):
            if n <= 5:
                A = bits_to_adj(trial, n)
            else:
                A = random_tournament(n, rng)

            # Compute all Omega dims
            allowed = {}
            for p in range(-1, n):
                allowed[p] = [] if p < 0 else enumerate_allowed_paths(A, n, p)

            omega_dims = []
            for p in range(n):
                od = compute_omega_dim(A, n, p, allowed[p], allowed[p-1])
                omega_dims.append(od)

            chi = sum((-1)**p * omega_dims[p] for p in range(n))
            chi_dist[chi] += 1

            # Track c3
            c3 = 0
            for a, b, c in combinations(range(n), 3):
                if A[a][b] and A[b][c] and A[c][a]: c3 += 1
                if A[a][c] and A[c][b] and A[b][a]: c3 += 1
            chi_by_c3[c3].append(chi)

            sc = is_strongly_connected(A, n)
            chi_by_sc[sc].append(chi)

            # Compute H for small n
            if n <= 7:
                H = len(enumerate_allowed_paths(A, n, n-1))
                H_by_chi[chi].append(H)

            cnt += 1

        print(f"\n  n={n} (d={n-1}), {cnt} samples:")
        print(f"    chi distribution:")
        for chi_val in sorted(chi_dist.keys()):
            print(f"      chi={chi_val}: {chi_dist[chi_val]} ({100*chi_dist[chi_val]/cnt:.1f}%)")

        # chi vs c3
        print(f"    chi vs c3:")
        for c3_val in sorted(chi_by_c3.keys()):
            vals = chi_by_c3[c3_val]
            chi_vals_set = sorted(set(vals))
            if len(vals) > 3:
                avg_chi = np.mean(vals)
                print(f"      c3={c3_val}: avg_chi={avg_chi:.3f}, chi in {chi_vals_set}")

        # chi by SC
        for sc_val in [True, False]:
            if chi_by_sc[sc_val]:
                vals = chi_by_sc[sc_val]
                print(f"    SC={sc_val}: avg_chi={np.mean(vals):.3f}, range [{min(vals)},{max(vals)}]")

        # H by chi
        if H_by_chi:
            print(f"    H statistics by chi:")
            for chi_val in sorted(H_by_chi.keys()):
                vals = H_by_chi[chi_val]
                if vals:
                    print(f"      chi={chi_val}: avg H={np.mean(vals):.1f}, range [{min(vals)},{max(vals)}]")

    # Part 2: Alternating sum of |A_p| vs chi
    # chi_A = sum (-1)^p |A_p| — is this related to chi?
    print("\n" + "=" * 70)
    print("Part 2: Alternating sum of |A_p| vs chi(Omega)")
    print("chi_A = sum (-1)^p |A_p|")
    print("=" * 70)

    for n in range(3, 9):
        rng = np.random.RandomState(42 + n)
        N = min(200, 2**(n*(n-1)//2))

        chi_A_dist = defaultdict(int)
        chi_pairs = []  # (chi_Omega, chi_A)
        cnt = 0

        for trial in range(N):
            if n <= 5:
                A = bits_to_adj(trial, n)
            else:
                A = random_tournament(n, rng)

            allowed = {}
            for p in range(-1, n):
                allowed[p] = [] if p < 0 else enumerate_allowed_paths(A, n, p)

            omega_dims = []
            Ap_sizes = []
            for p in range(n):
                ap = allowed[p]
                Ap_sizes.append(len(ap))
                od = compute_omega_dim(A, n, p, ap, allowed[p-1])
                omega_dims.append(od)

            chi_Omega = sum((-1)**p * omega_dims[p] for p in range(n))
            chi_A = sum((-1)**p * Ap_sizes[p] for p in range(n))
            chi_A_dist[chi_A] += 1
            chi_pairs.append((chi_Omega, chi_A))
            cnt += 1

        print(f"\n  n={n}: chi_A distribution:")
        for val in sorted(chi_A_dist.keys()):
            print(f"    chi_A={val}: {chi_A_dist[val]} ({100*chi_A_dist[val]/cnt:.1f}%)")

        # Check if chi_A is always the same
        chi_A_vals = [p[1] for p in chi_pairs]
        if len(set(chi_A_vals)) == 1:
            print(f"    chi_A = {chi_A_vals[0]} ALWAYS! (tournament-independent)")
        else:
            print(f"    chi_A range: [{min(chi_A_vals)}, {max(chi_A_vals)}]")

        # Check if chi_Omega = chi_A always
        agree = sum(1 for p in chi_pairs if p[0] == p[1])
        print(f"    chi_Omega = chi_A: {agree}/{cnt}")

    # Part 3: chi(T) formula attempt
    # For the simplicial complex: chi(Delta_{n-1}) = 1 (contractible)
    # chi(simplex) = sum (-1)^p C(n,p+1) = ... = 1
    # Verify: sum (-1)^p C(n,p+1) for p=0..n-1
    print("\n" + "=" * 70)
    print("Part 3: Simplicial Euler characteristic")
    print("chi(Delta_{n-1}) = sum (-1)^p C(n,p+1) = 1")
    print("=" * 70)

    for n in range(2, 10):
        chi_simplex = sum((-1)**p * comb(n, p+1) for p in range(n))
        print(f"  n={n}: chi(Delta_{n-1}) = sum (-1)^p C({n},p+1) = {chi_simplex}")

    # Part 4: The DEFECT chi(T) - chi(simplex) = chi(T) - 1
    # When is chi(T) != 1?
    print("\n" + "=" * 70)
    print("Part 4: chi(T) - 1 distribution (topological defect)")
    print("chi(T) = 1: topologically simplex-like (contractible)")
    print("chi(T) = 0: homologically non-trivial")
    print("=" * 70)

    for n in range(3, 9):
        rng = np.random.RandomState(42 + n)
        N = min(500, 2**(n*(n-1)//2))

        defect_dist = defaultdict(int)
        cnt = 0

        for trial in range(N):
            if n <= 5:
                A = bits_to_adj(trial, n)
            else:
                A = random_tournament(n, rng)

            allowed = {}
            for p in range(-1, n):
                allowed[p] = [] if p < 0 else enumerate_allowed_paths(A, n, p)

            omega_dims = []
            for p in range(n):
                od = compute_omega_dim(A, n, p, allowed[p], allowed[p-1])
                omega_dims.append(od)

            chi = sum((-1)**p * omega_dims[p] for p in range(n))
            defect = chi - 1
            defect_dist[defect] += 1
            cnt += 1

        print(f"\n  n={n}: defect = chi - 1:")
        for d in sorted(defect_dist.keys()):
            pct = 100 * defect_dist[d] / cnt
            label = ""
            if d == 0: label = " (simplex-like)"
            elif d == -1: label = " (beta_1=1 or beta_3=1)"
            print(f"    defect={d:+d}: {defect_dist[d]} ({pct:.1f}%){label}")

if __name__ == '__main__':
    main()
