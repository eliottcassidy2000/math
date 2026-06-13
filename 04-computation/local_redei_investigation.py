"""
Local Rédei investigation: Every subtournament has odd # of Ham paths.

Observation from filling_ratio_formula.py: The number of directed
Hamiltonian paths on any subset of vertices is ALWAYS ODD.

This is a direct consequence of Rédei's theorem applied to each
induced subtournament: every tournament has an odd number of
Hamiltonian paths.

But the DISTRIBUTION of these counts carries structural information:
- Transitive: always 1 (minimum)
- Cyclic 3-tournament: always 3
- Near-transitive: usually 1 or 3
- Highly cyclic: larger odd numbers

Questions:
1. What is the distribution of H(T[S]) for random subtournament S of size k?
2. How does this relate to Omega_p dimensions?
3. Is there a formula for dim(Omega_p) in terms of {H(T[S]) : |S|=p+1}?
4. The SURPLUS dim(Omega_p) - C(n,p+1) = function of cyclic content?
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

def ham_paths_on_subset(A, subset):
    """Count Hamiltonian paths on induced subtournament."""
    sub = sorted(subset)
    k = len(sub)
    if k <= 1: return 1
    # Build sub-adjacency
    sub_A = np.zeros((k, k), dtype=int)
    for i, u in enumerate(sub):
        for j, v in enumerate(sub):
            sub_A[i][j] = A[u][v]
    # Count via bitmask DP
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
    print("LOCAL RÉDEI: Hamiltonian path counts on subtournaments")
    print("Every subtournament has ODD # of Ham paths (by Rédei)")
    print("=" * 70)

    # Part 1: Distribution of H(T[S]) for subsets of size k
    print("\n--- Part 1: H(T[S]) distribution for random tournaments ---")

    for n in range(4, 9):
        rng = np.random.RandomState(42 + n)
        N = min(100, 2**(n*(n-1)//2))

        for k in range(3, min(n+1, 7)):
            h_dist = defaultdict(int)
            total_subsets = 0
            for _ in range(N):
                A = random_tournament(n, rng)
                for S in combinations(range(n), k):
                    h = ham_paths_on_subset(A, S)
                    h_dist[h] += 1
                    total_subsets += 1

            c = comb(n, k)
            print(f"\n  n={n}, k={k}: {total_subsets} subsets ({N} tours x C({n},{k})={c})")
            for h_val in sorted(h_dist.keys()):
                pct = 100 * h_dist[h_val] / total_subsets
                print(f"    H={h_val}: {h_dist[h_val]} ({pct:.1f}%)")

    # Part 2: Correlation between surplus Omega_p and sum of H(T[S])
    print("\n" + "=" * 70)
    print("Part 2: dim(Omega_p) vs sum_{|S|=p+1} H(T[S])")
    print("If Omega_p = C(n,p+1), all subsets have H=1 (transitive-like)")
    print("Surplus = dim(Omega_p) - C(n,p+1)")
    print("Total H = sum_{|S|=p+1} H(T[S])")
    print("Excess H = Total H - C(n,p+1) = sum of (H(T[S])-1)")
    print("=" * 70)

    for n in [5, 6, 7]:
        rng = np.random.RandomState(42 + n)
        N = min(200, 2**(n*(n-1)//2))

        for p in [2, 3]:
            k = p + 1
            if k > n: continue

            data = []  # (omega_surplus, h_excess) pairs
            for _ in range(N):
                A = random_tournament(n, rng)

                # Compute Omega_p
                ap = enumerate_allowed_paths(A, n, p)
                apm1 = enumerate_allowed_paths(A, n, p-1) if p > 0 else []
                od = compute_omega_dim(A, n, p, ap, apm1)
                surplus = od - comb(n, k)

                # Compute total H on all k-subsets
                total_h = sum(ham_paths_on_subset(A, S) for S in combinations(range(n), k))
                excess_h = total_h - comb(n, k)

                data.append((surplus, excess_h))

            surpluses = [d[0] for d in data]
            excesses = [d[1] for d in data]

            if len(set(surpluses)) > 1 and len(set(excesses)) > 1:
                corr = np.corrcoef(surpluses, excesses)[0, 1]
            else:
                corr = float('nan')

            print(f"\n  n={n}, p={p}: correlation(surplus, excess_H) = {corr:.4f}")
            print(f"    surplus range: [{min(surpluses)}, {max(surpluses)}]")
            print(f"    excess_H range: [{min(excesses)}, {max(excesses)}]")

            # Check if there's a formula
            # Is surplus = some function of excess_H?
            unique_pairs = defaultdict(list)
            for s, e in data:
                unique_pairs[e].append(s)

            if len(unique_pairs) <= 20:
                print(f"    (excess_H -> surplus) mapping:")
                for e in sorted(unique_pairs.keys()):
                    vals = unique_pairs[e]
                    if len(set(vals)) == 1:
                        print(f"      excess_H={e}: surplus={vals[0]} (deterministic)")
                    else:
                        print(f"      excess_H={e}: surplus in {sorted(set(vals))} ({len(vals)} samples)")

    # Part 3: Does dim(Omega_p) = |A_p| - rank(constraint matrix)?
    # The constraint matrix has rows = non-allowed (p-1)-paths
    # rank = number of independent constraints
    # dim(Omega_p) = |A_p| - rank
    # So surplus = |A_p| - rank - C(n,p+1)
    #           = (|A_p| - C(n,p+1)) - rank
    #           = excess_H_total - rank
    # Therefore: surplus = excess_paths - rank_constraints
    print("\n" + "=" * 70)
    print("Part 3: Omega_p = |A_p| - rank(constraints)")
    print("surplus = excess_paths - rank = (|A_p| - C) - rank")
    print("=" * 70)

    for n in [5, 6, 7]:
        rng = np.random.RandomState(42 + n)
        N = min(100, 2**(n*(n-1)//2))

        for p in [2, 3]:
            k = p + 1
            if k > n: continue

            decomp = []
            for _ in range(N):
                A = random_tournament(n, rng)
                ap = enumerate_allowed_paths(A, n, p)
                apm1 = enumerate_allowed_paths(A, n, p-1) if p > 0 else []
                od = compute_omega_dim(A, n, p, ap, apm1)

                Ap = len(ap)
                C = comb(n, k)
                rank = Ap - od  # constraint rank
                excess_paths = Ap - C
                surplus = od - C

                decomp.append((excess_paths, rank, surplus))

            ep = [d[0] for d in decomp]
            rk = [d[1] for d in decomp]
            su = [d[2] for d in decomp]

            print(f"\n  n={n}, p={p}:")
            print(f"    avg |A_p| - C = {np.mean(ep):.1f}, avg rank = {np.mean(rk):.1f}, avg surplus = {np.mean(su):.1f}")
            print(f"    Identity: surplus = excess_paths - rank: {all(d[2] == d[0] - d[1] for d in decomp)}")

            # What fraction of excess paths survive?
            survival = [d[2] / d[0] if d[0] > 0 else float('nan') for d in decomp]
            valid_survival = [s for s in survival if not np.isnan(s)]
            if valid_survival:
                print(f"    avg survival = surplus/excess = {np.mean(valid_survival):.4f}")

if __name__ == '__main__':
    main()
