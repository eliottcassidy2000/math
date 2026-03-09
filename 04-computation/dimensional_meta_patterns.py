"""
Meta-patterns across dimensions: tournament = binary relation on (n-1)-simplex.

Dimension d = n-1. Track how homological invariants scale with d.

Key quantities to track per dimension:
1. Fraction with beta_1 > 0 (1-hole rate)
2. Fraction with beta_3 > 0 (3-hole rate)
3. Fraction with ALL beta_p = 0 for p >= 1 (homologically trivial)
4. Average dim(Omega_p) / C(n, p+1) — "filling fraction"
5. Euler characteristic distribution
6. Maximum beta values

Framing: Omega_p lives in the space of allowed (p+1)-vertex paths.
The simplex has C(n, p+1) subsets of size p+1.
Each subset can support between 0 and (p+1)! directed paths.
In a tournament, each subset supports EXACTLY some number of directed
Hamiltonian paths on those vertices.

For a TRANSITIVE tournament: |A_p| = C(n, p+1) exactly (one path per subset).
For CYCLIC subtournaments: |A_p| > C(n, p+1) (multiple paths per subset).
"""
import numpy as np
from itertools import combinations
from collections import defaultdict
import sys

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

def compute_omega_basis(A, n, p, allowed_p, allowed_pm1):
    dim_Ap = len(allowed_p)
    if dim_Ap == 0: return np.zeros((0, 0))
    if p == 0: return np.eye(dim_Ap)
    allowed_pm1_set = set(allowed_pm1)
    non_allowed = {}
    na_count = 0
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if len(set(face)) == len(face) and face not in allowed_pm1_set:
                if face not in non_allowed:
                    non_allowed[face] = na_count
                    na_count += 1
    if na_count == 0: return np.eye(dim_Ap)
    P = np.zeros((na_count, dim_Ap))
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if face in non_allowed:
                P[non_allowed[face], j] += sign
    U, S, Vt = np.linalg.svd(P, full_matrices=True)
    rank = sum(s > 1e-10 for s in S)
    ns = Vt[rank:].T
    return ns if ns.shape[1] > 0 else np.zeros((dim_Ap, 0))

def build_boundary_matrix(allowed_p, allowed_pm1):
    if not allowed_p or not allowed_pm1:
        return np.zeros((max(len(allowed_pm1), 0), max(len(allowed_p), 0)))
    idx = {path: i for i, path in enumerate(allowed_pm1)}
    M = np.zeros((len(allowed_pm1), len(allowed_p)))
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if face in idx:
                M[idx[face], j] += sign
    return M

def path_betti_numbers(A, n, max_dim=None):
    if max_dim is None: max_dim = n - 1
    allowed = {}
    for p in range(-1, max_dim + 2):
        allowed[p] = [] if p < 0 else enumerate_allowed_paths(A, n, p)
    omega = {}
    omega_dims = []
    for p in range(max_dim + 2):
        omega[p] = compute_omega_basis(A, n, p, allowed[p], allowed[p-1])
        omega_dims.append(omega[p].shape[1] if omega[p].ndim == 2 else 0)
    betti = []
    for p in range(max_dim + 1):
        dim_om = omega_dims[p]
        if dim_om == 0:
            betti.append(0)
            continue
        bd_p = build_boundary_matrix(allowed[p], allowed[p-1])
        bd_p_om = bd_p @ omega[p]
        if bd_p_om.size > 0:
            sv = np.linalg.svd(bd_p_om, compute_uv=False)
            rk = sum(s > 1e-8 for s in sv)
        else:
            rk = 0
        ker = dim_om - rk
        dim_om1 = omega_dims[p+1]
        if dim_om1 > 0:
            bd1 = build_boundary_matrix(allowed[p+1], allowed[p])
            bd1_om = bd1 @ omega[p+1]
            sv1 = np.linalg.svd(bd1_om, compute_uv=False)
            im = sum(s > 1e-8 for s in sv1)
        else:
            im = 0
        betti.append(ker - im)
    return betti, omega_dims[:max_dim+1]

def fast_omega_dims_and_Ap(A, n, max_dim=None):
    """Compute just |A_p| and dim(Omega_p), no Betti."""
    if max_dim is None: max_dim = n - 1
    allowed = {}
    for p in range(-1, max_dim + 1):
        allowed[p] = [] if p < 0 else enumerate_allowed_paths(A, n, p)
    Ap_sizes = []
    omega_dims = []
    for p in range(max_dim + 1):
        Ap_sizes.append(len(allowed[p]))
        od = compute_omega_dim(A, n, p, allowed[p], allowed[p-1])
        omega_dims.append(od)
    return Ap_sizes, omega_dims

def main():
    print("=" * 70)
    print("DIMENSIONAL META-PATTERNS: Tournament = Binary Relation on Simplex")
    print("Dimension d = n-1. Tracking homological invariants vs d.")
    print("=" * 70)

    # Part 1: Exhaustive at small n, sampling at larger n
    for n in range(3, 9):
        d = n - 1
        num_edges = n * (n - 1) // 2
        total = 2 ** num_edges

        print(f"\n{'='*60}")
        print(f"d = {d} (n = {n}): {num_edges} edges, 2^{num_edges} = {total} tournaments")
        print(f"{'='*60}")

        if n <= 6:
            # Exhaustive
            sample_size = total
            tournaments = range(total)
            is_exhaustive = True
        else:
            # Sample
            sample_size = min(500, total)
            rng = np.random.RandomState(42 + n)
            is_exhaustive = False

        betti_dist = defaultdict(int)
        chi_dist = defaultdict(int)
        Ap_totals = defaultdict(float)
        omega_totals = defaultdict(float)
        max_betti = defaultdict(int)
        count = 0

        for trial_idx in range(sample_size):
            if is_exhaustive:
                A = bits_to_adj(trial_idx, n)
            else:
                A = random_tournament(n, rng)

            # Compute Betti
            max_dim = min(n - 1, 7)  # cap at 7 for speed
            betti, odims = path_betti_numbers(A, n, max_dim=max_dim)

            bv = tuple(int(b) for b in betti)
            betti_dist[bv] += 1

            chi = sum((-1)**p * betti[p] for p in range(len(betti)))
            chi_dist[int(chi)] += 1

            # Track max Betti per dimension
            for p, b in enumerate(betti):
                max_betti[p] = max(max_betti[p], int(b))

            # Track omega dims and |A_p|
            Ap_list = []
            for p in range(max_dim + 1):
                Ap_list.append(len(enumerate_allowed_paths(A, n, p)) if p <= max_dim else 0)
            for p in range(min(len(Ap_list), len(odims))):
                Ap_totals[p] += Ap_list[p]
                omega_totals[p] += odims[p]

            count += 1

        # Report
        mode = "exhaustive" if is_exhaustive else f"sampled ({sample_size})"
        print(f"  Mode: {mode}")

        # Betti vector distribution (top 5)
        print(f"\n  Betti vector distribution (top 5):")
        for bv, cnt in sorted(betti_dist.items(), key=lambda x: -x[1])[:5]:
            pct = 100 * cnt / count
            chi = sum((-1)**p * bv[p] for p in range(len(bv)))
            print(f"    {list(bv)}: {cnt} ({pct:.1f}%), chi={chi}")

        # Chi distribution
        print(f"\n  Euler characteristic distribution:")
        for chi_val in sorted(chi_dist.keys()):
            cnt = chi_dist[chi_val]
            print(f"    chi={chi_val}: {cnt} ({100*cnt/count:.1f}%)")

        # Max Betti per dimension
        print(f"\n  Max beta_p across all tournaments:")
        for p in sorted(max_betti.keys()):
            if max_betti[p] > 0 or p <= 4:
                print(f"    max beta_{p} = {max_betti[p]}")

        # Fraction with nonzero beta_p
        for p in range(1, min(6, len(list(betti_dist.keys())[0]))):
            n_pos = sum(cnt for bv, cnt in betti_dist.items() if len(bv) > p and bv[p] > 0)
            if n_pos > 0 or p <= 3:
                print(f"    P(beta_{p} > 0) = {n_pos}/{count} ({100*n_pos/count:.2f}%)")

        # Average filling fractions
        print(f"\n  Average dim(Omega_p) / C(n,p+1):")
        from math import comb
        for p in range(min(6, max(len(bv) for bv in betti_dist.keys()))):
            c = comb(n, p + 1)
            if c > 0:
                avg_omega = omega_totals[p] / count
                avg_Ap = Ap_totals[p] / count
                filling = avg_omega / c if c > 0 else 0
                path_mult = avg_Ap / c if c > 0 else 0
                print(f"    p={p}: avg|A_p|={avg_Ap:.1f} ({path_mult:.2f}*C), "
                      f"avg Omega_p={avg_omega:.1f} ({filling:.3f}*C)")

    # Part 2: Summary table — rates across dimensions
    print("\n" + "=" * 70)
    print("SUMMARY: Homological rates across dimensions")
    print("=" * 70)
    print("d  | P(b1>0)  | P(b3>0)  | P(trivial) | chi range  | max b1 | max b3")
    print("---|----------|----------|------------|------------|--------|-------")

    # Re-run quick computation for summary
    for n in range(3, 9):
        d = n - 1
        num_edges = n * (n - 1) // 2
        total = 2 ** num_edges

        if n <= 6:
            sample_size = total
            is_exhaustive = True
        else:
            sample_size = 500
            is_exhaustive = False

        rng = np.random.RandomState(42 + n)
        b1_pos = 0
        b3_pos = 0
        trivial = 0
        chi_vals = set()
        max_b1 = 0
        max_b3 = 0
        cnt = 0

        for trial_idx in range(sample_size):
            if is_exhaustive:
                A = bits_to_adj(trial_idx, n)
            else:
                A = random_tournament(n, rng)

            betti, _ = path_betti_numbers(A, n, max_dim=min(n-1, 5))
            chi = sum((-1)**p * betti[p] for p in range(len(betti)))
            chi_vals.add(int(chi))

            b1 = int(betti[1]) if len(betti) > 1 else 0
            b3 = int(betti[3]) if len(betti) > 3 else 0
            if b1 > 0: b1_pos += 1
            if b3 > 0: b3_pos += 1
            if all(b == 0 for b in betti[1:]): trivial += 1
            max_b1 = max(max_b1, b1)
            max_b3 = max(max_b3, b3)
            cnt += 1

        chi_range = sorted(chi_vals)
        print(f"{d:2d} | {100*b1_pos/cnt:7.2f}% | {100*b3_pos/cnt:7.2f}% | "
              f"{100*trivial/cnt:9.2f}% | {chi_range} | {max_b1:6d} | {max_b3:5d}")

if __name__ == '__main__':
    main()
