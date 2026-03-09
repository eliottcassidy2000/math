"""
Check beta_5 for MULTIPLE n=8 H-maximizers (H=661).

The first sample had beta_5=0. At n=6, only half of H-maximizers had
beta_3>0 (the Paley deletion half). So we need multiple H=661 tournaments.

Score class: (3,3,3,3,4,4,4,4) — the SC class at n=8.
"""
import numpy as np
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

def ham_path_count(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        bits = bin(mask).count('1')
        if bits == 1: continue
        for v in range(n):
            if not (mask & (1 << v)): continue
            prev_mask = mask ^ (1 << v)
            for u in range(n):
                if (prev_mask & (1 << u)) and A[u][v]:
                    if (mask, v) not in dp:
                        dp[(mask, v)] = 0
                    dp[(mask, v)] += dp.get((prev_mask, u), 0)
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

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

def betti_single(A, n, target_p):
    """Compute just beta_{target_p} efficiently."""
    allowed = {}
    for p in [target_p - 1, target_p, target_p + 1]:
        if p < 0:
            allowed[p] = []
        else:
            allowed[p] = enumerate_allowed_paths(A, n, p)

    # Need omega at target_p and target_p+1
    omega_p = compute_omega_basis(A, n, target_p, allowed[target_p], allowed[target_p-1])
    dim_om = omega_p.shape[1] if omega_p.ndim == 2 else 0
    if dim_om == 0:
        return 0, dim_om

    # Kernel of d_p
    bd_p = build_boundary_matrix(allowed[target_p], allowed[target_p-1])
    bd_p_om = bd_p @ omega_p
    if bd_p_om.size > 0:
        sv = np.linalg.svd(bd_p_om, compute_uv=False)
        rk = sum(s > 1e-8 for s in sv)
    else:
        rk = 0
    ker = dim_om - rk

    # Need omega at target_p+1 for image
    omega_p1 = compute_omega_basis(A, n, target_p+1, allowed[target_p+1], allowed[target_p])
    dim_om1 = omega_p1.shape[1] if omega_p1.ndim == 2 else 0
    if dim_om1 > 0:
        bd1 = build_boundary_matrix(allowed[target_p+1], allowed[target_p])
        bd1_om = bd1 @ omega_p1
        sv1 = np.linalg.svd(bd1_om, compute_uv=False)
        im = sum(s > 1e-8 for s in sv1)
    else:
        im = 0

    return ker - im, dim_om

def main():
    rng = np.random.RandomState(123)
    n = 8
    target_score = (3,3,3,3,4,4,4,4)

    print("=" * 70)
    print(f"Finding multiple H=661 tournaments at n={n}")
    print(f"Checking beta_5 (= beta_{{n-3}}) for each")
    print("=" * 70)

    h661_found = []
    h661_beta5 = []
    h661_odim5 = []

    trials = 0
    while len(h661_found) < 20 and trials < 200000:
        A = random_tournament(n, rng)
        trials += 1
        scores = tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))
        if scores != target_score:
            continue
        H = ham_path_count(A, n)
        if H != 661:
            continue

        # Found an H=661 tournament! Compute beta_5
        b5, odim5 = betti_single(A, n, 5)
        h661_found.append(A.copy())
        h661_beta5.append(b5)
        h661_odim5.append(odim5)
        print(f"  #{len(h661_found)}: beta_5={b5}, dim(Omega_5)={odim5} (trial {trials})")

    print(f"\nTotal trials: {trials}")
    print(f"H=661 found: {len(h661_found)}")
    print(f"beta_5 values: {h661_beta5}")
    print(f"dim(Omega_5) values: {h661_odim5}")

    b5_dist = defaultdict(int)
    for b in h661_beta5:
        b5_dist[b] += 1
    print(f"beta_5 distribution: {dict(sorted(b5_dist.items()))}")

    odim5_dist = defaultdict(int)
    for o in h661_odim5:
        odim5_dist[o] += 1
    print(f"Omega_5 dim distribution: {dict(sorted(odim5_dist.items()))}")

    # Also check other Betti for the first few
    print("\n" + "=" * 70)
    print("Full Betti for first 3 H=661 tournaments")
    print("=" * 70)

    for idx in range(min(3, len(h661_found))):
        A = h661_found[idx]
        betti_full = []
        odims_full = []
        for p in range(n):
            bp, op = betti_single(A, n, p)
            betti_full.append(bp)
            odims_full.append(op)
        print(f"  #{idx+1}: betti = {betti_full}")
        print(f"        omega = {odims_full}")
        chi = sum((-1)**p * betti_full[p] for p in range(n))
        print(f"        chi = {chi}")

    # Check: what H values do SC-score tournaments have?
    print("\n" + "=" * 70)
    print(f"H distribution for score {target_score}")
    print("=" * 70)

    h_dist = defaultdict(int)
    rng2 = np.random.RandomState(456)
    for trial in range(50000):
        A = random_tournament(n, rng2)
        scores = tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))
        if scores != target_score:
            continue
        H = ham_path_count(A, n)
        h_dist[H] += 1

    print(f"H values found (from {sum(h_dist.values())} SC-score tournaments):")
    for H in sorted(h_dist.keys()):
        marker = " <-- MAX" if H == 661 else ""
        print(f"  H={H}: {h_dist[H]}{marker}")

if __name__ == '__main__':
    main()
