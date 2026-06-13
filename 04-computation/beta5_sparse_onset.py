"""
Beta_5 onset search using sparse SVD methods.

Previous attempt: Paley T_11 hit OOM with dense 53790x53790 matrix.
Solution: Use scipy.sparse for the constraint and boundary matrices.

Targets:
1. Paley T_7: known beta_5=0 (dense computation confirmed)
2. Paley T_11: needs sparse methods (n=11, p=5 has ~53k allowed 6-paths)
3. Random n=10 tournaments
4. n=9 specifically (previous search found 0 beta_5>0 in 200 random)

Strategy for sparse:
- Build constraint matrix P as scipy.sparse.csc_matrix
- Use scipy.sparse.linalg.svds for rank computation
- For boundary matrices, same approach
"""

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import svds
from math import comb
import time

def random_tournament(n, rng):
    A = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def paley_tournament(p):
    n = p
    A = np.zeros((n, n), dtype=int)
    qr = {k*k % n for k in range(1, n)}
    for i in range(n):
        for j in range(n):
            if i != j and ((j-i) % n) in qr:
                A[i][j] = 1
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

def compute_omega_dim_sparse(n, p, allowed_p, allowed_pm1):
    """Compute dim(Omega_p) using sparse matrices."""
    dim_Ap = len(allowed_p)
    if dim_Ap == 0: return 0
    if p == 0: return dim_Ap

    allowed_pm1_set = set(allowed_pm1)
    non_allowed = {}
    na_count = 0

    # First pass: identify non-allowed faces
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if len(set(face)) == len(face) and face not in allowed_pm1_set:
                if face not in non_allowed:
                    non_allowed[face] = na_count
                    na_count += 1

    if na_count == 0: return dim_Ap

    # Build sparse constraint matrix
    rows, cols, vals = [], [], []
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if face in non_allowed:
                rows.append(non_allowed[face])
                cols.append(j)
                vals.append(float(sign))

    P = sparse.csr_matrix((vals, (rows, cols)), shape=(na_count, dim_Ap))

    # Compute rank using SVD
    # For large matrices, use iterative SVD
    min_dim = min(na_count, dim_Ap)
    if min_dim == 0:
        return dim_Ap

    if min_dim <= 500:
        # Dense SVD for small matrices
        P_dense = P.toarray()
        sv = np.linalg.svd(P_dense, compute_uv=False)
        rank = sum(s > 1e-10 for s in sv)
    else:
        # Iterative approach: binary search for rank using svds
        # Start by trying to find rank via QR decomposition
        try:
            from scipy.sparse.linalg import LinearOperator
            # Use QR for rank
            P_dense = P.toarray()  # might still be feasible if shape is moderate
            sv = np.linalg.svd(P_dense, compute_uv=False)
            rank = sum(s > 1e-10 for s in sv)
        except MemoryError:
            # Fall back to iterative SVD
            # Estimate rank by probing
            k_try = min(min_dim - 1, 100)
            sv = svds(P.astype(float), k=k_try, return_singular_vectors=False)
            rank_lower = sum(s > 1e-10 for s in sv)
            if rank_lower < k_try:
                rank = rank_lower
            else:
                # Need more singular values
                k_try = min(min_dim - 1, 500)
                sv = svds(P.astype(float), k=k_try, return_singular_vectors=False)
                rank = sum(s > 1e-10 for s in sv)
                print(f"    WARNING: may underestimate rank (tried {k_try} sv, got rank {rank})")

    return dim_Ap - rank

def betti_sparse(A, n, target_p, verbose=False):
    """Compute single Betti number using sparse methods."""
    t0 = time.time()

    allowed = {}
    for p in [target_p - 1, target_p, target_p + 1]:
        if p < 0:
            allowed[p] = []
        elif p < n:
            t1 = time.time()
            allowed[p] = enumerate_allowed_paths(A, n, p)
            if verbose:
                print(f"    |A_{p}| = {len(allowed[p])} ({time.time()-t1:.1f}s)")
        else:
            allowed[p] = []

    # dim(Omega_p)
    t1 = time.time()
    dim_om = compute_omega_dim_sparse(n, target_p, allowed[target_p], allowed[target_p-1])
    if verbose:
        print(f"    dim(Omega_{target_p}) = {dim_om} ({time.time()-t1:.1f}s)")

    if dim_om == 0:
        return 0

    # Build boundary map d_p: Omega_p -> Omega_{p-1} (restricted to allowed)
    # Actually need: ker(d_p on Omega_p) - im(d_{p+1} on Omega_{p+1})
    # Use the same approach but with the Omega basis

    # For sparse computation, compute Betti via:
    # beta_p = dim(Omega_p) - rank(d_p|Omega_p) - rank(d_{p+1}|Omega_{p+1}) where these compose
    # Actually: beta_p = ker(d_p|Omega_p) - im(d_{p+1}|Omega_{p+1})

    # Compute rank(d_p restricted to Omega_p)
    # This is the boundary map from Omega_p to Omega_{p-1}
    t1 = time.time()
    allowed_pm1_idx = {path: i for i, path in enumerate(allowed[target_p-1])}
    allowed_p_idx = {path: i for i, path in enumerate(allowed[target_p])}

    # Build boundary matrix (in allowed basis)
    rows_bd, cols_bd, vals_bd = [], [], []
    for j, path in enumerate(allowed[target_p]):
        for sign, face in boundary_coeffs(path):
            if face in allowed_pm1_idx:
                rows_bd.append(allowed_pm1_idx[face])
                cols_bd.append(j)
                vals_bd.append(float(sign))

    if rows_bd:
        bd_p = sparse.csr_matrix((vals_bd, (rows_bd, cols_bd)),
                                  shape=(len(allowed[target_p-1]), len(allowed[target_p])))
    else:
        ker_om = dim_om
        bd_p = None

    # Need Omega_p basis first to restrict
    # For large n, computing the full Omega basis is expensive
    # Alternative: compute Betti directly via double rank computation

    # beta_p = dim(Omega_p) - rank(d_p on allowed_p | constraint=0) - im(d_{p+1})
    # This is getting complex. Use the direct formula:
    # Omega_p = ker(constraint_p), d_p|Omega_p = d_p|ker(constraint_p)
    # beta_p = dim(ker(d_p) intersect Omega_p) - dim(im(d_{p+1}) intersect Omega_p)

    # Simpler approach: just compute omega dims at p-1, p, p+1 and boundary ranks
    # within the Omega subspaces.

    # For now, fall back to the dense approach for dimensions that fit in memory
    dim_Ap = len(allowed[target_p])
    if dim_Ap > 20000:
        print(f"    dim(A_{target_p}) = {dim_Ap}, too large for dense Omega basis")
        print(f"    Using rank-based estimation...")

        # Compute ranks of augmented system
        # [constraint_p; boundary_p] to get rank of boundary restricted to Omega
        # This is still large...

        # Alternative: just report dim(Omega_p) and flag as incomplete
        print(f"    dim(Omega_{target_p}) = {dim_om}")
        print(f"    Need specialized sparse Betti computation for full answer")
        return None

    # Dense approach for moderate sizes
    from numpy.linalg import svd as np_svd

    # Compute Omega basis
    allowed_pm1_set = set(allowed[target_p-1])
    non_allowed = {}
    na_count = 0
    for j, path in enumerate(allowed[target_p]):
        for sign, face in boundary_coeffs(path):
            if len(set(face)) == len(face) and face not in allowed_pm1_set:
                if face not in non_allowed:
                    non_allowed[face] = na_count
                    na_count += 1

    if na_count == 0:
        omega_basis = np.eye(dim_Ap)
    else:
        P = np.zeros((na_count, dim_Ap))
        for j, path in enumerate(allowed[target_p]):
            for sign, face in boundary_coeffs(path):
                if face in non_allowed:
                    P[non_allowed[face], j] += sign
        U, S, Vt = np_svd(P, full_matrices=True)
        rank = sum(s > 1e-10 for s in S)
        omega_basis = Vt[rank:].T
        if omega_basis.shape[1] == 0:
            return 0

    dim_omega = omega_basis.shape[1]
    if verbose:
        print(f"    Omega basis shape: {omega_basis.shape}")

    # Boundary d_p restricted to Omega_p
    if bd_p is not None:
        bd_om = bd_p.toarray() @ omega_basis
    else:
        bd_om = np.zeros((len(allowed[target_p-1]), dim_omega))

    sv = np_svd(bd_om, compute_uv=False)
    rk_dp = sum(s > 1e-8 for s in sv)
    ker_dp = dim_omega - rk_dp

    # im(d_{p+1} on Omega_{p+1})
    dim_Ap1 = len(allowed[target_p+1])
    if dim_Ap1 == 0:
        im_dp1 = 0
    elif dim_Ap1 > 20000:
        print(f"    dim(A_{target_p+1}) = {dim_Ap1}, too large for dense")
        return None
    else:
        # Compute Omega_{p+1} basis
        allowed_p_set = set(allowed[target_p])
        non_allowed_p1 = {}
        na_count_p1 = 0
        for j, path in enumerate(allowed[target_p+1]):
            for sign, face in boundary_coeffs(path):
                if len(set(face)) == len(face) and face not in allowed_p_set:
                    if face not in non_allowed_p1:
                        non_allowed_p1[face] = na_count_p1
                        na_count_p1 += 1

        if na_count_p1 == 0:
            omega_p1_basis = np.eye(dim_Ap1)
        else:
            P1 = np.zeros((na_count_p1, dim_Ap1))
            for j, path in enumerate(allowed[target_p+1]):
                for sign, face in boundary_coeffs(path):
                    if face in non_allowed_p1:
                        P1[non_allowed_p1[face], j] += sign
            U1, S1, Vt1 = np_svd(P1, full_matrices=True)
            rank1 = sum(s > 1e-10 for s in S1)
            omega_p1_basis = Vt1[rank1:].T
            if omega_p1_basis.shape[1] == 0:
                im_dp1 = 0
                return ker_dp - im_dp1

        dim_omega_p1 = omega_p1_basis.shape[1]

        # Boundary d_{p+1}: A_{p+1} -> A_p
        rows_bd1, cols_bd1, vals_bd1 = [], [], []
        for j, path in enumerate(allowed[target_p+1]):
            for sign, face in boundary_coeffs(path):
                if face in allowed_p_idx:
                    rows_bd1.append(allowed_p_idx[face])
                    cols_bd1.append(j)
                    vals_bd1.append(float(sign))

        if rows_bd1:
            bd_p1 = sparse.csr_matrix((vals_bd1, (rows_bd1, cols_bd1)),
                                       shape=(dim_Ap, dim_Ap1))
            bd_p1_om = bd_p1.toarray() @ omega_p1_basis
        else:
            bd_p1_om = np.zeros((dim_Ap, dim_omega_p1))

        sv1 = np_svd(bd_p1_om, compute_uv=False)
        im_dp1 = sum(s > 1e-8 for s in sv1)

    betti = ker_dp - im_dp1
    if verbose:
        print(f"    ker(d_{target_p}) = {ker_dp}, im(d_{target_p+1}) = {im_dp1}")
        print(f"    beta_{target_p} = {betti}")
        print(f"    Total time: {time.time()-t0:.1f}s")

    return betti

def main():
    print("=" * 70)
    print("BETA_5 ONSET SEARCH (SPARSE METHODS)")
    print("=" * 70)

    # Part 1: Verify known results
    print("\n--- Part 1: Verify known results ---")

    # Paley T_7: beta_5 should be 0
    A7 = paley_tournament(7)
    print("\nPaley T_7:")
    for p in [1, 3, 5]:
        b = betti_sparse(A7, 7, p, verbose=True)
        print(f"  beta_{p} = {b}")

    # Part 2: Random search at n=9 (expanded from 200 to 1000)
    print("\n--- Part 2: Random n=9, beta_5 ---")
    n = 9
    rng = np.random.RandomState(12345)
    N = 500
    b5_pos = 0
    cnt = 0

    for trial in range(N):
        A = random_tournament(n, rng)
        b5 = betti_sparse(A, n, 5)
        if b5 is not None and b5 > 0:
            b5_pos += 1
            print(f"  FOUND beta_5={b5} at trial {trial}!")
        cnt += 1
        if cnt % 100 == 0:
            print(f"  {cnt}/{N}: beta_5>0 = {b5_pos}/{cnt}")

    print(f"\n  n=9: beta_5>0 = {b5_pos}/{cnt}")

    # Part 3: Try n=10
    print("\n--- Part 3: Random n=10, beta_5 ---")
    n = 10
    rng = np.random.RandomState(54321)
    N = 100
    b5_pos = 0
    cnt = 0

    for trial in range(N):
        A = random_tournament(n, rng)
        # At n=10, |A_5| could be very large
        # First check: how many allowed 6-paths?
        allowed_5 = enumerate_allowed_paths(A, n, 5)
        if len(allowed_5) > 15000:
            if trial < 5:
                print(f"  Trial {trial}: |A_5| = {len(allowed_5)}, skipping (too large)")
            continue

        b5 = betti_sparse(A, n, 5, verbose=(trial < 3))
        if b5 is not None and b5 > 0:
            b5_pos += 1
            print(f"  FOUND beta_5={b5} at trial {trial}!")
        cnt += 1
        if cnt > 0 and cnt % 20 == 0:
            print(f"  {cnt} computed: beta_5>0 = {b5_pos}/{cnt}")

    print(f"\n  n=10: beta_5>0 = {b5_pos}/{cnt}")

    # Part 4: Paley T_11 with path count estimation
    print("\n--- Part 4: Paley T_11 path counts ---")
    A11 = paley_tournament(11)
    print("Paley T_11:")
    for p in [1, 3, 5]:
        t0 = time.time()
        allowed = enumerate_allowed_paths(A11, 11, p)
        t1 = time.time()
        print(f"  |A_{p}| = {len(allowed)} ({t1-t0:.1f}s)")

    # Can we compute dim(Omega_5) at least?
    print("\n  Computing dim(Omega_5) for Paley T_11...")
    for p in [4, 5, 6]:
        t0 = time.time()
        ap = enumerate_allowed_paths(A11, 11, p)
        apm1 = enumerate_allowed_paths(A11, 11, p-1)
        print(f"  |A_{p}| = {len(ap)}, |A_{p-1}| = {len(apm1)}", end="")
        if len(ap) <= 20000:
            dim_om = compute_omega_dim_sparse(11, p, ap, apm1)
            print(f" -> dim(Omega_{p}) = {dim_om} ({time.time()-t0:.1f}s)")
        else:
            print(f" -> too large for dense ({time.time()-t0:.1f}s)")

if __name__ == '__main__':
    main()
