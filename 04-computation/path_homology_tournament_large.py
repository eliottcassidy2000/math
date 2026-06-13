"""
Sample path homology of tournaments at n=7,8,9 to check beta_2=0 universality.
Also look at higher Betti numbers.
"""

import numpy as np
from numpy.linalg import matrix_rank
import random
import time

def get_allowed_paths(n_vertices, length, adj):
    paths = []
    vertices = list(range(n_vertices))
    def gen(current_path):
        if len(current_path) == length + 1:
            paths.append(tuple(current_path))
            return
        for v in vertices:
            if v != current_path[-1] and adj[current_path[-1]][v]:
                gen(current_path + [v])
    for v in vertices:
        gen([v])
    return paths

def boundary_coefficients(path):
    result = {}
    for j in range(len(path)):
        face = path[:j] + path[j+1:]
        sign = (-1)**j
        if face in result:
            result[face] += sign
        else:
            result[face] = sign
    result = {k: v for k, v in result.items() if v != 0}
    return result

def compute_path_homology(adj, max_dim=4, verbose=False):
    n = len(adj)
    results = {}

    allowed = {}
    for dim in range(max_dim + 2):
        allowed[dim] = get_allowed_paths(n, dim, adj)
        if len(allowed[dim]) == 0 and dim >= 2:
            # No paths beyond this dimension
            for d in range(dim, max_dim + 2):
                allowed[d] = []
            break

    path_index = {}
    allowed_set = {}
    for dim in range(max_dim + 2):
        if dim not in allowed:
            allowed[dim] = []
        path_index[dim] = {p: i for i, p in enumerate(allowed[dim])}
        allowed_set[dim] = set(allowed[dim])

    omega_bases = {}
    omega_bases[0] = np.eye(len(allowed[0])) if len(allowed[0]) > 0 else np.zeros((0, 0))

    for dim in range(1, max_dim + 2):
        num_paths = len(allowed[dim])
        if num_paths == 0:
            omega_bases[dim] = np.zeros((num_paths, 0))
            continue

        non_allowed_faces = set()
        boundary_data = []
        for p in allowed[dim]:
            bd = boundary_coefficients(p)
            boundary_data.append(bd)
            for face in bd:
                if face not in allowed_set[dim-1]:
                    non_allowed_faces.add(face)

        if not non_allowed_faces:
            omega_bases[dim] = np.eye(num_paths)
        else:
            na_list = sorted(non_allowed_faces)
            na_index = {f: i for i, f in enumerate(na_list)}
            M = np.zeros((len(na_list), num_paths))
            for j, bd in enumerate(boundary_data):
                for face, coeff in bd.items():
                    if face in na_index:
                        M[na_index[face], j] += coeff

            try:
                U, S, Vt = np.linalg.svd(M, full_matrices=True)
                tol = 1e-10
                rank = np.sum(S > tol)
                if rank < num_paths:
                    kernel = Vt[rank:].T
                    omega_bases[dim] = kernel
                else:
                    omega_bases[dim] = np.zeros((num_paths, 0))
            except np.linalg.LinAlgError:
                # Fallback: use QR on M^T for kernel
                Q, R = np.linalg.qr(M.T, mode='complete')
                tol = 1e-10
                rank = np.sum(np.abs(np.diag(R[:min(M.shape), :])) > tol) if min(M.shape) > 0 else 0
                if rank < num_paths:
                    omega_bases[dim] = Q[:, rank:]
                else:
                    omega_bases[dim] = np.zeros((num_paths, 0))

    for dim in range(max_dim + 1):
        odim = omega_bases[dim].shape[1] if (omega_bases[dim].ndim == 2 and omega_bases[dim].shape[1] > 0) else 0
        if odim == 0:
            results[dim] = 0
            continue

        if dim == 0:
            im_rank = 0
            odim1 = omega_bases[1].shape[1] if (omega_bases[1].ndim == 2 and omega_bases[1].shape[1] > 0) else 0
            if odim1 > 0:
                D1 = np.zeros((len(allowed[0]), len(allowed[1])))
                for j, p in enumerate(allowed[1]):
                    bd = boundary_coefficients(p)
                    for face, coeff in bd.items():
                        if face in path_index[0]:
                            D1[path_index[0][face], j] += coeff
                D1_omega = D1 @ omega_bases[1]
                im_rank = matrix_rank(D1_omega, tol=1e-10)
            beta = odim - im_rank
            results[dim] = beta
            continue

        D = np.zeros((len(allowed[dim-1]), len(allowed[dim])))
        for j, p in enumerate(allowed[dim]):
            bd = boundary_coefficients(p)
            for face, coeff in bd.items():
                if face in path_index[dim-1]:
                    D[path_index[dim-1][face], j] += coeff

        D_omega = D @ omega_bases[dim]

        omega_prev_dim = omega_bases[dim-1].shape[1] if (omega_bases[dim-1].ndim == 2 and omega_bases[dim-1].shape[1] > 0) else 0
        if omega_prev_dim == 0:
            ker_rank = odim
        else:
            D_in_omega = np.linalg.lstsq(omega_bases[dim-1], D_omega, rcond=None)[0]
            ker_rank = odim - matrix_rank(D_in_omega, tol=1e-10)

        im_rank = 0
        next_dim = dim + 1
        if next_dim <= max_dim + 1:
            odim_next = omega_bases[next_dim].shape[1] if (omega_bases[next_dim].ndim == 2 and omega_bases[next_dim].shape[1] > 0) else 0
            if odim_next > 0:
                D_next = np.zeros((len(allowed[dim]), len(allowed[next_dim])))
                for j, p in enumerate(allowed[next_dim]):
                    bd = boundary_coefficients(p)
                    for face, coeff in bd.items():
                        if face in path_index[dim]:
                            D_next[path_index[dim][face], j] += coeff
                D_next_omega = D_next @ omega_bases[next_dim]
                D_next_in_omega = np.linalg.lstsq(omega_bases[dim], D_next_omega, rcond=None)[0]
                im_rank = matrix_rank(D_next_in_omega, tol=1e-10)

        beta = ker_rank - im_rank
        results[dim] = beta

    return results

def random_tournament(n):
    adj = [[False]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                adj[i][j] = True
            else:
                adj[j][i] = True
    return adj


if __name__ == "__main__":
    random.seed(12345)

    for n in [7, 8, 9]:
        n_samples = 500 if n <= 8 else 100
        max_d = min(n-1, 4)  # Keep max_dim reasonable

        print("=" * 60)
        print(f"n={n}: Random sample ({n_samples} tournaments, max_dim={max_d})")
        print("=" * 60)

        betti_counts = {}
        t0 = time.time()

        for i in range(n_samples):
            adj = random_tournament(n)
            h = compute_path_homology(adj, max_dim=max_d, verbose=False)
            key = tuple(h.get(j, 0) for j in range(max_d + 1))
            if key not in betti_counts:
                betti_counts[key] = 0
            betti_counts[key] += 1

            if (i + 1) % 100 == 0:
                elapsed = time.time() - t0
                print(f"  Progress: {i+1}/{n_samples} ({elapsed:.1f}s)")

        elapsed = time.time() - t0
        print(f"\nCompleted in {elapsed:.1f}s")
        print(f"\nBetti number distributions:")
        for key in sorted(betti_counts.keys()):
            print(f"  betas = {key}: {betti_counts[key]} tournaments")

        # Check beta_2 = 0
        all_b2_zero = all(key[2] == 0 for key in betti_counts.keys())
        print(f"\nAll sampled n={n} tournaments have beta_2 = 0: {all_b2_zero}")

        # Check any higher betti nonzero
        max_nonzero = 0
        for key in betti_counts.keys():
            for j in range(len(key)):
                if key[j] > 0:
                    max_nonzero = max(max_nonzero, j)
        print(f"Highest dimension with nonzero beta: {max_nonzero}")
        print()
