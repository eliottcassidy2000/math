"""
Compute GLMY path homology of tournaments at n=6 and sample n=7.
Check if beta_2 = 0 universally.
"""

import numpy as np
from numpy.linalg import matrix_rank
import random
import time

def get_allowed_paths(n_vertices, length, adj):
    """Get all allowed paths of given length (number of edges)."""
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
    """Compute d(path) as a dict: (n-1)-tuple -> coefficient."""
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
    """Compute GLMY path homology of a digraph given by adjacency matrix."""
    n = len(adj)
    results = {}

    allowed = {}
    for dim in range(max_dim + 2):
        allowed[dim] = get_allowed_paths(n, dim, adj)

    path_index = {}
    allowed_set = {}
    for dim in range(max_dim + 2):
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

            U, S, Vt = np.linalg.svd(M, full_matrices=True)
            tol = 1e-10
            rank = np.sum(S > tol)
            if rank < num_paths:
                kernel = Vt[rank:].T
                omega_bases[dim] = kernel
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

def make_tournament(n, bits):
    """Create tournament adjacency matrix from bit encoding."""
    adj = [[False]*n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << k):
                adj[i][j] = True
            else:
                adj[j][i] = True
            k += 1
    return adj

def check_multisquare_free(adj):
    """Check if digraph is multisquare-free."""
    n = len(adj)
    for x in range(n):
        for y in range(n):
            if x == y:
                continue
            count = 0
            for z in range(n):
                if z != x and z != y and adj[x][z] and adj[z][y]:
                    count += 1
            if count >= 3:
                return False, (x, y, count)
    return True, None

def random_tournament(n):
    """Generate a random tournament on n vertices."""
    adj = [[False]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                adj[i][j] = True
            else:
                adj[j][i] = True
    return adj


if __name__ == "__main__":
    # n=6: exhaustive (2^15 = 32768 tournaments)
    print("=" * 60)
    print("n=6: Exhaustive enumeration (32768 tournaments)")
    print("=" * 60)

    num_edges_6 = 6*5//2  # = 15
    betti_counts_6 = {}
    msf_count_6 = 0
    t0 = time.time()

    for bits in range(2**num_edges_6):
        adj = make_tournament(6, bits)
        msf, info = check_multisquare_free(adj)
        if msf:
            msf_count_6 += 1
        h = compute_path_homology(adj, max_dim=3, verbose=False)
        key = tuple(h.get(i, 0) for i in range(4))
        if key not in betti_counts_6:
            betti_counts_6[key] = {"count": 0, "msf": 0}
        betti_counts_6[key]["count"] += 1
        if msf:
            betti_counts_6[key]["msf"] += 1

        if (bits + 1) % 5000 == 0:
            elapsed = time.time() - t0
            print(f"  Progress: {bits+1}/{2**num_edges_6} ({elapsed:.1f}s)")

    elapsed = time.time() - t0
    print(f"\nCompleted in {elapsed:.1f}s")
    print(f"Total 6-vertex tournaments: {2**num_edges_6}")
    print(f"Multisquare-free: {msf_count_6}")
    print(f"\nBetti number distributions:")
    for key in sorted(betti_counts_6.keys()):
        info = betti_counts_6[key]
        print(f"  (beta_0,...,beta_3) = {key}: {info['count']} tournaments ({info['msf']} MSF)")

    all_b2_zero = all(key[2] == 0 for key in betti_counts_6.keys())
    print(f"\nAll n=6 tournaments have beta_2 = 0: {all_b2_zero}")

    # n=7: Sample
    print("\n" + "=" * 60)
    print("n=7: Random sample (1000 tournaments)")
    print("=" * 60)

    random.seed(42)
    betti_counts_7 = {}
    msf_count_7 = 0
    n_samples = 1000
    t0 = time.time()

    for i in range(n_samples):
        adj = random_tournament(7)
        msf, info = check_multisquare_free(adj)
        if msf:
            msf_count_7 += 1
        h = compute_path_homology(adj, max_dim=3, verbose=False)
        key = tuple(h.get(j, 0) for j in range(4))
        if key not in betti_counts_7:
            betti_counts_7[key] = {"count": 0, "msf": 0}
        betti_counts_7[key]["count"] += 1
        if msf:
            betti_counts_7[key]["msf"] += 1

        if (i + 1) % 200 == 0:
            elapsed = time.time() - t0
            print(f"  Progress: {i+1}/{n_samples} ({elapsed:.1f}s)")

    elapsed = time.time() - t0
    print(f"\nCompleted in {elapsed:.1f}s")
    print(f"Sampled {n_samples} random 7-vertex tournaments")
    print(f"Multisquare-free: {msf_count_7}")
    print(f"\nBetti number distributions:")
    for key in sorted(betti_counts_7.keys()):
        info = betti_counts_7[key]
        print(f"  (beta_0,...,beta_3) = {key}: {info['count']} tournaments ({info['msf']} MSF)")

    all_b2_zero_7 = all(key[2] == 0 for key in betti_counts_7.keys())
    print(f"\nAll sampled n=7 tournaments have beta_2 = 0: {all_b2_zero_7}")
