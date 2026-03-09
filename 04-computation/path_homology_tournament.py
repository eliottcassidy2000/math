"""
Compute GLMY path homology of tournaments.

Definitions (Grigoryan-Lin-Muranov-Yau):
- An n-path is a sequence (v0, v1, ..., vn) of vertices (consecutive distinct)
- The path is "allowed" if each (v_i, v_{i+1}) is an arc
- Lambda_n = free module on allowed n-paths
- boundary: d(v0,...,vn) = sum_{j=0}^{n} (-1)^j (v0,...,hat{v_j},...,vn)
  where hat means omit v_j
- Omega_n = {omega in Lambda_n : d(omega) in Lambda_{n-1}}
  i.e., all face (n-1)-tuples appearing in d(omega) must be allowed paths
- H_n = ker(d: Omega_n -> Omega_{n-1}) / im(d: Omega_{n+1} -> Omega_n)
"""

import numpy as np
from numpy.linalg import matrix_rank

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
    # Remove zero entries
    result = {k: v for k, v in result.items() if v != 0}
    return result

def compute_path_homology(adj, max_dim=4, verbose=True):
    """Compute GLMY path homology of a digraph given by adjacency matrix."""
    n = len(adj)
    results = {}

    # Compute allowed paths for each dimension
    allowed = {}
    for dim in range(max_dim + 2):
        allowed[dim] = get_allowed_paths(n, dim, adj)
        if verbose:
            print(f"  dim {dim}: {len(allowed[dim])} allowed {dim}-paths")

    # Index allowed paths
    path_index = {}
    allowed_set = {}
    for dim in range(max_dim + 2):
        path_index[dim] = {p: i for i, p in enumerate(allowed[dim])}
        allowed_set[dim] = set(allowed[dim])

    # Compute Omega_n for each dimension
    omega_bases = {}
    omega_bases[0] = np.eye(len(allowed[0])) if len(allowed[0]) > 0 else np.zeros((0, 0))

    for dim in range(1, max_dim + 2):
        num_paths = len(allowed[dim])
        if num_paths == 0:
            omega_bases[dim] = np.zeros((num_paths, 0))
            continue

        # For each allowed dim-path, compute boundary
        # Find non-allowed faces
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

            # Omega_dim = kernel of M
            U, S, Vt = np.linalg.svd(M, full_matrices=True)
            tol = 1e-10
            rank = np.sum(S > tol)
            if rank < num_paths:
                kernel = Vt[rank:].T
                omega_bases[dim] = kernel
            else:
                omega_bases[dim] = np.zeros((num_paths, 0))

        if verbose:
            odim = omega_bases[dim].shape[1] if omega_bases[dim].ndim == 2 else 0
            print(f"  Omega_{dim}: dimension {odim}")

    # Compute H_n = ker(d: Omega_n -> Omega_{n-1}) / im(d: Omega_{n+1} -> Omega_n)
    for dim in range(max_dim + 1):
        odim = omega_bases[dim].shape[1] if (omega_bases[dim].ndim == 2 and omega_bases[dim].shape[1] > 0) else 0
        if odim == 0:
            results[dim] = 0
            if verbose:
                print(f"  H_{dim} = 0 (Omega_{dim} = 0)")
            continue

        if dim == 0:
            # H_0: ker(d_0) / im(d_1|Omega_1)
            # d_0 = 0, so ker = all of Omega_0
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
            if verbose:
                print(f"  H_{dim}: beta_{dim} = {beta}")
            continue

        # Build boundary d_dim: Lambda_dim -> Lambda_{dim-1}
        D = np.zeros((len(allowed[dim-1]), len(allowed[dim])))
        for j, p in enumerate(allowed[dim]):
            bd = boundary_coefficients(p)
            for face, coeff in bd.items():
                if face in path_index[dim-1]:
                    D[path_index[dim-1][face], j] += coeff

        # Restrict to Omega_dim
        D_omega = D @ omega_bases[dim]

        # Compute kernel rank
        omega_prev_dim = omega_bases[dim-1].shape[1] if (omega_bases[dim-1].ndim == 2 and omega_bases[dim-1].shape[1] > 0) else 0
        if omega_prev_dim == 0:
            ker_rank = odim
        else:
            D_in_omega = np.linalg.lstsq(omega_bases[dim-1], D_omega, rcond=None)[0]
            ker_rank = odim - matrix_rank(D_in_omega, tol=1e-10)

        # Image of d_{dim+1}
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
        if verbose:
            print(f"  H_{dim}: ker={ker_rank}, im={im_rank}, beta_{dim} = {beta}")

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


if __name__ == "__main__":
    print("=" * 60)
    print("PATH HOMOLOGY OF SMALL TOURNAMENTS")
    print("=" * 60)

    # n=3: Transitive
    print("\n--- n=3: Transitive tournament ---")
    adj_trans3 = [[False, True, True], [False, False, True], [False, False, False]]
    print("Adjacency: 0->1, 0->2, 1->2")
    msf, info = check_multisquare_free(adj_trans3)
    print(f"Multisquare-free: {msf}")
    h = compute_path_homology(adj_trans3, max_dim=3)
    print(f"Betti numbers: {h}")

    # n=3: Cyclic
    print("\n--- n=3: Cyclic tournament ---")
    adj_cyc3 = [[False, True, False], [False, False, True], [True, False, False]]
    print("Adjacency: 0->1, 1->2, 2->0")
    msf, info = check_multisquare_free(adj_cyc3)
    print(f"Multisquare-free: {msf}")
    h = compute_path_homology(adj_cyc3, max_dim=3)
    print(f"Betti numbers: {h}")

    # n=4: All tournaments
    print("\n--- n=4: All tournaments ---")
    num_edges_4 = 4*3//2
    betti_counts = {}
    msf_count = 0
    for bits in range(2**num_edges_4):
        adj = make_tournament(4, bits)
        msf, info = check_multisquare_free(adj)
        if msf:
            msf_count += 1
        h = compute_path_homology(adj, max_dim=3, verbose=False)
        key = tuple(h.get(i, 0) for i in range(4))
        if key not in betti_counts:
            betti_counts[key] = {"count": 0, "msf": 0, "example_bits": bits}
        betti_counts[key]["count"] += 1
        if msf:
            betti_counts[key]["msf"] += 1

    print(f"Total 4-vertex tournaments: {2**num_edges_4}")
    print(f"Multisquare-free: {msf_count}")
    print(f"\nBetti number distributions:")
    for key in sorted(betti_counts.keys()):
        info = betti_counts[key]
        print(f"  (beta_0,...,beta_3) = {key}: {info['count']} tournaments ({info['msf']} MSF)")

    # n=5: All tournaments
    print("\n--- n=5: All tournaments ---")
    num_edges_5 = 5*4//2
    betti_counts_5 = {}
    msf_count_5 = 0
    for bits in range(2**num_edges_5):
        adj = make_tournament(5, bits)
        msf, info = check_multisquare_free(adj)
        if msf:
            msf_count_5 += 1
        h = compute_path_homology(adj, max_dim=3, verbose=False)
        key = tuple(h.get(i, 0) for i in range(4))
        if key not in betti_counts_5:
            betti_counts_5[key] = {"count": 0, "msf": 0}
        betti_counts_5[key]["count"] += 1
        if msf:
            betti_counts_5[key]["msf"] += 1

    print(f"Total 5-vertex tournaments: {2**num_edges_5}")
    print(f"Multisquare-free: {msf_count_5}")
    print(f"\nBetti number distributions:")
    for key in sorted(betti_counts_5.keys()):
        info = betti_counts_5[key]
        print(f"  (beta_0,...,beta_3) = {key}: {info['count']} tournaments ({info['msf']} MSF)")

    all_b2_zero_4 = all(key[2] == 0 for key in betti_counts.keys())
    all_b2_zero_5 = all(key[2] == 0 for key in betti_counts_5.keys())
    print(f"\nAll n=4 tournaments have beta_2 = 0: {all_b2_zero_4}")
    print(f"All n=5 tournaments have beta_2 = 0: {all_b2_zero_5}")
