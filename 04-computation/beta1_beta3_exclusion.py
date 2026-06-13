"""
WHY are beta_1 and beta_3 mutually exclusive?

Key observations:
- beta_1 > 0 means the tournament has a "1-cycle" in homology (kernel of d_1 not in image of d_2)
- beta_3 > 0 means there's a "3-cycle" in homology
- They never cooccur (confirmed n=6 exhaustive, n=7 sampled)

Structural investigation:
1. What tournaments have beta_1 > 0? Characterize them.
2. What tournaments have beta_3 > 0? Characterize them.
3. What structural property precludes both?

At n=6 (exhaustive):
- beta_1=1: 4800/32768 (14.6%) — related to connectivity?
- beta_3=1: 320/32768 (1.0%) — related to H-maximization?

beta_1 counts 1-dimensional "holes": closed allowed paths (in Omega_1)
that are not boundaries of 2-chains (in Omega_2).

Investigate: Is beta_1 related to strong connectivity properties?
Is beta_3 related to having many directed cycles?
"""
import numpy as np
from itertools import combinations
from collections import defaultdict

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

def is_strongly_connected(A, n):
    """Check if tournament is strongly connected."""
    # BFS from vertex 0 using forward edges
    visited = set()
    stack = [0]
    while stack:
        v = stack.pop()
        if v in visited: continue
        visited.add(v)
        for u in range(n):
            if A[v][u] and u not in visited:
                stack.append(u)
    if len(visited) != n:
        return False
    # BFS from vertex 0 using reverse edges
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

def strong_components(A, n):
    """Compute number of strong components via Tarjan-like approach."""
    # Use Kosaraju's algorithm
    # Pass 1: DFS on original graph, record finish order
    visited = set()
    finish = []
    def dfs1(v):
        visited.add(v)
        for u in range(n):
            if A[v][u] and u not in visited:
                dfs1(u)
        finish.append(v)
    for v in range(n):
        if v not in visited:
            dfs1(v)
    # Pass 2: DFS on reverse graph in reverse finish order
    visited2 = set()
    num_components = 0
    def dfs2(v):
        visited2.add(v)
        for u in range(n):
            if A[u][v] and u not in visited2:
                dfs2(u)
    for v in reversed(finish):
        if v not in visited2:
            dfs2(v)
            num_components += 1
    return num_components

def vertex_connectivity(A, n):
    """Approximate vertex connectivity: min number of vertices whose removal
    disconnects the strongly connected tournament."""
    if not is_strongly_connected(A, n):
        return 0
    # Check kappa >= 1: remove each vertex
    for v in range(n):
        keep = [i for i in range(n) if i != v]
        B = np.zeros((n-1, n-1), dtype=int)
        for i, ki in enumerate(keep):
            for j, kj in enumerate(keep):
                B[i][j] = A[ki][kj]
        if not is_strongly_connected(B, n-1):
            return 1
    return 2  # kappa >= 2

def count_3cycles(A, n):
    count = 0
    for a, b, c in combinations(range(n), 3):
        if A[a][b] and A[b][c] and A[c][a]: count += 1
        if A[a][c] and A[c][b] and A[b][a]: count += 1
    return count

def main():
    n = 6
    total = 2**(n*(n-1)//2)

    print("=" * 70)
    print(f"Structural characterization of beta_1>0 and beta_3>0 at n={n}")
    print("=" * 70)

    b1_data = []
    b3_data = []
    trivial_data = []

    for bits in range(total):
        A = bits_to_adj(bits, n)
        betti, odims = path_betti_numbers(A, n)

        sc = is_strongly_connected(A, n)
        kappa = vertex_connectivity(A, n) if sc else 0
        c3 = count_3cycles(A, n)
        scores = tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))

        entry = {
            'bits': bits,
            'sc': sc,
            'kappa': kappa,
            'c3': c3,
            'scores': scores,
            'betti': list(betti),
            'odims': list(odims),
        }

        if betti[1] > 0:
            b1_data.append(entry)
        elif len(betti) > 3 and betti[3] > 0:
            b3_data.append(entry)
        else:
            trivial_data.append(entry)

    print(f"\nTotal: {total}")
    print(f"beta_1>0: {len(b1_data)} ({100*len(b1_data)/total:.1f}%)")
    print(f"beta_3>0: {len(b3_data)} ({100*len(b3_data)/total:.1f}%)")
    print(f"trivial: {len(trivial_data)} ({100*len(trivial_data)/total:.1f}%)")

    # Strong connectivity
    print("\n--- Strong connectivity ---")
    b1_sc = sum(1 for e in b1_data if e['sc'])
    b3_sc = sum(1 for e in b3_data if e['sc'])
    triv_sc = sum(1 for e in trivial_data if e['sc'])
    print(f"beta_1>0 SC: {b1_sc}/{len(b1_data)} ({100*b1_sc/len(b1_data):.1f}%)")
    print(f"beta_3>0 SC: {b3_sc}/{len(b3_data)} ({100*b3_sc/len(b3_data):.1f}%)")
    print(f"trivial SC: {triv_sc}/{len(trivial_data)} ({100*triv_sc/len(trivial_data):.1f}%)")

    # Vertex connectivity
    print("\n--- Vertex connectivity (kappa) ---")
    for name, data in [("beta_1>0", b1_data), ("beta_3>0", b3_data), ("trivial", trivial_data)]:
        kappa_dist = defaultdict(int)
        for e in data:
            kappa_dist[e['kappa']] += 1
        print(f"{name}: kappa dist = {dict(sorted(kappa_dist.items()))}")

    # 3-cycle counts
    print("\n--- 3-cycle counts ---")
    for name, data in [("beta_1>0", b1_data), ("beta_3>0", b3_data), ("trivial", trivial_data)]:
        c3_vals = [e['c3'] for e in data]
        print(f"{name}: c3 range [{min(c3_vals)}, {max(c3_vals)}], mean={np.mean(c3_vals):.1f}")

    # Score distributions
    print("\n--- Score distributions ---")
    for name, data in [("beta_1>0", b1_data), ("beta_3>0", b3_data)]:
        score_dist = defaultdict(int)
        for e in data:
            score_dist[e['scores']] += 1
        print(f"\n{name} score classes:")
        for s in sorted(score_dist.keys()):
            print(f"  {s}: {score_dist[s]}")

    # Omega dims comparison
    print("\n--- Omega dims ---")
    for name, data in [("beta_1>0", b1_data), ("beta_3>0", b3_data)]:
        odim_dist = defaultdict(int)
        for e in data:
            odim_dist[tuple(e['odims'])] += 1
        print(f"\n{name} Omega dim patterns (top 5):")
        for od, cnt in sorted(odim_dist.items(), key=lambda x: -x[1])[:5]:
            print(f"  {list(od)}: {cnt}")

    # Key question: what's the COMPLEMENTARY structural difference?
    print("\n" + "=" * 70)
    print("Key structural comparison: beta_1 vs beta_3")
    print("=" * 70)

    # Can we find a single invariant that separates them?
    # Try: number of non-allowed 1-paths (edges) — these are "missing edges"
    # A non-allowed 1-path (a,b) means a->b exists but a is not reachable from b
    # (i.e., there's no edge b->a, so b doesn't have both an in-edge from a and out-edge to a)
    # Wait, in a tournament a->b means A[a][b]=1, and since it's a tournament, A[b][a]=0.
    # So the 1-path (a,b) is allowed iff A[a][b]=1 — ALL directed edges are allowed 1-paths!
    # There are NO non-allowed 1-paths in a tournament.
    # Therefore dim(Omega_1) = dim(A_1) = number of edges with the correct direction.

    # What about 2-paths? (a,b,c) allowed iff A[a][b]=1 and A[b][c]=1 and a,b,c distinct.
    # A 2-path has boundary: (b,c) - (a,c) + (a,b).
    # All faces are edges of the tournament (or reversed), so all are allowed UNLESS
    # the face reverses an edge. Face (b,c): edge b->c? Yes iff A[b][c]=1.
    # If A[b][c]=0 (i.e., c->b), then face (b,c) is NOT allowed.

    # So non-allowed faces of 2-paths come from "non-transitive" triples.

    # Check: for beta_1>0, count how many non-allowed (p-1)-faces each p-path has
    for name, data in [("beta_1>0", b1_data[:5]), ("beta_3>0", b3_data[:5])]:
        print(f"\n{name} examples:")
        for e in data[:3]:
            A = bits_to_adj(e['bits'], n)
            # Omega dims
            print(f"  bits={e['bits']}: scores={e['scores']}, c3={e['c3']}, kappa={e['kappa']}")
            print(f"    betti={e['betti']}, omega={e['odims']}")

            # Count non-transitive triples
            nt = 0
            for a, b, c in combinations(range(n), 3):
                if A[a][b] and A[b][c] and A[c][a]: nt += 1
                if A[a][c] and A[c][b] and A[b][a]: nt += 1
            print(f"    directed 3-cycles: {nt}")

if __name__ == '__main__':
    main()
