"""
Check even Betti numbers (beta_2, beta_4, beta_6) for tournaments.

Uses the correct GLMY path homology computation where Omega_p is the
subspace of allowed p-paths whose boundary has no non-allowed face terms.

Key question: Is beta_4 always 0? OPEN-Q-024 claims yes, but THM-099
says Paley T_7 has beta_4 = 6. This script resolves the contradiction.
"""
import numpy as np
from itertools import permutations, combinations
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
    """All directed paths of length p with distinct vertices."""
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
    """Alternating face map: remove each vertex with sign."""
    return [((-1)**i, path[:i] + path[i+1:]) for i in range(len(path))]

def build_boundary_matrix(allowed_p, allowed_pm1):
    """Full boundary matrix (only terms landing in allowed_{p-1})."""
    if not allowed_p or not allowed_pm1:
        return np.zeros((max(len(allowed_pm1), 0), max(len(allowed_p), 0)))
    idx = {path: i for i, path in enumerate(allowed_pm1)}
    M = np.zeros((len(allowed_pm1), len(allowed_p)))
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if face in idx:
                M[idx[face], j] += sign
    return M

def compute_omega_basis(A, n, p, allowed_p, allowed_pm1):
    """
    Omega_p = subspace of allowed p-paths whose boundary has NO non-allowed faces.
    Returns basis as columns of a matrix (rows = allowed_p dimension).
    """
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

def path_betti_numbers(A, n, max_dim=None):
    """Compute Betti numbers beta_0, ..., beta_{max_dim}."""
    if max_dim is None: max_dim = n - 1
    allowed = {}
    for p in range(-1, max_dim + 2):
        if p < 0:
            allowed[p] = []
        else:
            allowed[p] = enumerate_allowed_paths(A, n, p)
    omega = {}
    for p in range(max_dim + 2):
        omega[p] = compute_omega_basis(A, n, p, allowed[p], allowed[p-1])
    betti = []
    for p in range(max_dim + 1):
        dim_om = omega[p].shape[1] if omega[p].ndim == 2 else 0
        if dim_om == 0:
            betti.append(0)
            continue
        # ker(d_p restricted to Omega_p)
        bd_p = build_boundary_matrix(allowed[p], allowed[p-1])
        bd_p_om = bd_p @ omega[p]
        if bd_p_om.size > 0:
            sv = np.linalg.svd(bd_p_om, compute_uv=False)
            rk = sum(s > 1e-8 for s in sv)
        else:
            rk = 0
        ker = dim_om - rk
        # im(d_{p+1} restricted to Omega_{p+1})
        dim_om1 = omega[p+1].shape[1] if omega[p+1].ndim == 2 else 0
        if dim_om1 > 0:
            bd1 = build_boundary_matrix(allowed[p+1], allowed[p])
            bd1_om = bd1 @ omega[p+1]
            sv1 = np.linalg.svd(bd1_om, compute_uv=False)
            im = sum(s > 1e-8 for s in sv1)
        else:
            im = 0
        betti.append(ker - im)
    return betti

def main():
    rng = np.random.RandomState(42)

    # Test 1: n=3 3-cycle
    print("=" * 70)
    print("TEST 1: 3-cycle tournament (a->b->c->a)")
    print("=" * 70)
    n = 3
    A = np.array([[0,1,0],[0,0,1],[1,0,0]])
    betti = path_betti_numbers(A, n)
    print(f"3-cycle: betti = {betti}")

    # Transitive tournament a->b, a->c, b->c
    A2 = np.array([[0,1,1],[0,0,1],[0,0,0]])
    betti2 = path_betti_numbers(A2, n)
    print(f"Transitive: betti = {betti2}")

    # Test 2: n=5 exhaustive
    print("\n" + "=" * 70)
    print("TEST 2: n=5 exhaustive")
    print("=" * 70)
    n = 5
    total = 2**(n*(n-1)//2)
    betti_dist = defaultdict(int)
    for bits in range(total):
        A = bits_to_adj(bits, n)
        betti = path_betti_numbers(A, n)
        betti_dist[tuple(betti)] += 1

    print(f"n={n}: {total} tournaments")
    for bv, count in sorted(betti_dist.items(), key=lambda x: -x[1]):
        print(f"  {list(bv)}: {count} ({100*count/total:.1f}%)")

    # Test 3: Paley T_7
    print("\n" + "=" * 70)
    print("TEST 3: Paley T_7")
    print("=" * 70)
    n = 7
    A = np.zeros((n, n), dtype=int)
    qr = {1, 2, 4}
    for i in range(n):
        for j in range(n):
            if i != j and ((j - i) % n) in qr:
                A[i][j] = 1
    betti = path_betti_numbers(A, n)
    print(f"Paley T_7: betti = {betti}")

    # Also test dimension by dimension
    for p in range(n):
        allowed_p = enumerate_allowed_paths(A, n, p)
        omega_p = compute_omega_basis(A, n, p, allowed_p,
                                      enumerate_allowed_paths(A, n, p-1) if p > 0 else [])
        dim_A = len(allowed_p)
        dim_Om = omega_p.shape[1] if omega_p.ndim == 2 and omega_p.shape[1] > 0 else 0
        print(f"  p={p}: |A_p|={dim_A}, dim(Omega_p)={dim_Om}, beta_p={betti[p]}")

    # Test 4: n=7 sampling for even Betti numbers
    print("\n" + "=" * 70)
    print("TEST 4: n=7 random sampling — even Betti check")
    print("=" * 70)
    n = 7
    nsamp = 300  # Small sample due to n=7 being slow
    beta2_vals = defaultdict(int)
    beta4_vals = defaultdict(int)
    beta6_vals = defaultdict(int)

    for idx in range(nsamp):
        A = random_tournament(n, rng)
        betti = path_betti_numbers(A, n)
        beta2_vals[betti[2]] += 1
        beta4_vals[betti[4]] += 1
        if len(betti) > 6:
            beta6_vals[betti[6]] += 1

        if idx % 50 == 0:
            print(f"  {idx}/{nsamp}...")

    print(f"\nn=7 ({nsamp} random):")
    print(f"  beta_2 distribution: {dict(sorted(beta2_vals.items()))}")
    print(f"  beta_4 distribution: {dict(sorted(beta4_vals.items()))}")
    print(f"  beta_6 distribution: {dict(sorted(beta6_vals.items()))}")

    beta2_nonzero = sum(v for k, v in beta2_vals.items() if k != 0)
    beta4_nonzero = sum(v for k, v in beta4_vals.items() if k != 0)
    beta6_nonzero = sum(v for k, v in beta6_vals.items() if k != 0)
    print(f"\n  beta_2 != 0: {beta2_nonzero}/{nsamp}")
    print(f"  beta_4 != 0: {beta4_nonzero}/{nsamp}")
    print(f"  beta_6 != 0: {beta6_nonzero}/{nsamp}")

if __name__ == '__main__':
    main()
