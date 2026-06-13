"""
Deeper investigation of beta_4 > 0 ↔ H-maximizer connection.

Finding from beta4_characterization.py: At n=7, beta_4 = 6 iff H = 189 (Paley class).

Questions:
1. Do ALL Paley-class tournaments have beta_4 = 6? (Yes, 6/6 tested)
2. Is this related to the Omega palindromic dimensions?
3. At n=6, does beta_4 behavior exist?
4. What is the algebraic structure? Euler char, Betti sum, etc.
"""
import numpy as np
from itertools import permutations, combinations
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

def main():
    rng = np.random.RandomState(42)

    # Part 1: Exhaustive n=6 (check beta_4 at n=6)
    print("=" * 70)
    print("PART 1: n=6 exhaustive — all Betti numbers")
    print("=" * 70)
    n = 6
    total = 2**(n*(n-1)//2)
    betti_dist = defaultdict(int)
    for bits in range(total):
        A = bits_to_adj(bits, n)
        betti, _ = path_betti_numbers(A, n)
        betti_dist[tuple(betti)] += 1
    print(f"n={n}: {total} tournaments")
    for bv, count in sorted(betti_dist.items(), key=lambda x: -x[1]):
        print(f"  {list(bv)}: {count} ({100*count/total:.1f}%)")

    # Part 2: Euler characteristic at n=7
    print("\n" + "=" * 70)
    print("PART 2: n=7 Paley — Euler characteristic and Poincare polynomial")
    print("=" * 70)
    n = 7
    A = np.zeros((n, n), dtype=int)
    qr = {1, 2, 4}
    for i in range(n):
        for j in range(n):
            if i != j and ((j - i) % n) in qr:
                A[i][j] = 1
    betti, odims = path_betti_numbers(A, n)
    chi = sum((-1)**p * betti[p] for p in range(n))
    total_betti = sum(betti)
    print(f"Paley T_7:")
    print(f"  betti = {betti}")
    print(f"  omega_dims = {odims}")
    print(f"  Euler char chi = {chi}")
    print(f"  total Betti = {total_betti}")
    print(f"  Poincare poly P(t) = {' + '.join(f'{betti[p]}*t^{p}' for p in range(n) if betti[p])}")

    # Part 3: Omega dims comparison across regular classes
    print("\n" + "=" * 70)
    print("PART 3: Omega dims for the 3 regular classes at n=7")
    print("=" * 70)

    # Generate one from each class
    classes = {}
    for trial in range(100000):
        A = random_tournament(n, rng)
        scores = tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))
        if scores != (3,3,3,3,3,3,3):
            continue
        H = ham_path_count(A, n)
        if H not in classes:
            betti, odims = path_betti_numbers(A, n)
            classes[H] = (list(betti), list(odims), A.copy())
            print(f"  H={H}: betti={betti}, omega={odims}")
        if len(classes) == 3:
            break

    # Part 4: Detailed analysis of WHY beta_4 = 6 for Paley
    print("\n" + "=" * 70)
    print("PART 4: Detailed ker/im analysis for Paley T_7 at p=4")
    print("=" * 70)

    A = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(n):
            if i != j and ((j - i) % n) in qr:
                A[i][j] = 1

    allowed = {}
    for p in range(-1, n + 1):
        allowed[p] = [] if p < 0 else enumerate_allowed_paths(A, n, p)

    omega = {}
    for p in range(n + 1):
        omega[p] = compute_omega_basis(A, n, p, allowed[p], allowed[p-1])

    for p in range(n):
        dim_om = omega[p].shape[1] if omega[p].ndim == 2 else 0
        if dim_om == 0:
            print(f"  p={p}: dim(Omega)=0, beta=0")
            continue

        # ker(d_p)
        bd_p = build_boundary_matrix(allowed[p], allowed[p-1])
        bd_p_om = bd_p @ omega[p]
        if bd_p_om.size > 0:
            sv = np.linalg.svd(bd_p_om, compute_uv=False)
            rk = sum(s > 1e-8 for s in sv)
        else:
            rk = 0
        ker = dim_om - rk

        # im(d_{p+1})
        if p + 1 < n + 1:
            dim_omp1 = omega[p+1].shape[1] if omega[p+1].ndim == 2 else 0
            if dim_omp1 > 0:
                bd1 = build_boundary_matrix(allowed[p+1], allowed[p])
                bd1_om = bd1 @ omega[p+1]
                sv1 = np.linalg.svd(bd1_om, compute_uv=False)
                im = sum(s > 1e-8 for s in sv1)
            else:
                im = 0
        else:
            im = 0

        beta = ker - im
        print(f"  p={p}: dim(Omega)={dim_om}, rank(d)={rk}, ker={ker}, im(d+1)={im}, beta={beta}")

    # Part 5: Compare with H=171 (other regular)
    print("\n" + "=" * 70)
    print("PART 5: Same analysis for H=171 regular tournament")
    print("=" * 70)
    if 171 in classes:
        A = classes[171][2]
        allowed = {}
        for p in range(-1, n + 1):
            allowed[p] = [] if p < 0 else enumerate_allowed_paths(A, n, p)

        omega = {}
        for p in range(n + 1):
            omega[p] = compute_omega_basis(A, n, p, allowed[p], allowed[p-1])

        for p in range(n):
            dim_om = omega[p].shape[1] if omega[p].ndim == 2 else 0
            if dim_om == 0:
                print(f"  p={p}: dim(Omega)=0, beta=0")
                continue

            bd_p = build_boundary_matrix(allowed[p], allowed[p-1])
            bd_p_om = bd_p @ omega[p]
            if bd_p_om.size > 0:
                sv = np.linalg.svd(bd_p_om, compute_uv=False)
                rk = sum(s > 1e-8 for s in sv)
            else:
                rk = 0
            ker = dim_om - rk

            if p + 1 < n + 1:
                dim_omp1 = omega[p+1].shape[1] if omega[p+1].ndim == 2 else 0
                if dim_omp1 > 0:
                    bd1 = build_boundary_matrix(allowed[p+1], allowed[p])
                    bd1_om = bd1 @ omega[p+1]
                    sv1 = np.linalg.svd(bd1_om, compute_uv=False)
                    im = sum(s > 1e-8 for s in sv1)
                else:
                    im = 0
            else:
                im = 0

            beta = ker - im
            print(f"  p={p}: dim(Omega)={dim_om}, rank(d)={rk}, ker={ker}, im(d+1)={im}, beta={beta}")

if __name__ == '__main__':
    main()
