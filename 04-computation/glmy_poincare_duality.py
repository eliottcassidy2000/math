"""
Investigate Poincare duality in GLMY path homology of tournaments.

For Paley T_7, Omega dimensions are palindromic: 7, 21, 42, 63, 63, 42, 21.
And Betti numbers are palindromic: 1, 0, 0, 0, 6, 0, 0.

Questions:
1. Is dim(Omega_p) always palindromic for tournaments?
2. Is the Betti vector always palindromic?
3. What is the Euler characteristic chi = sum (-1)^p beta_p?
4. Does Poincare duality hold: beta_p = beta_{n-1-p}?
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

def main():
    rng = np.random.RandomState(42)

    # PART 1: Poincare duality at n=5 (exhaustive)
    print("=" * 70)
    print("PART 1: Poincare duality at n=5 (exhaustive)")
    print("=" * 70)

    n = 5
    total = 2**(n*(n-1)//2)
    pd_holds = 0
    pd_fails = 0
    omega_palindrome = 0
    omega_not_palindrome = 0
    chi_dist = defaultdict(int)

    for bits in range(total):
        A = bits_to_adj(bits, n)
        betti, omega_dims = path_betti_numbers(A, n)

        # Check Poincare duality: beta_p = beta_{n-1-p}
        pd = all(betti[p] == betti[n-1-p] for p in range(n))
        if pd:
            pd_holds += 1
        else:
            pd_fails += 1

        # Check Omega palindrome
        od = omega_dims
        op = all(od[p] == od[n-1-p] for p in range(n))
        if op:
            omega_palindrome += 1
        else:
            omega_not_palindrome += 1

        # Euler characteristic
        chi = sum((-1)**p * betti[p] for p in range(n))
        chi_dist[chi] += 1

    print(f"n={n}: {total} tournaments")
    print(f"  Betti palindrome (PD): {pd_holds}, not: {pd_fails}")
    print(f"  Omega palindrome: {omega_palindrome}, not: {omega_not_palindrome}")
    print(f"  Euler char chi: {dict(sorted(chi_dist.items()))}")

    # PART 2: n=6 (exhaustive for PD check, smaller batch for full Betti)
    print("\n" + "=" * 70)
    print("PART 2: n=6 sampled")
    print("=" * 70)

    n = 6
    nsamp = 500
    pd_holds = 0
    pd_fails = 0
    omega_palindrome = 0
    omega_not_palindrome = 0
    chi_dist = defaultdict(int)

    for _ in range(nsamp):
        A = random_tournament(n, rng)
        betti, omega_dims = path_betti_numbers(A, n)

        pd = all(betti[p] == betti[n-1-p] for p in range(n))
        if pd: pd_holds += 1
        else:
            pd_fails += 1
            if pd_fails <= 3:
                print(f"  PD FAILS: betti={betti}, omega={omega_dims}")

        od = omega_dims
        op = all(od[p] == od[n-1-p] for p in range(n))
        if op: omega_palindrome += 1
        else: omega_not_palindrome += 1

        chi = sum((-1)**p * betti[p] for p in range(n))
        chi_dist[chi] += 1

    print(f"n={n} ({nsamp} random):")
    print(f"  Betti palindrome (PD): {pd_holds}, not: {pd_fails}")
    print(f"  Omega palindrome: {omega_palindrome}, not: {omega_not_palindrome}")
    print(f"  Euler char chi: {dict(sorted(chi_dist.items()))}")

    # PART 3: Specific tournaments at n=7
    print("\n" + "=" * 70)
    print("PART 3: Specific n=7 tournaments")
    print("=" * 70)

    n = 7
    # Paley T_7
    A = np.zeros((n, n), dtype=int)
    qr = {1, 2, 4}
    for i in range(n):
        for j in range(n):
            if i != j and ((j - i) % n) in qr:
                A[i][j] = 1
    betti, omega_dims = path_betti_numbers(A, n)
    chi = sum((-1)**p * betti[p] for p in range(n))
    pd = all(betti[p] == betti[n-1-p] for p in range(n))
    op = all(omega_dims[p] == omega_dims[n-1-p] for p in range(n))
    print(f"Paley T_7: betti={betti}, omega={omega_dims}")
    print(f"  PD: {pd}, Omega palindrome: {op}, chi: {chi}")

    # Other circulant at n=7
    for S_list in [{1,2,3}, {1,3}]:
        A = np.zeros((n, n), dtype=int)
        for i in range(n):
            for j in range(n):
                if i != j and ((j-i) % n) in S_list:
                    A[i][j] = 1
        ok = all(A[i][j] + A[j][i] == 1 for i in range(n) for j in range(i+1, n))
        if not ok: continue
        betti, omega_dims = path_betti_numbers(A, n)
        chi = sum((-1)**p * betti[p] for p in range(n))
        pd = all(betti[p] == betti[n-1-p] for p in range(n))
        op = all(omega_dims[p] == omega_dims[n-1-p] for p in range(n))
        H = len([p for p in __import__('itertools').permutations(range(n))
                if all(A[p[i]][p[i+1]] for i in range(n-1))])
        print(f"\nCirculant S={S_list}: H={H}")
        print(f"  betti={betti}, omega={omega_dims}")
        print(f"  PD: {pd}, Omega palindrome: {op}, chi: {chi}")

    # Random n=7 tournaments
    print("\nRandom n=7 tournaments:")
    for trial in range(30):
        A = random_tournament(n, rng)
        betti, omega_dims = path_betti_numbers(A, n)
        chi = sum((-1)**p * betti[p] for p in range(n))
        pd = all(betti[p] == betti[n-1-p] for p in range(n))
        if not pd:
            print(f"  PD FAILS: betti={betti}")
            break
    else:
        print(f"  PD holds for all 30 samples")

    # PART 4: Euler characteristic analysis
    print("\n" + "=" * 70)
    print("PART 4: Euler characteristic formula")
    print("=" * 70)
    print("chi = sum (-1)^p beta_p. For Paley T_7: chi = 1+0-0+0+6-0+0 = 7 = n")
    print("Does chi = 1 or chi = n or something else?")

    for n in [5]:
        total = 2**(n*(n-1)//2)
        chi_vals = []
        for bits in range(total):
            A = bits_to_adj(bits, n)
            betti, _ = path_betti_numbers(A, n)
            chi = sum((-1)**p * betti[p] for p in range(n))
            chi_vals.append(chi)
        print(f"n={n}: chi values: {dict(sorted(defaultdict(int, ((c, chi_vals.count(c)) for c in set(chi_vals))).items()))}")

if __name__ == '__main__':
    main()
