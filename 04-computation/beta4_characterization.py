"""
Characterize which tournaments have beta_4 > 0 in GLMY path homology.

Known: Paley T_7 has beta_4 = 6 (THM-099). Random n=7 tournaments give beta_4 = 0.
Question: Is beta_4 > 0 only for Paley / H-maximizers, or broader?

Strategy:
1. Test all regular n=7 tournaments (score (3,3,3,3,3,3,3))
2. Test near-regular n=7 tournaments
3. Test the circulant S={1,3} tournament at n=7
4. Look for structural characterization
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

def count_cycles(A, n, k):
    """Count directed k-cycles."""
    count = 0
    for subset in combinations(range(n), k):
        for perm in permutations(subset):
            if all(A[perm[i]][perm[(i+1) % k]] for i in range(k)):
                count += 1
    return count // k  # each cycle counted k times

def ham_path_count(A, n):
    """Count Hamiltonian paths via DP."""
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

def is_isomorphic(A, B, n):
    """Quick check via score sequence + cycle counts."""
    sa = tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))
    sb = tuple(sorted(sum(B[i][j] for j in range(n)) for i in range(n)))
    if sa != sb: return False
    ca = count_cycles(A, n, 3)
    cb = count_cycles(B, n, 3)
    return ca == cb  # rough check

def main():
    rng = np.random.RandomState(42)
    n = 7

    # Part 1: Paley T_7
    print("=" * 70)
    print("PART 1: Paley T_7 (known beta_4 = 6)")
    print("=" * 70)
    A = np.zeros((n, n), dtype=int)
    qr = {1, 2, 4}
    for i in range(n):
        for j in range(n):
            if i != j and ((j - i) % n) in qr:
                A[i][j] = 1
    betti, odims = path_betti_numbers(A, n)
    H = ham_path_count(A, n)
    c3 = count_cycles(A, n, 3)
    c5 = count_cycles(A, n, 5)
    print(f"Paley T_7: H={H}, c3={c3}, c5={c5}")
    print(f"  betti = {betti}")
    print(f"  omega_dims = {odims}")

    # Part 2: All circulant tournaments at n=7
    print("\n" + "=" * 70)
    print("PART 2: All circulant tournaments at n=7")
    print("=" * 70)

    for S_bits in range(1, 2**3):
        S = set()
        for k in range(1, 4):
            if S_bits & (1 << (k-1)):
                S.add(k)
        # Must be a valid tournament: k in S => n-k not in S
        valid = True
        for k in S:
            if (n - k) % n in S:
                valid = False
                break
        if not valid or len(S) != 3:  # regular = out-degree 3
            continue

        A = np.zeros((n, n), dtype=int)
        for i in range(n):
            for j in range(n):
                if i != j and ((j - i) % n) in S:
                    A[i][j] = 1
        # Verify tournament
        ok = all(A[i][j] + A[j][i] == 1 for i in range(n) for j in range(i+1, n))
        if not ok: continue

        betti, odims = path_betti_numbers(A, n)
        H = ham_path_count(A, n)
        c3 = count_cycles(A, n, 3)
        c5 = count_cycles(A, n, 5)
        is_paley = (S == qr)
        print(f"Circulant S={S}: H={H}, c3={c3}, c5={c5}, Paley={is_paley}")
        print(f"  betti = {betti}")

    # Part 3: Sample regular n=7 tournaments
    print("\n" + "=" * 70)
    print("PART 3: Regular n=7 tournaments (score 3,3,3,3,3,3,3)")
    print("=" * 70)

    n_regular = 0
    n_beta4_pos = 0
    beta4_vals = defaultdict(int)
    beta4_pos_examples = []

    for trial in range(20000):
        A = random_tournament(n, rng)
        scores = tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))
        if scores != (3, 3, 3, 3, 3, 3, 3):
            continue
        n_regular += 1
        betti, odims = path_betti_numbers(A, n)
        b4 = betti[4]
        beta4_vals[b4] += 1
        if b4 > 0:
            n_beta4_pos += 1
            H = ham_path_count(A, n)
            c3 = count_cycles(A, n, 3)
            c5 = count_cycles(A, n, 5)
            if len(beta4_pos_examples) < 5:
                beta4_pos_examples.append((H, c3, c5, list(betti), list(odims)))
            print(f"  FOUND beta_4={b4}: H={H}, c3={c3}, c5={c5}, betti={betti}")

        if n_regular % 50 == 0:
            print(f"  {n_regular} regular so far, {n_beta4_pos} with beta_4>0")

    print(f"\nRegular n=7: {n_regular} tested, {n_beta4_pos} with beta_4>0")
    print(f"beta_4 distribution: {dict(sorted(beta4_vals.items()))}")

    # Part 4: Does beta_4 > 0 correlate with H value?
    print("\n" + "=" * 70)
    print("PART 4: beta_4 by H value among regular tournaments")
    print("=" * 70)

    # The three rigid classes at n=7 regular:
    # H=189 (Paley/BIBD): 240 tours
    # H=175 (LTT): 720 tours
    # H=171 (other): 1680 tours

    h_to_beta4 = defaultdict(lambda: defaultdict(int))
    n_tested = 0
    for trial in range(50000):
        A = random_tournament(n, rng)
        scores = tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))
        if scores != (3, 3, 3, 3, 3, 3, 3):
            continue
        n_tested += 1
        H = ham_path_count(A, n)
        betti, _ = path_betti_numbers(A, n)
        b4 = betti[4]
        h_to_beta4[H][b4] += 1

        if n_tested % 100 == 0:
            print(f"  {n_tested} regular tested...")
        if n_tested >= 500:
            break

    print(f"\nH -> beta_4 distribution:")
    for H_val in sorted(h_to_beta4.keys()):
        print(f"  H={H_val}: {dict(sorted(h_to_beta4[H_val].items()))}")

if __name__ == '__main__':
    main()
