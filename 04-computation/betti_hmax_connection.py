"""
Test the conjecture: nonzero higher Betti numbers <=> H-maximizer within score class.

At n=6:
  - beta_3=1 iff H = max within score class (confirmed for both score classes)

At n=7:
  - beta_4=6 iff H=189 (Paley/BIBD class = global max among regular)

At n=5: Test same pattern.

Refined conjecture: beta_{n-3} > 0 iff H = max within score class.
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

def main():
    # Part 1: n=5 exhaustive
    print("=" * 70)
    print("PART 1: n=5 exhaustive — Betti vs H-max")
    print("=" * 70)
    n = 5
    total = 2**(n*(n-1)//2)

    # Collect data
    data = []
    for bits in range(total):
        A = bits_to_adj(bits, n)
        scores = tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))
        H = ham_path_count(A, n)
        betti, odims = path_betti_numbers(A, n)
        data.append((bits, scores, H, tuple(betti), tuple(odims)))

    # Find max H per score class
    max_H_by_score = defaultdict(int)
    for _, scores, H, _, _ in data:
        max_H_by_score[scores] = max(max_H_by_score[scores], H)

    # Check: does any higher Betti > 0 correlate with max H?
    print(f"\nn={n}: {total} tournaments, {len(set(s for _, s, _, _, _ in data))} score classes")

    # All unique Betti vectors
    betti_dist = defaultdict(int)
    for _, _, _, betti, _ in data:
        betti_dist[betti] += 1
    print(f"\nBetti vector distribution:")
    for bv, count in sorted(betti_dist.items(), key=lambda x: -x[1]):
        print(f"  {list(bv)}: {count}")

    # Check beta_2 (which we proved is 0)
    # At n=5, beta_{n-3} = beta_2 = 0 always

    # Check: is beta_1 = 1 related to H-max?
    print(f"\nbeta_1 = 1 vs H-max within score class:")
    for scores in sorted(set(s for _, s, _, _, _ in data)):
        max_H = max_H_by_score[scores]
        with_b1 = [(bits, H, betti) for bits, s, H, betti, _ in data if s == scores and betti[1] == 1]
        without_b1 = [(bits, H, betti) for bits, s, H, betti, _ in data if s == scores and betti[1] == 0]
        if with_b1:
            H_b1 = set(H for _, H, _ in with_b1)
            H_no = set(H for _, H, _ in without_b1)
            is_max = max_H in H_b1
            print(f"  Score {scores}: max_H={max_H}, beta_1=1 H={sorted(H_b1)}, beta_1=0 H={sorted(H_no)}, max in b1? {is_max}")

    # Part 2: n=4 exhaustive
    print("\n" + "=" * 70)
    print("PART 2: n=4 exhaustive")
    print("=" * 70)
    n = 4
    total = 2**(n*(n-1)//2)
    data4 = []
    for bits in range(total):
        A = bits_to_adj(bits, n)
        scores = tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))
        H = ham_path_count(A, n)
        betti, odims = path_betti_numbers(A, n)
        data4.append((bits, scores, H, tuple(betti), tuple(odims)))

    betti_dist4 = defaultdict(int)
    for _, _, _, betti, _ in data4:
        betti_dist4[betti] += 1
    print(f"\nn={n}: Betti vector distribution:")
    for bv, count in sorted(betti_dist4.items(), key=lambda x: -x[1]):
        print(f"  {list(bv)}: {count}")

    # Part 3: n=3 exhaustive
    print("\n" + "=" * 70)
    print("PART 3: n=3 exhaustive")
    print("=" * 70)
    n = 3
    total = 2**(n*(n-1)//2)
    for bits in range(total):
        A = bits_to_adj(bits, n)
        scores = tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))
        H = ham_path_count(A, n)
        betti, odims = path_betti_numbers(A, n)
        print(f"  bits={bits}: scores={scores}, H={H}, betti={list(betti)}, omega={list(odims)}")

    # Part 4: Summary table
    print("\n" + "=" * 70)
    print("SUMMARY: Higher Betti and H-maximization")
    print("=" * 70)
    print("n | beta_{n-3} | relation to H-max?")
    print("3 | beta_0     | trivially 1 always")
    print("4 | beta_1     | = 1 iff has 3-cycle (= H > 1)")
    print("5 | beta_2     | always 0 (PROVED)")
    print("6 | beta_3     | = 1 iff H = max in score class (320/32768)")
    print("7 | beta_4     | = 6 iff H = 189 (Paley class, global max among regular)")
    print("")
    print("Pattern: beta_{n-3} > 0 picks out H-maximizers!")
    print("But beta_2 = 0 breaks the sequence at n=5...")
    print("")
    print("Alternative: beta_{2*floor((n-1)/2)} > 0 picks out H-maximizers?")
    print("n=3: beta_2 = 0 (transitive) or beta_2 = 0 (3-cycle)... no")
    print("n=6: beta_3 = beta_{2*2} = beta_4? No, beta_3")
    print("")
    print("Or: the HIGHEST nonzero even Betti (after beta_0)?")
    print("n=3: beta_0=1 always")
    print("n=5: beta_0=1 always, all even betti = 0")
    print("n=6: beta_0=1, beta_2=0, beta_4=0. Odd: beta_1,beta_3")
    print("n=7: beta_0=1, beta_2=0, beta_4=6 or 0, beta_6=0")

if __name__ == '__main__':
    main()
