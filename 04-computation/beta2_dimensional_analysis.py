"""
Dimensional analysis of beta_2=0.

beta_2 = ker(d_2|Omega_2) - im(d_3|Omega_3)

If we can show dim(ker d_2) = dim(im d_3) always, beta_2=0 follows.

key quantities:
  dim(ker d_2) = dim(Omega_2) - rank(d_2|Omega_2)
  dim(im d_3) = rank(d_3|Omega_3)

So beta_2 = 0 iff dim(Omega_2) - rank(d_2) = rank(d_3).

Let's track these numbers and see if there's a formula.

Also: since beta_0=1 always and beta_1 in {0,1}:
  rank(d_1) = dim(Omega_1) - ker(d_1) = C(n,2) - ker(d_1)
  beta_0 = ker(d_0) - im(d_1) = n - rank(d_1) = 1
  So rank(d_1) = n-1
  ker(d_1) = C(n,2) - (n-1)

  im(d_2) = rank(d_2|Omega_2) = dim(Omega_2) - ker(d_2)
  beta_1 = ker(d_1) - im(d_2) = [C(n,2)-(n-1)] - [dim(Omega_2) - ker(d_2)]

  If beta_2=0: ker(d_2) = im(d_3) = rank(d_3)
  So: im(d_2) = dim(Omega_2) - rank(d_3)
  And: beta_1 = [C(n,2)-(n-1)] - dim(Omega_2) + rank(d_3)

This gives: beta_1 + dim(Omega_2) = C(n,2) - n + 1 + rank(d_3)

We know beta_1 in {0,1} and dim(Omega_2) varies.
So rank(d_3) = dim(Omega_2) - C(n,2) + n - 1 + beta_1
            = dim(Omega_2) - C(n,2) + n - 1 + beta_1

Can we verify this formula and see if rank(d_3) is always this?
"""

import numpy as np
from math import comb
from itertools import combinations

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

def compute_omega_dim(A, n, p, allowed_p, allowed_pm1):
    dim_Ap = len(allowed_p)
    if dim_Ap == 0: return 0
    if p == 0: return dim_Ap
    allowed_pm1_set = set(allowed_pm1)
    non_allowed = {}
    na_count = 0
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if len(set(face)) == len(face) and face not in allowed_pm1_set:
                if face not in non_allowed:
                    non_allowed[face] = na_count
                    na_count += 1
    if na_count == 0: return dim_Ap
    P = np.zeros((na_count, dim_Ap))
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if face in non_allowed:
                P[non_allowed[face], j] += sign
    sv = np.linalg.svd(P, compute_uv=False)
    rank = sum(s > 1e-10 for s in sv)
    return dim_Ap - rank

def main():
    print("=" * 70)
    print("DIMENSIONAL ANALYSIS OF BETA_2 = 0")
    print("=" * 70)

    # Track: dim(Omega_2), |A_2|, #constraints, rank(constraints),
    # and the derived formula for rank(d_3)
    for n in [4, 5, 6]:
        print(f"\n--- n={n} ---")
        N = 2**(n*(n-1)//2)

        for trial in range(N):
            A = bits_to_adj(trial, n)

            allowed = {}
            for p in range(-1, n):
                allowed[p] = [] if p < 0 else enumerate_allowed_paths(A, n, p)

            # Omega dimensions
            od = {}
            for p in range(min(n, 5)):
                od[p] = compute_omega_dim(A, n, p, allowed[p], allowed[p-1])

            # |A_p| - dim(Omega_p) = rank of constraint matrix at level p
            rank_constraint = {}
            for p in range(1, min(n, 5)):
                rank_constraint[p] = len(allowed[p]) - od[p]

            # Formula check: if beta_2=0 and beta_1 in {0,1},
            # then rank(d_3) = dim(Omega_2) - C(n,2) + n - 1 + beta_1
            # But rank(d_3) = dim(Omega_3) - ker(d_3)
            # and beta_3 = ker(d_3) - im(d_4)

            # Actually let me just compute the Euler characteristic
            # chi = sum (-1)^p dim(Omega_p) (this is always 1 for tournaments)
            chi_omega = sum((-1)**p * od.get(p, 0) for p in range(n))

            # For tournament: chi_omega = 1 always? Let me check.
            # No — chi_omega is NOT the Betti number chi. It's the OMEGA chi.
            # chi = sum (-1)^p dim(Omega_p) = sum (-1)^p beta_p (by rank-nullity chain)
            # Wait, that's only for simplicial complexes.
            # For chain complexes: chi = sum (-1)^p dim(C_p) = sum (-1)^p beta_p
            # So chi_omega = sum (-1)^p dim(Omega_p) = sum (-1)^p beta_p.

            if trial < 3 or (n == 5 and trial < 10):
                c3 = sum(1 for i,j,k in combinations(range(n), 3)
                         if (A[i][j] and A[j][k] and A[k][i]) or
                            (A[i][k] and A[k][j] and A[j][i]))

                od_list = [od.get(p, 0) for p in range(min(n, 5))]
                rk_list = [rank_constraint.get(p, 0) for p in range(1, min(n, 5))]
                print(f"  trial {trial}: c3={c3}, Omega={od_list}, "
                      f"rank_constr={rk_list}, chi_omega={chi_omega}")

        # Does chi_omega = 1 always?
        all_chi1 = True
        for trial in range(N):
            A = bits_to_adj(trial, n)
            allowed = {}
            for p in range(-1, n):
                allowed[p] = [] if p < 0 else enumerate_allowed_paths(A, n, p)
            od = {}
            for p in range(min(n, 5)):
                od[p] = compute_omega_dim(A, n, p, allowed[p], allowed[p-1])
            chi = sum((-1)**p * od.get(p, 0) for p in range(n))
            if chi != 1:
                all_chi1 = False
                break

        print(f"  chi_omega = 1 for ALL {N} tournaments: {all_chi1}")

    # Part 2: The key relationship
    # chi = sum (-1)^p dim(Omega_p) = 1 (if true)
    # This means: n - C(n,2) + dim(Omega_2) - dim(Omega_3) + dim(Omega_4) - ... = 1
    # => dim(Omega_2) - dim(Omega_3) + dim(Omega_4) - ... = 1 - n + C(n,2)
    #    = 1 - n + n(n-1)/2 = (n^2 - 3n + 2)/2 = (n-1)(n-2)/2

    print("\n" + "=" * 70)
    print("VERIFICATION: chi(Omega) = 1 implies constraints on Omega dims")
    print("=" * 70)
    for n in [4, 5, 6, 7]:
        target = (n-1)*(n-2)//2
        print(f"\n  n={n}: sum_{p>=2} (-1)^p dim(Omega_p) should be {target}")

        rng = np.random.RandomState(42 + n)
        N = min(2**(n*(n-1)//2), 200)
        passes = 0
        for trial in range(N):
            if n <= 6:
                A = bits_to_adj(trial, n)
            else:
                A = random_tournament(n, rng)

            allowed = {}
            for p in range(-1, n):
                allowed[p] = [] if p < 0 else enumerate_allowed_paths(A, n, p)
            od = {}
            for p in range(n):
                od[p] = compute_omega_dim(A, n, p, allowed[p], allowed[p-1])

            # Check: n - C(n,2) + sum_{p>=2} (-1)^p od[p] = 1
            rest = sum((-1)**p * od.get(p, 0) for p in range(2, n))
            if n - comb(n, 2) + rest == 1:
                passes += 1
            elif trial < 5:
                print(f"    FAIL: trial {trial}, n-C(n,2)+rest = {n - comb(n,2) + rest}")

        print(f"    {passes}/{N} pass")

    # Part 3: Does chi_omega=1 hold at n=8 despite beta_4>0?
    # chi_omega = sum (-1)^p dim(Omega_p)
    # Since dim(Omega_p) = sum_{q<=p} (-1)^{p-q} dim(A_q) * (something...)
    # Actually chi_omega = sum (-1)^p beta_p.
    # At n=8 with beta_4=1: chi_omega = 1 - 0 + 0 - 0 + 1 - 0 + 0 - 0 = 2 ≠ 1!
    # But wait — chi_omega computed from Omega dims might still be 1
    # because chi = sum (-1)^p dim(Omega_p) EQUALS sum (-1)^p beta_p.
    # So if chi ≠ 1 at n=8 (which we saw: chi=2 for beta_4=1), then
    # chi_omega = sum (-1)^p dim(Omega_p) ≠ 1 at n=8!

    print("\n" + "=" * 70)
    print("CHI_OMEGA AT n=8")
    print("=" * 70)
    n = 8
    rng = np.random.RandomState(43)  # seed with beta_4>0 cases
    for trial in range(200):
        A = random_tournament(n, rng)
        allowed = {}
        for p in range(-1, n):
            allowed[p] = [] if p < 0 else enumerate_allowed_paths(A, n, p)
        od = {}
        for p in range(n):
            od[p] = compute_omega_dim(A, n, p, allowed[p], allowed[p-1])

        chi = sum((-1)**p * od.get(p, 0) for p in range(n))
        if chi != 1:
            od_list = [od.get(p, 0) for p in range(n)]
            print(f"  trial {trial}: chi_omega = {chi}, Omega = {od_list}")

if __name__ == '__main__':
    main()
