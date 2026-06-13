"""
Omega parity structure: What determines the parity of dim(Omega_p)?

Known:
- dim(Omega_0) = n (always)
- dim(Omega_1) = C(n,2) (always)
- dim(Omega_2) depends on tournament but has specific constraints

Questions:
1. Is dim(Omega_p) always determined mod 2 by n and p?
2. What is sum (-1)^p dim(Omega_p) mod 2?
3. Is there a mod-2 formula for dim(Omega_p)?
4. How does the PARITY of dim(Omega_p) relate to beta_p?

Also: the H(T_4) = 2*c3 + 1 identity suggests that mod-2 structure
of path counts is clean. Let's explore.
"""
import numpy as np
from math import comb
from collections import defaultdict
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
    print("OMEGA PARITY STRUCTURE: dim(Omega_p) mod 2")
    print("=" * 70)

    # Part 1: Parity of dim(Omega_p) at each n
    for n in range(3, 8):
        if n <= 6:
            total = 2**(n*(n-1)//2)
            N = total
        else:
            N = 300

        rng = np.random.RandomState(42 + n)

        parity_dist = defaultdict(lambda: defaultdict(int))
        omega_mod2 = defaultdict(list)
        cnt = 0

        for trial in range(N):
            if n <= 6:
                A = bits_to_adj(trial, n)
            else:
                A = random_tournament(n, rng)

            allowed = {}
            for p in range(-1, n):
                allowed[p] = [] if p < 0 else enumerate_allowed_paths(A, n, p)

            od_vec = []
            for p in range(n):
                od = compute_omega_dim(A, n, p, allowed[p], allowed[p-1])
                od_vec.append(od)
                parity_dist[p][od % 2] += 1

            omega_mod2_key = tuple(od % 2 for od in od_vec)
            omega_mod2[omega_mod2_key] = omega_mod2.get(omega_mod2_key, 0) + 1
            cnt += 1

        print(f"\n  n={n} (d={n-1}), {cnt} samples:")
        print(f"    dim(Omega_p) mod 2 distribution:")
        for p in range(n):
            even = parity_dist[p].get(0, 0)
            odd = parity_dist[p].get(1, 0)
            c = comb(n, p+1)
            print(f"      p={p}: even={even} ({100*even/cnt:.1f}%), odd={odd} ({100*odd/cnt:.1f}%), "
                  f"C(n,p+1)={c} ({'even' if c%2==0 else 'odd'})")

        print(f"    dim(Omega_p) mod 2 vector distribution:")
        for vec, count in sorted(omega_mod2.items(), key=lambda x: -x[1])[:5]:
            pct = 100 * count / cnt
            print(f"      {list(vec)}: {count} ({pct:.1f}%)")

    # Part 2: |A_p| mod 2 — is this always determined?
    # Since |A_p| = sum_{|S|=p+1} H(T[S]) and each H is ODD (Redei),
    # |A_p| mod 2 = C(n, p+1) mod 2 (sum of C(n,p+1) odd numbers)
    print("\n" + "=" * 70)
    print("Part 2: |A_p| mod 2 = C(n,p+1) mod 2 (from local Redei)")
    print("Proof: |A_p| = sum of C(n,p+1) ODD numbers => |A_p| = C(n,p+1) mod 2")
    print("=" * 70)

    for n in range(3, 8):
        rng = np.random.RandomState(42 + n)
        N = min(200, 2**(n*(n-1)//2))

        violations = 0
        cnt = 0
        for trial in range(N):
            if n <= 5:
                A = bits_to_adj(trial, n)
            else:
                A = random_tournament(n, rng)

            for p in range(n):
                ap = enumerate_allowed_paths(A, n, p)
                c = comb(n, p+1)
                if len(ap) % 2 != c % 2:
                    violations += 1
            cnt += 1

        print(f"  n={n}: {violations} violations out of {cnt*n} checks ({'PASS' if violations == 0 else 'FAIL'})")

    # Part 3: rank(constraint_matrix) mod 2
    # dim(Omega_p) = |A_p| - rank(P_p)
    # dim(Omega_p) mod 2 = (|A_p| - rank(P_p)) mod 2 = (C(n,p+1) - rank(P_p)) mod 2
    # So the parity of dim(Omega_p) is determined by C(n,p+1) and rank(P_p) mod 2
    print("\n" + "=" * 70)
    print("Part 3: rank(P_p) mod 2 determines dim(Omega_p) mod 2")
    print("dim(Omega_p) mod 2 = (C(n,p+1) + rank(P_p)) mod 2")
    print("=" * 70)

    for n in range(3, 8):
        rng = np.random.RandomState(42 + n)
        N = min(200, 2**(n*(n-1)//2))

        rank_parity = defaultdict(lambda: defaultdict(int))

        for trial in range(N):
            if n <= 5:
                A = bits_to_adj(trial, n)
            else:
                A = random_tournament(n, rng)

            allowed = {}
            for p in range(-1, n):
                allowed[p] = [] if p < 0 else enumerate_allowed_paths(A, n, p)

            for p in range(1, n):
                od = compute_omega_dim(A, n, p, allowed[p], allowed[p-1])
                rank_p = len(allowed[p]) - od
                rank_parity[p][rank_p % 2] += 1

        print(f"\n  n={n}:")
        for p in range(1, n):
            even = rank_parity[p].get(0, 0)
            odd = rank_parity[p].get(1, 0)
            total_checks = even + odd
            if total_checks > 0:
                print(f"    p={p}: rank even={even} ({100*even/total_checks:.1f}%), "
                      f"rank odd={odd} ({100*odd/total_checks:.1f}%)")

    # Part 4: The key question — is rank(P_2) mod 2 always equal to something simple?
    # rank(P_2) = |A_2| - dim(Omega_2)
    # |A_2| = C(n,3) + 2*c3 => |A_2| mod 2 = C(n,3) mod 2
    # So rank(P_2) mod 2 = (C(n,3) - dim(Omega_2)) mod 2
    print("\n" + "=" * 70)
    print("Part 4: Is rank(P_2) mod 2 constant for each n?")
    print("=" * 70)

    for n in range(3, 8):
        if n <= 6:
            total = 2**(n*(n-1)//2)
            N = total
        else:
            N = 500
        rng = np.random.RandomState(42 + n)

        rank2_parities = set()
        for trial in range(N):
            if n <= 6:
                A = bits_to_adj(trial, n)
            else:
                A = random_tournament(n, rng)

            a2 = enumerate_allowed_paths(A, n, 2)
            a1 = enumerate_allowed_paths(A, n, 1)
            od = compute_omega_dim(A, n, 2, a2, a1)
            rank2 = len(a2) - od
            rank2_parities.add(rank2 % 2)

        constant = len(rank2_parities) == 1
        print(f"  n={n}: rank(P_2) mod 2 in {rank2_parities} - {'CONSTANT' if constant else 'VARIES'}")

if __name__ == '__main__':
    main()
