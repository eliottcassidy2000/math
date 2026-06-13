"""
Poincare polynomial of the path complex: P(T,t) = sum dim(Omega_p) * t^p

This encodes all Omega dimensions. For transitive: P = sum C(n,p+1)*t^p = ((1+t)^n - 1)/t.

For a general tournament, P(T,t) is a "deformed" version of the simplex Poincare polynomial.

The RATIO P(T,t) / P(trans,t) is the "filling polynomial" — its coefficients
measure the relative inflation/deflation at each dimension.

Also: the independence polynomial I(Omega(T), x) satisfies H(T) = I(Omega(T), 2).
I(Omega(T), -1) = sum (-1)^k i_k = alt sum of independence numbers.

There's a deep connection:
- P(T,t) encodes the CHAIN COMPLEX structure
- I(Omega(T),x) encodes the CONFLICT GRAPH structure
- Both are attached to the same tournament T

Question: Is there a functional equation relating P and I?

Also: the Poincare SERIES P(T,t)/(1-t) or P(T,t)/(1+t) might have nice properties.
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

def compute_conflict_graph_IP(A, n):
    """Compute independence polynomial of Omega(T) conflict graph."""
    # Get 3-cycles as vertices of Omega graph
    cycles = []
    for a, b, c in combinations(range(n), 3):
        if A[a][b] and A[b][c] and A[c][a]:
            cycles.append(frozenset([a,b,c]))
        if A[a][c] and A[c][b] and A[b][a]:
            cycles.append(frozenset([a,b,c]))

    m = len(cycles)
    if m == 0: return [1]

    # Build conflict graph
    adj = [[False]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1, m):
            if cycles[i] & cycles[j]:
                adj[i][j] = adj[j][i] = True

    # Independence polynomial by inclusion-exclusion on independent sets
    # For small m, enumerate via bitmask
    if m > 20: return None  # too large

    ip = [0] * (m + 1)
    for mask in range(1 << m):
        bits = bin(mask).count('1')
        is_indep = True
        vs = [i for i in range(m) if mask & (1 << i)]
        for i in range(len(vs)):
            for j in range(i+1, len(vs)):
                if adj[vs[i]][vs[j]]:
                    is_indep = False
                    break
            if not is_indep:
                break
        if is_indep:
            ip[bits] += 1

    # Trim trailing zeros
    while ip and ip[-1] == 0:
        ip.pop()
    return ip

def main():
    print("=" * 70)
    print("POINCARE POLYNOMIAL: P(T,t) = sum dim(Omega_p) * t^p")
    print("Transitive: P = ((1+t)^n - 1)/t")
    print("=" * 70)

    # Part 1: Compute P(T,t) for specific tournaments
    print("\n--- Part 1: Poincare polynomials for named tournaments ---")

    for n in [3, 5, 7]:
        configs = []
        # Transitive
        A_trans = np.array([[1 if j > i else 0 for j in range(n)] for i in range(n)])
        configs.append(("Transitive", A_trans))
        # Paley
        if n in [3, 7]:
            A_paley = np.zeros((n, n), dtype=int)
            qr = {k*k % n for k in range(1, n)}
            for i in range(n):
                for j in range(n):
                    if i != j and ((j-i) % n) in qr:
                        A_paley[i][j] = 1
            configs.append(("Paley", A_paley))
        # Random
        rng = np.random.RandomState(42)
        A_rand = random_tournament(n, rng)
        configs.append(("Random", A_rand))

        for label, A in configs:
            allowed = {}
            for p in range(-1, n):
                allowed[p] = [] if p < 0 else enumerate_allowed_paths(A, n, p)

            od = []
            for p in range(n):
                d = compute_omega_dim(A, n, p, allowed[p], allowed[p-1])
                od.append(d)

            # Transitive baseline
            trans_od = [comb(n, p+1) for p in range(n)]

            # Ratios
            ratios = [od[p] / trans_od[p] if trans_od[p] > 0 else 0 for p in range(n)]

            # P(T, -1) = sum (-1)^p * dim(Omega_p) = chi(T)
            chi = sum((-1)**p * od[p] for p in range(n))

            # P(T, 1) = sum dim(Omega_p)
            P1 = sum(od)

            print(f"\n  {label} T_{n}:")
            print(f"    Omega dims: {od}")
            print(f"    Transitive: {trans_od}")
            print(f"    Ratios: {[f'{r:.3f}' for r in ratios]}")
            print(f"    P(T,-1) = chi = {chi}")
            print(f"    P(T,1) = sum dims = {P1}")
            print(f"    P(trans,1) = 2^n - 1 = {2**n - 1}")

    # Part 2: P(T,1) distribution (total Omega dimension)
    print("\n" + "=" * 70)
    print("Part 2: P(T,1) = sum dim(Omega_p) distribution")
    print("Transitive: P(trans,1) = 2^n - 1")
    print("=" * 70)

    for n in range(3, 8):
        rng = np.random.RandomState(42 + n)
        N = min(200, 2**(n*(n-1)//2))
        trans_total = 2**n - 1

        P1_vals = []
        for trial in range(N):
            if n <= 5:
                A = bits_to_adj(trial, n)
            else:
                A = random_tournament(n, rng)

            allowed = {}
            for p in range(-1, n):
                allowed[p] = [] if p < 0 else enumerate_allowed_paths(A, n, p)

            total = sum(compute_omega_dim(A, n, p, allowed[p], allowed[p-1]) for p in range(n))
            P1_vals.append(total)

        print(f"\n  n={n}: P(trans,1) = {trans_total}")
        print(f"    avg P(T,1) = {np.mean(P1_vals):.1f}")
        print(f"    min = {min(P1_vals)}, max = {max(P1_vals)}")
        print(f"    ratio avg/trans = {np.mean(P1_vals)/trans_total:.4f}")

    # Part 3: Independence polynomial vs Poincare polynomial
    print("\n" + "=" * 70)
    print("Part 3: I(Omega(T), x) vs P(T, t) — are they related?")
    print("=" * 70)

    for n in [5, 6, 7]:
        rng = np.random.RandomState(42 + n)
        N = min(30, 2**(n*(n-1)//2))

        for trial in range(min(5, N)):
            if n <= 5:
                A = bits_to_adj(trial, n)
            else:
                A = random_tournament(n, rng)

            # Poincare polynomial
            allowed = {}
            for p in range(-1, n):
                allowed[p] = [] if p < 0 else enumerate_allowed_paths(A, n, p)
            od = [compute_omega_dim(A, n, p, allowed[p], allowed[p-1]) for p in range(n)]

            # Independence polynomial
            ip = compute_conflict_graph_IP(A, n)

            # H = I(Omega, 2)
            H = len(allowed[n-1]) if n-1 in allowed else 0
            if ip:
                I_at_2 = sum(c * 2**k for k, c in enumerate(ip))
            else:
                I_at_2 = "too large"

            # P(-1) = chi, I(-1) = alternating sum of i_k
            chi = sum((-1)**p * od[p] for p in range(n))
            if ip:
                I_at_m1 = sum(c * (-1)**k for k, c in enumerate(ip))
            else:
                I_at_m1 = "too large"

            print(f"\n  n={n}, trial {trial}:")
            print(f"    P(t) coeffs = {od}, P(-1) = {chi}")
            if ip:
                print(f"    I(x) coeffs = {ip}, I(-1) = {I_at_m1}")
                print(f"    H = {H}, I(2) = {I_at_2}, match: {H == I_at_2}")

    # Part 4: Does P(T,1)/P(trans,1) approach a limit?
    print("\n" + "=" * 70)
    print("Part 4: P(T,1)/(2^n - 1) ratio scaling")
    print("=" * 70)

    for n in range(3, 9):
        rng = np.random.RandomState(42 + n)
        N = min(200, 2**(n*(n-1)//2))
        trans_total = 2**n - 1

        ratios = []
        for trial in range(N):
            if n <= 5:
                A = bits_to_adj(trial, n)
            else:
                A = random_tournament(n, rng)

            allowed = {}
            for p in range(-1, n):
                allowed[p] = [] if p < 0 else enumerate_allowed_paths(A, n, p)

            total = sum(compute_omega_dim(A, n, p, allowed[p], allowed[p-1]) for p in range(n))
            ratios.append(total / trans_total)

        print(f"  n={n}: avg ratio = {np.mean(ratios):.4f}, "
              f"range [{min(ratios):.4f}, {max(ratios):.4f}]")

if __name__ == '__main__':
    main()
