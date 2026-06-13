"""
Betti rate scaling: P(beta_p > 0) vs dimension d = n-1.

From dimensional_meta_patterns data:
P(beta_1>0): 25%, 37.5%, 25%, 15.8%, 5.4%, 1.7% at n=3,...,8
P(beta_3>0): -, -, 0%, 0.98%, 8.7%, 18.7% at n=5,...,8

Questions:
1. Does P(beta_1>0) decay exponentially with n? Like C * exp(-alpha * n)?
2. Does P(beta_3>0) grow toward some limit?
3. What about P(beta_5>0) and higher?
4. Is there a universal scaling for P(beta_p>0) as a function of p/n?

For larger n (9, 10), we can only sample and compute beta_1 efficiently.
Computing beta_3 at n=9 is expensive but feasible.
"""
import numpy as np
from math import comb
from collections import defaultdict
from itertools import combinations

def random_tournament(n, rng):
    A = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

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

def betti_single(A, n, target_p):
    allowed = {}
    for p in [target_p - 1, target_p, target_p + 1]:
        if p < 0:
            allowed[p] = []
        else:
            allowed[p] = enumerate_allowed_paths(A, n, p)
    omega_p = compute_omega_basis(A, n, target_p, allowed[target_p], allowed[target_p-1])
    dim_om = omega_p.shape[1] if omega_p.ndim == 2 else 0
    if dim_om == 0:
        return 0
    bd_p = build_boundary_matrix(allowed[target_p], allowed[target_p-1])
    bd_p_om = bd_p @ omega_p
    if bd_p_om.size > 0:
        sv = np.linalg.svd(bd_p_om, compute_uv=False)
        rk = sum(s > 1e-8 for s in sv)
    else:
        rk = 0
    ker = dim_om - rk
    omega_p1 = compute_omega_basis(A, n, target_p+1, allowed[target_p+1], allowed[target_p])
    dim_om1 = omega_p1.shape[1] if omega_p1.ndim == 2 else 0
    if dim_om1 > 0:
        bd1 = build_boundary_matrix(allowed[target_p+1], allowed[target_p])
        bd1_om = bd1 @ omega_p1
        sv1 = np.linalg.svd(bd1_om, compute_uv=False)
        im = sum(s > 1e-8 for s in sv1)
    else:
        im = 0
    return ker - im

def main():
    print("=" * 70)
    print("BETTI RATE SCALING: P(beta_p > 0) vs dimension")
    print("=" * 70)

    # Part 1: beta_1 rate at n=3..10
    print("\n--- Part 1: P(beta_1 > 0) scaling ---")

    b1_rates = {}
    for n in range(3, 11):
        rng = np.random.RandomState(42 + n)
        if n <= 5:
            N = 2**(n*(n-1)//2)
        elif n <= 7:
            N = 500
        else:
            N = 300

        b1_pos = 0
        cnt = 0
        for trial in range(N):
            if n <= 5:
                A = bits_to_adj(trial, n)
            else:
                A = random_tournament(n, rng)

            b1 = betti_single(A, n, 1)
            if b1 > 0: b1_pos += 1
            cnt += 1

        rate = b1_pos / cnt
        b1_rates[n] = rate
        print(f"  n={n} (d={n-1}): P(beta_1>0) = {100*rate:.2f}% ({b1_pos}/{cnt})")

    # Fit exponential decay
    from math import log
    xs = sorted(b1_rates.keys())
    xs_fit = [n for n in xs if b1_rates[n] > 0]
    if len(xs_fit) >= 3:
        log_rates = [log(b1_rates[n]) for n in xs_fit]
        ns = [float(n) for n in xs_fit]
        # Linear fit: log(rate) = a*n + b => rate ~ exp(a*n + b)
        coeffs = np.polyfit(ns, log_rates, 1)
        print(f"\n  Exponential fit: P(beta_1>0) ~ exp({coeffs[0]:.3f}*n + {coeffs[1]:.3f})")
        print(f"  Half-life in n: {-log(2)/coeffs[0]:.2f}")
        for n in xs_fit:
            pred = np.exp(coeffs[0]*n + coeffs[1])
            print(f"    n={n}: actual={100*b1_rates[n]:.2f}%, predicted={100*pred:.2f}%")

    # Part 2: beta_3 rate at n=5..9
    print("\n--- Part 2: P(beta_3 > 0) scaling ---")

    b3_rates = {}
    for n in range(5, 10):
        rng = np.random.RandomState(42 + n)
        if n <= 5:
            N = 2**(n*(n-1)//2)
        elif n <= 7:
            N = 500
        elif n == 8:
            N = 200
        else:
            N = 100  # n=9 is expensive

        b3_pos = 0
        cnt = 0
        for trial in range(N):
            if n <= 5:
                A = bits_to_adj(trial, n)
            else:
                A = random_tournament(n, rng)

            b3 = betti_single(A, n, 3)
            if b3 > 0: b3_pos += 1
            cnt += 1

            if cnt % 50 == 0:
                print(f"    n={n}: {cnt}/{N} done, P(b3>0) so far = {100*b3_pos/cnt:.1f}%")

        rate = b3_pos / cnt
        b3_rates[n] = rate
        print(f"  n={n} (d={n-1}): P(beta_3>0) = {100*rate:.2f}% ({b3_pos}/{cnt})")

    # Part 3: beta_5 rate at n=7..9
    print("\n--- Part 3: P(beta_5 > 0) scaling ---")

    b5_rates = {}
    for n in range(7, 10):
        rng = np.random.RandomState(42 + n)
        N = 100 if n <= 8 else 50

        b5_pos = 0
        cnt = 0
        for trial in range(N):
            A = random_tournament(n, rng)
            b5 = betti_single(A, n, 5)
            if b5 > 0: b5_pos += 1
            cnt += 1

        rate = b5_pos / cnt
        b5_rates[n] = rate
        print(f"  n={n} (d={n-1}): P(beta_5>0) = {100*rate:.2f}% ({b5_pos}/{cnt})")

    # Part 4: Summary table
    print("\n" + "=" * 70)
    print("SUMMARY: Betti rates across dimensions")
    print("=" * 70)
    print(f"{'n':>3} | {'d':>3} | {'P(b1>0)':>10} | {'P(b3>0)':>10} | {'P(b5>0)':>10}")
    print(f"{'-'*3}-+-{'-'*3}-+-{'-'*10}-+-{'-'*10}-+-{'-'*10}")
    for n in range(3, 11):
        d = n - 1
        b1 = f"{100*b1_rates.get(n, float('nan')):.2f}%" if n in b1_rates else "-"
        b3 = f"{100*b3_rates.get(n, float('nan')):.2f}%" if n in b3_rates else "-"
        b5 = f"{100*b5_rates.get(n, float('nan')):.2f}%" if n in b5_rates else "-"
        print(f"{n:3d} | {d:3d} | {b1:>10} | {b3:>10} | {b5:>10}")

    # Part 5: Defect rate = P(chi != 1) = P(beta_1>0 OR beta_3>0 OR ...)
    # Since beta_1 and beta_3 are mutually exclusive and each <= 1:
    # P(chi != 1) = P(b1>0) + P(b3>0) (approximately, ignoring rare higher betti)
    print("\n--- Part 5: Defect rate P(chi != 1) ---")
    for n in range(5, 10):
        b1 = b1_rates.get(n, 0)
        b3 = b3_rates.get(n, 0)
        total = b1 + b3
        print(f"  n={n}: P(b1>0)+P(b3>0) = {100*b1:.2f}% + {100*b3:.2f}% = {100*total:.2f}%")

if __name__ == '__main__':
    main()
