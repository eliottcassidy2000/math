"""
opus-2026-05-22-S4: Fast 3-cycle independence polynomial for all-0 staircase T_k.

Key structure: ALL 3-cycles in T_k come from exactly 2 pairs.
For pairs a < b (0-indexed), two 3-cycles exist:
  A(a,b): (2a+1) -> (2a) -> (2b) -> (2a+1)
  B(a,b): (2b+1) -> (2b) -> (2a+1) -> (2b+1)

This gives 2*C(k,2) = k(k-1) total 3-cycles.
Conflict graph Omega_3 built directly from vertex overlap.
"""
import sys
from fractions import Fraction
from itertools import combinations
import time

def compute_ip_3cycle(k, verbose=False):
    """Compute 3-cycle independence polynomial for T_k."""
    t0 = time.time()
    # Label cycles: A(a,b) = 2*(a*k+b) — no, simpler:
    # Index: A(a,b) -> idx 2*pair_idx, B(a,b) -> idx 2*pair_idx+1
    # pair_idx = a*(2k-a-1)//2 + (b-a-1) for 0<=a<b<k (C(k,2) pairs)

    # Build cycle list: (type, a, b, vertices)
    cycles = []
    for a in range(k):
        for b in range(a+1, k):
            # A(a,b): vertices {2a, 2a+1, 2b}
            cycles.append(('A', a, b, frozenset([2*a, 2*a+1, 2*b])))
            # B(a,b): vertices {2a+1, 2b, 2b+1}
            cycles.append(('B', a, b, frozenset([2*a+1, 2*b, 2*b+1])))

    m = len(cycles)  # = k*(k-1)

    # Build adjacency: two cycles conflict iff they share a vertex
    adj = [set() for _ in range(m)]
    for i in range(m):
        for j in range(i+1, m):
            if cycles[i][3] & cycles[j][3]:  # share a vertex
                adj[i].add(j)
                adj[j].add(i)

    # Compute independence polynomial via backtracking
    # ip[s] = number of independent sets of size s
    ip = [0] * (m + 1)
    ip[0] = 1

    def backtrack(start, size, forbidden):
        if size > 0:
            ip[size] += 1
        for v in range(start, m):
            if v not in forbidden:
                backtrack(v+1, size+1, forbidden | adj[v])

    # Better: use DP on vertices in order
    # dp[mask] is expensive; use recursive IS counting with pruning
    def count_is(idx, current_size, excluded):
        """Count all IS starting from idx, with excluded set."""
        ip[current_size] += 1
        for i in range(idx, m):
            if i not in excluded:
                count_is(i+1, current_size+1, excluded | adj[i])

    # Run the counting
    count_is(0, 0, set())

    elapsed = time.time() - t0
    # Trim trailing zeros
    while len(ip) > 1 and ip[-1] == 0:
        ip.pop()

    if verbose:
        print(f"k={k}: {ip} [{elapsed:.1f}s]")
    return ip, elapsed

def extract_coefficients(data, m, num_terms):
    """
    Given data = [(k, alpha_m(k))], express alpha_m as
    sum_{j=0}^{num_terms-1} c_j * C(k, 2m-j).
    Uses Gaussian elimination with fractions.
    """
    from fractions import Fraction

    def C(n, r):
        if r < 0 or r > n: return 0
        from math import comb
        return comb(n, r)

    # Build system: for each k in data, sum_j c_j * C(k, 2m-j) = alpha_m(k)
    n_eq = len(data)
    n_var = num_terms

    A = [[Fraction(C(k, 2*m - j)) for j in range(n_var)] for k, _ in data]
    b = [Fraction(val) for _, val in data]

    # Gaussian elimination
    for col in range(min(n_eq, n_var)):
        # Find pivot
        pivot = None
        for row in range(col, n_eq):
            if A[row][col] != 0:
                pivot = row
                break
        if pivot is None:
            continue
        A[col], A[pivot] = A[pivot], A[col]
        b[col], b[pivot] = b[pivot], b[col]

        for row in range(n_eq):
            if row != col and A[row][col] != 0:
                factor = A[row][col] / A[col][col]
                for c2 in range(n_var):
                    A[row][c2] -= factor * A[col][c2]
                b[row] -= factor * b[col]

    coeffs = [b[i] / A[i][i] if A[i][i] != 0 else Fraction(0) for i in range(min(n_eq, n_var))]
    return coeffs

if __name__ == '__main__':
    print("=== 3-CYCLE IP FOR ALL-0 STAIRCASE T_k ===\n")

    # Compute for k=2..18
    results = {}
    for k in range(2, 19):
        ip, t = compute_ip_3cycle(k)
        d = len(ip) - 1
        print(f"k={k}: d={d}, I_3={ip} [{t:.1f}s]")
        results[k] = ip

    print("\n=== EXTRACTING alpha_m FORMULAS ===")

    # For each m, collect data points where d >= m
    from math import floor, comb

    def num_terms(m):
        return m // 2 + 1

    for m in range(1, 13):
        data = [(k, results[k][m]) for k in results if len(results[k]) > m]
        if len(data) < num_terms(m) + 1:
            print(f"alpha_{m}: insufficient data (have {len(data)}, need {num_terms(m)+1})")
            continue

        nt = num_terms(m)
        coeffs = extract_coefficients(data[:nt+2], m, nt)

        terms = []
        for j, c in enumerate(coeffs):
            if c != 0:
                terms.append(f"({c})*C(k,{2*m-j})")
        formula = " + ".join(terms) if terms else "0"
        print(f"alpha_{m} = {formula}")

        # Verify on all available data
        errors = 0
        for k, val in data:
            pred = sum(c * comb(k, 2*m - j) for j, c in enumerate(coeffs))
            if pred != val:
                print(f"  ERROR at k={k}: predicted {pred}, actual {val}")
                errors += 1
        if errors == 0:
            print(f"  Verified on all {len(data)} data points ✓")

    print("\n=== LAST COEFFICIENTS (diagonal) ===")
    for m in range(1, 13):
        data = [(k, results[k][m]) for k in results if len(results[k]) > m]
        if data:
            k_min, val = min(data, key=lambda x: x[0])
            print(f"c_{{{m},{m//2}}} = alpha_{m}(k={k_min}) = {val}")
