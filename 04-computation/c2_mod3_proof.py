"""
c2_mod3_proof.py
kind-pasteur-2026-03-07-S36

ALGEBRAIC PROOF that S_r = 0 mod 3 for all n >= 5.

Write F(T,x) = sum_k c_k (x-1)^k (Taylor expansion around x=1).
Over F_3: x^3 - 1 = (x-1)^3. So S_r = 0 mod 3 for all r iff (x-1)^3 | F(T,x) mod 3.
This requires c_0, c_1, c_2 all = 0 mod 3.

  c_0 = F(T,1) = n! (always = 0 mod 3 for n >= 3)

  c_1 = F'(T,1) = sum_P asc(P) = n!*(n-1)/2
    (tournament-INDEPENDENT! each position contributes n!/(n(n-1))*C(n,2) by symmetry)
    = 0 mod 3 for n >= 3 (factor of 3 from n!)

  c_2 = F''(T,1)/2 = sum_P C(asc(P), 2) = (sum_P asc^2 - sum_P asc) / 2
    Decompose into overlapping and non-overlapping position pairs:
    c_2 = A_non + (n-2)! * dp(T)
    where:
      A_non = C(n-1,2)-(n-2)) * (n-4)! * n(n-1)(n-2)(n-3)/4
            = [(n-2)(n-3)/2] * (n-4)! * [n(n-1)(n-2)(n-3)/4]
      dp(T) = sum_v outdeg(v)*indeg(v) = sum_v s_v*(n-1-s_v)  (2-path count)

    Key: A_non = 0 mod 3 for n >= 5 (product contains factor of 3)
         (n-2)! = 0 mod 3 for n >= 5 (since n-2 >= 3)

    Therefore c_2 = 0 mod 3 for ALL T at n >= 5, INDEPENDENT OF dp(T).

COMBINED: c_0, c_1, c_2 all = 0 mod 3 for n >= 5
=> (x-1)^3 | F(T,x) mod 3
=> S_r = 0 mod 3 for all r
=> F(T,omega) = 0 mod 9 when additionally v_3(n!) >= 2, i.e., n >= 6.

This script verifies the formula and the theorem.
"""

import os; os.environ['PYTHONIOENCODING'] = 'utf-8'
from collections import defaultdict
import random


def tournament_from_bits(n, bits):
    adj = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            idx += 1
    return adj


def compute_F_dp(adj, n):
    """Compute F(T,x) coefficients using DP."""
    dp = [[[0] * n for _ in range(n)] for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v][0] = 1

    for mask in range(1, 1 << n):
        for last in range(n):
            if not (mask & (1 << last)):
                continue
            for fwd in range(n):
                if dp[mask][last][fwd] == 0:
                    continue
                for nxt in range(n):
                    if mask & (1 << nxt):
                        continue
                    new_mask = mask | (1 << nxt)
                    if adj[last][nxt]:
                        dp[new_mask][nxt][fwd + 1] += dp[mask][last][fwd]
                    else:
                        dp[new_mask][nxt][fwd] += dp[mask][last][fwd]

    full = (1 << n) - 1
    F = [0] * n
    for last in range(n):
        for fwd in range(n):
            F[fwd] += dp[full][last][fwd]
    return F


def factorial(n):
    r = 1
    for i in range(1, n+1):
        r *= i
    return r


def v3(n):
    if n == 0:
        return float('inf')
    v = 0
    n = abs(n)
    while n % 3 == 0:
        v += 1
        n //= 3
    return v


def dp_count(adj, n):
    """Count directed 2-paths: sum_v outdeg(v)*indeg(v)."""
    total = 0
    for v in range(n):
        out_v = sum(adj[v])
        in_v = sum(adj[u][v] for u in range(n))
        total += out_v * in_v
    return total


def compute_c012(F, n):
    """Compute Taylor coefficients c_0, c_1, c_2 of F(T,x) around x=1."""
    c0 = sum(F)  # F(T,1) = n!
    c1 = sum(k * F[k] for k in range(n))  # F'(T,1)
    c2 = sum(k * (k-1) * F[k] for k in range(n)) // 2  # F''(T,1)/2
    return c0, c1, c2


def A_non_formula(n):
    """Compute A_non = [(n-2)(n-3)/2] * (n-4)! * [n(n-1)(n-2)(n-3)/4]."""
    if n < 4:
        return 0
    # Number of non-overlapping position pairs
    num_pairs = (n-2)*(n-3)//2
    # Number of vertex-disjoint arc pairs (constant for all tournaments)
    arc_pairs = n*(n-1)*(n-2)*(n-3)//4
    # Fill remaining positions
    fill = factorial(n-4) if n >= 4 else 1
    return num_pairs * fill * arc_pairs


# Verify the decomposition c_2 = A_non + (n-2)! * dp(T)
print("=" * 60)
print("Verification: c_2(T) = A_non + (n-2)! * dp(T)")
print("=" * 60)

for n in [4, 5, 6, 7]:
    m = n * (n - 1) // 2
    num_T = 1 << m
    a_non = A_non_formula(n)
    fac_nm2 = factorial(n - 2)
    errors = 0

    # Test all (n<=6) or sample (n=7)
    if n <= 6:
        indices = range(num_T)
    else:
        random.seed(42)
        indices = [random.randint(0, num_T-1) for _ in range(2000)]

    for bits in indices:
        adj = tournament_from_bits(n, bits)
        F = compute_F_dp(adj, n)
        c0, c1, c2_actual = compute_c012(F, n)
        dp = dp_count(adj, n)
        c2_formula = a_non + fac_nm2 * dp
        if c2_actual != c2_formula:
            errors += 1
            if errors <= 3:
                print(f"  MISMATCH at n={n}, bits={bits}: c2={c2_actual}, formula={c2_formula}")

    tested = num_T if n <= 6 else len(indices)
    print(f"  n={n}: {tested-errors}/{tested} match  (A_non={a_non}, (n-2)!={fac_nm2})")
    print(f"    v_3(A_non) = {v3(a_non) if a_non > 0 else 'inf'}, v_3((n-2)!) = {v3(fac_nm2)}")


# Print the proof summary
print("\n" + "=" * 60)
print("PROOF SUMMARY")
print("=" * 60)

for n in range(3, 11):
    c0 = factorial(n)
    c1 = factorial(n) * (n-1) // 2
    a_non = A_non_formula(n) if n >= 4 else 0
    fac_nm2 = factorial(n-2) if n >= 2 else 1
    v3_c0 = v3(c0)
    v3_c1 = v3(c1)
    v3_anon = v3(a_non) if a_non > 0 else 'inf'
    v3_facnm2 = v3(fac_nm2) if fac_nm2 > 0 else 0

    c2_always_0_mod3 = (v3_anon != 'inf' and v3_anon >= 1 and v3_facnm2 >= 1) or \
                       (v3_anon == 'inf' and v3_facnm2 >= 1) or \
                       (isinstance(v3_anon, int) and v3_anon >= 1 and v3_facnm2 >= 1)

    sr_all_0 = v3_c0 >= 1 and v3_c1 >= 1 and c2_always_0_mod3
    mod9 = sr_all_0 and v3_c0 >= 2

    print(f"  n={n:2d}: c0={c0:>8d} (v3={v3_c0}), c1={c1:>8d} (v3={v3_c1}), "
          f"A_non v3={str(v3_anon):>3s}, (n-2)! v3={v3_facnm2}")
    print(f"         3|c_2 for ALL T: {'YES' if c2_always_0_mod3 else 'NO'}  "
          f"S_r=0 mod 3: {'YES' if sr_all_0 else 'NO'}  "
          f"F(omega)=0 mod 9: {'YES' if mod9 else 'NO'}")


# Final exhaustive verification at n=5,6
print("\n" + "=" * 60)
print("Exhaustive verification")
print("=" * 60)

for n in [5, 6]:
    m = n * (n-1) // 2
    num_T = 1 << m
    c2_failures = 0
    sr_failures = 0
    mod9_failures = 0

    for bits in range(num_T):
        adj = tournament_from_bits(n, bits)
        F = compute_F_dp(adj, n)
        c0, c1, c2 = compute_c012(F, n)

        if c2 % 3 != 0:
            c2_failures += 1

        S = [0, 0, 0]
        for k in range(n):
            S[k % 3] += F[k]
        if any(s % 3 != 0 for s in S):
            sr_failures += 1

        # Check mod 9
        a = S[0] - S[2]
        b = S[1] - S[2]
        if a % 9 != 0 or b % 9 != 0:
            mod9_failures += 1

    print(f"  n={n} ({num_T} tournaments):")
    print(f"    c_2 = 0 mod 3: {num_T-c2_failures}/{num_T}")
    print(f"    S_r = 0 mod 3: {num_T-sr_failures}/{num_T}")
    print(f"    F(omega) = 0 mod 9: {num_T-mod9_failures}/{num_T}")


print("\nDONE")
