#!/usr/bin/env python3
"""
Numerical verification of even-r-powers at n=7,8 with random tournaments.

For a c-tournament with t_ij = r + s_ij (s_ij = -s_ji), verify:
  M[a,b](r,s) = M[a,b](-r,s) for random s values and multiple r values.

This tests the even-r-powers conjecture efficiently for large n.

kind-pasteur-2026-03-06-S23b
"""
import random
from itertools import permutations
from collections import defaultdict

def compute_M(n, a, b, r, s):
    """Compute M[a,b] numerically for given r and skew-symmetric s."""
    U = [v for v in range(n) if v != a and v != b]

    def t(i, j):
        if i == j: return 0
        return r + s[i][j]

    total = 0.0
    for mask in range(1 << len(U)):
        S = [U[i] for i in range(len(U)) if mask & (1 << i)]
        R = [U[i] for i in range(len(U)) if not (mask & (1 << i))]
        sign = (-1)**len(S)

        S_set = list(sorted(set(S) | {a}))
        R_set = list(sorted(set(R) | {b}))

        # E_a: Ham paths through S_set ending at a
        ea = 0.0
        if len(S_set) == 1:
            ea = 1.0
        else:
            for p in permutations(S_set):
                if p[-1] != a: continue
                prod = 1.0
                for i in range(len(p)-1):
                    prod *= t(p[i], p[i+1])
                ea += prod

        # B_b: Ham paths through R_set starting at b
        bb = 0.0
        if len(R_set) == 1:
            bb = 1.0
        else:
            for p in permutations(R_set):
                if p[0] != b: continue
                prod = 1.0
                for i in range(len(p)-1):
                    prod *= t(p[i], p[i+1])
                bb += prod

        total += sign * ea * bb
    return total

def make_random_skew(n):
    s = [[0.0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            val = random.uniform(-2, 2)
            s[i][j] = val
            s[j][i] = -val
    return s

print("=" * 70)
print("NUMERIC VERIFICATION: EVEN r-POWERS for M[a,b]")
print("=" * 70)

for n in [5, 6, 7, 8]:
    print(f"\nn={n}:")
    num_tests = 50 if n <= 7 else 10
    max_diff = 0.0
    failures = 0

    for test in range(num_tests):
        s = make_random_skew(n)
        a, b = 0, 1
        r_val = random.uniform(0.5, 3.0)

        M_pos = compute_M(n, a, b, r_val, s)
        M_neg = compute_M(n, a, b, -r_val, s)

        diff = abs(M_pos - M_neg)
        max_diff = max(max_diff, diff)

        if diff > 1e-6:
            failures += 1
            print(f"  FAIL test {test}: M(r={r_val:.3f}) = {M_pos:.6f}, "
                  f"M(-r) = {M_neg:.6f}, diff = {diff:.2e}")

    if failures == 0:
        print(f"  ALL {num_tests} TESTS PASS. "
              f"Max |M(r)-M(-r)| = {max_diff:.2e}")
    else:
        print(f"  {failures}/{num_tests} FAILURES!")

    # Also check M[a,b] = M[b,a]
    sym_failures = 0
    max_sym_diff = 0.0
    for test in range(num_tests):
        s = make_random_skew(n)
        r_val = random.uniform(0.5, 3.0)
        a, b = random.sample(range(n), 2)

        M_ab = compute_M(n, a, b, r_val, s)
        M_ba = compute_M(n, b, a, r_val, s)

        diff = abs(M_ab - M_ba)
        max_sym_diff = max(max_sym_diff, diff)
        if diff > 1e-6:
            sym_failures += 1

    if sym_failures == 0:
        print(f"  SYMMETRY: ALL {num_tests} tests pass. "
              f"Max |M[a,b]-M[b,a]| = {max_sym_diff:.2e}")
    else:
        print(f"  SYMMETRY: {sym_failures}/{num_tests} FAILURES!")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
