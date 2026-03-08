#!/usr/bin/env python3
"""
UNIVERSAL WORPITZKY COEFFICIENTS

THEOREM (just discovered):
W(T,r) = (H/2^{n-1}) * sum_k c_k C(r+k, n-1)

where c_k depends ONLY on n, not on T!

c_k(n) values:
  n=3: [1, 6, 1]
  n=4: [1, 23, 23, 1]
  n=5: [1, 76, 230, 76, 1]

These are palindromic: c_k = c_{n-1-k}

QUESTION: What are these numbers?
Let's check OEIS and known sequences.

c_1(3) = 6
c_1(4) = 23
c_1(5) = 76
c_2(5) = 230

Differences: 23-6=17, 76-23=53. Not obvious.
Ratios: 23/6 ≈ 3.83, 76/23 ≈ 3.30. Not clean.

Try: c_1(n) = ? Let me compute more values.
"""
from itertools import permutations
from math import comb, factorial
import random
import numpy as np

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

def W_poly_coeffs_dp(A, n):
    """DP computation of W(T,r) coefficients, O(n^2 * 2^n)"""
    # dp[(mask, v)] = polynomial coefficients of sum over paths ending at v visiting mask
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = np.array([1.0])  # constant polynomial = 1

    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or (mask, v) not in dp:
                continue
            pv = dp[(mask, v)]
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v][u]:
                    s_e = A[v][u] - 0.5  # +0.5 for forward
                    # Multiply current polynomial by (r + s_e)
                    new_poly = np.zeros(len(pv) + 1)
                    new_poly[:len(pv)] += s_e * pv  # s_e * p(r)
                    new_poly[1:len(pv)+1] += pv     # r * p(r)

                    key = (mask | (1 << u), u)
                    if key in dp:
                        # Extend if needed
                        if len(dp[key]) < len(new_poly):
                            old = dp[key]
                            dp[key] = np.zeros(len(new_poly))
                            dp[key][:len(old)] = old
                        elif len(new_poly) < len(dp[key]):
                            tmp = np.zeros(len(dp[key]))
                            tmp[:len(new_poly)] = new_poly
                            new_poly = tmp
                        dp[key] += new_poly
                    else:
                        dp[key] = new_poly.copy()

    full = (1 << n) - 1
    W = np.zeros(n)
    for v in range(n):
        if (full, v) in dp:
            p = dp[(full, v)]
            for k in range(min(n, len(p))):
                W[k] += p[k]
    return W

def ham_path_count_dp(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or (mask, v) not in dp: continue
            for u in range(n):
                if mask & (1 << u): continue
                if A[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + dp[(mask, v)]
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def worpitzky_basis_matrix(d):
    n = d + 1
    M = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            M[i][j] = comb(i + j, d)
    return M

def expand_in_worpitzky(poly_coeffs):
    d = len(poly_coeffs) - 1
    n = d + 1
    vals = [sum(poly_coeffs[k] * x**k for k in range(n)) for x in range(n)]
    M = worpitzky_basis_matrix(d)
    w = np.linalg.solve(M, vals)
    return w

random.seed(42)

# ===== COMPUTE UNIVERSAL CONSTANTS AT n=6,7,8 =====
print("=" * 70)
print("COMPUTING UNIVERSAL CONSTANTS c_k(n)")
print("=" * 70)
print("W(T,r) = (H/2^{n-1}) * sum_k c_k(n) * C(r+k, n-1)")
print()

all_constants = {}

for n in range(3, 9):
    print(f"\nn={n}:")
    scale = 2**(n-1)

    # Use a few random tournaments to compute c_k
    constants = None
    num_tests = 20 if n <= 7 else 10

    for trial in range(num_tests):
        A = random_tournament(n)
        H = ham_path_count_dp(A, n)
        W = W_poly_coeffs_dp(A, n)
        w = expand_in_worpitzky(W)

        c = [round(w[k] * scale / H, 6) for k in range(n)]

        if constants is None:
            constants = c
        else:
            # Verify
            if any(abs(c[k] - constants[k]) > 0.01 for k in range(n)):
                print(f"  VIOLATION at trial {trial}! c = {c} vs {constants}")
                break

    all_constants[n] = [int(round(x)) for x in constants]
    print(f"  c_k({n}) = {all_constants[n]}")
    print(f"  sum = {sum(all_constants[n])}")
    print(f"  palindromic: {all(all_constants[n][k] == all_constants[n][n-1-k] for k in range(n))}")

# ===== IDENTIFY THE SEQUENCE =====
print("\n\n" + "=" * 70)
print("IDENTIFYING THE SEQUENCE")
print("=" * 70)

# c_0(n) = 1 for all n (verified)
# c_1 sequence: 6, 23, 76, ?, ?
# c_k(n): full triangles

print("\nc_0(n): ", [all_constants[n][0] for n in sorted(all_constants.keys())])
print("c_1(n): ", [all_constants[n][1] for n in sorted(all_constants.keys()) if n >= 4])

for k in range(5):
    vals = []
    for n in sorted(all_constants.keys()):
        if k < n:
            vals.append(all_constants[n][k])
    print(f"c_{k}: {vals}")

# Check: is c_k(n) = A(n-1,k) * f(n) for some function f?
print("\n\nRatio c_k(n) / A(n-1,k) (Eulerian number):")
def eulerian_number(n, k):
    return sum((-1)**j * comb(n+1, j) * (k+1-j)**n for j in range(k+2))

for n in sorted(all_constants.keys()):
    d = n - 1
    ratios = []
    for k in range(n):
        E = eulerian_number(d, k)
        c = all_constants[n][k]
        if E != 0:
            ratios.append(c / E)
        else:
            ratios.append(None)
    print(f"  n={n}: {[f'{r:.4f}' if r is not None else 'N/A' for r in ratios]}")

# Check: sum c_k(n) = ?
print("\n\nsum c_k(n):")
for n in sorted(all_constants.keys()):
    s = sum(all_constants[n])
    print(f"  n={n}: sum = {s}")
    # Check if sum = (2n-2)!! or similar
    print(f"    (2n-3)!! = {int(np.prod(range(1, 2*n-2, 2)))}")
    print(f"    (n-1)! * 2^{n-1} = {factorial(n-1) * 2**(n-1)}")

# Try to identify via generating function
print("\n\nFull triangle:")
for n in sorted(all_constants.keys()):
    print(f"  n={n}: {all_constants[n]}")

# Check OEIS sequences
print("\n\nLooking for c_1 sequence: ", [all_constants[n][1] for n in sorted(all_constants.keys()) if n >= 3])
print("Checking specific formulas:")
for n in sorted(all_constants.keys()):
    c1 = all_constants[n][1] if n >= 3 else None
    # Try: c_1 = (2^{n-1} - 1) * n / gcd?
    # n=3: c_1=6.  2^2-1=3. 3*2=6. YES!
    # n=4: c_1=23. 2^3-1=7. 7*?=23? No.
    # Try: c_1 = 3^{n-1} - 1
    # n=3: 9-1=8. No.
    # Try: c_1 = sum_{j=0}^{n-2} C(n-1, j) * (j+1)^{n-2} ... complicated
    if c1 is not None:
        print(f"  n={n}: c_1 = {c1}")
        # Various guesses
        for formula_name, formula_val in [
            ("3^{n-1} - 2^{n-1}", 3**(n-1) - 2**(n-1)),
            ("3^{n-1}", 3**(n-1)),
            ("2^{2n-3} - 1", 2**(2*n-3) - 1),
            ("(2n-2)!/(n-1)!", factorial(2*n-2)//factorial(n-1)),
            ("C(2n-2, n-1)", comb(2*n-2, n-1)),
        ]:
            match = "YES" if formula_val == c1 else ""
            if match:
                print(f"    {formula_name} = {formula_val} {match}")
