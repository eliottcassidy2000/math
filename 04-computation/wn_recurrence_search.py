#!/usr/bin/env python3
"""Search for recurrence relations for W(n) with polynomial coefficients.
Try W(n) = P(n)*W(n-1) + Q(n)*W(n-2) + R(n)*W(n-3) + ...
where P, Q, R are polynomials in n of varying degree.
"""
from fractions import Fraction as F
import itertools

wn = {1: 1, 2: 2, 3: 8, 4: 32, 5: 158, 6: 928, 7: 6350, 8: 49752, 9: 439670,
      10: 4327904, 11: 46963358, 12: 556953448, 13: 7166360054,
      14: 99428495088, 15: 1479600188798, 16: 23506712352248,
      17: 397095175477430, 18: 7107209383674112, 19: 134345623603516190,
      20: 2674426516381764744, 21: 55925062706620208438,
      22: 1225582324497129993488, 23: 28088043261491650347134,
      24: 671901551280362054066968, 25: 16746599265666151628174198,
      26: 434187457363955400414175008, 27: 11692423738081050318010736030}

def try_recurrence(depth, poly_degrees):
    """Try recurrence W(n) = sum_{i=1}^{depth} P_i(n)*W(n-i)
    where P_i has degree poly_degrees[i-1].
    Returns coefficients if found, None otherwise.
    """
    # Total unknowns
    num_unknowns = sum(d + 1 for d in poly_degrees)

    # Need num_unknowns equations, starting from n = depth+1
    start_n = depth + 1
    end_n = start_n + num_unknowns

    if end_n > 27:
        return None  # not enough data

    # Build system: for each n, sum of P_i(n)*W(n-i) = W(n)
    # Unknowns are coefficients of each P_i
    rows = []
    for n in range(start_n, end_n):
        row = []
        for i in range(depth):
            deg = poly_degrees[i]
            for p in range(deg + 1):  # n^p * W(n-i-1)
                row.append(F(n**p * wn[n - i - 1]))
            # row coefficients for P_i: c_0*W(n-i-1) + c_1*n*W(n-i-1) + ... + c_d*n^d*W(n-i-1)
        rhs = F(wn[n])
        rows.append((row, rhs))

    # Gaussian elimination
    mat = [list(row) + [rhs] for row, rhs in rows]
    m = num_unknowns
    for col in range(m):
        pivot = None
        for r in range(col, m):
            if mat[r][col] != 0:
                pivot = r
                break
        if pivot is None:
            return None
        mat[col], mat[pivot] = mat[pivot], mat[col]
        for r in range(m):
            if r != col and mat[r][col] != 0:
                factor = F(mat[r][col], mat[col][col])
                for c in range(m + 1):
                    mat[r][c] -= factor * mat[col][c]

    sol = [F(mat[i][m], mat[i][i]) for i in range(m)]

    # Verify on remaining data points
    verify_start = end_n
    verify_end = min(28, verify_start + 5)

    for n in range(verify_start, verify_end):
        if n not in wn:
            break
        predicted = F(0)
        idx = 0
        for i in range(depth):
            deg = poly_degrees[i]
            poly_val = F(0)
            for p in range(deg + 1):
                poly_val += sol[idx] * F(n**p)
                idx += 1
            predicted += poly_val * F(wn[n - i - 1])

        if predicted != wn[n]:
            return None

    # Extract polynomials
    polys = []
    idx = 0
    for i in range(depth):
        deg = poly_degrees[i]
        coeffs = [sol[idx + p] for p in range(deg + 1)]
        polys.append(coeffs)
        idx += deg + 1

    return polys

# Systematic search
print("=== Searching for W(n) recurrence ===\n")

# 2-term recurrences: W(n) = P(n)*W(n-1) + Q(n)*W(n-2)
for dp in range(1, 6):
    for dq in range(0, dp + 1):
        result = try_recurrence(2, [dp, dq])
        if result is not None:
            print(f"FOUND: depth=2, deg(P)={dp}, deg(Q)={dq}")
            for i, coeffs in enumerate(result):
                print(f"  P_{i+1}(n) = {' + '.join(f'{float(c):.6f}*n^{p}' for p, c in enumerate(coeffs))}")
            print()

# 3-term recurrences
for dp in range(1, 4):
    for dq in range(0, dp + 1):
        for dr in range(0, dq + 1):
            result = try_recurrence(3, [dp, dq, dr])
            if result is not None:
                print(f"FOUND: depth=3, deg(P)={dp}, deg(Q)={dq}, deg(R)={dr}")
                for i, coeffs in enumerate(result):
                    print(f"  P_{i+1}(n) = {' + '.join(f'{float(c):.6f}*n^{p}' for p, c in enumerate(coeffs))}")
                print()

# 4-term with low-degree polys
for dp in range(1, 3):
    for dq in range(0, dp + 1):
        for dr in range(0, dq + 1):
            for ds in range(0, dr + 1):
                result = try_recurrence(4, [dp, dq, dr, ds])
                if result is not None:
                    print(f"FOUND: depth=4, deg(P)={dp}, deg(Q)={dq}, deg(R)={dr}, deg(S)={ds}")
                    for i, coeffs in enumerate(result):
                        print(f"  P_{i+1}(n) = {' + '.join(f'{float(c):.6f}*n^{p}' for p, c in enumerate(coeffs))}")
                    print()

# Also try NUD(n) recurrence analog. NUD(n) = (n-1)*NUD(n-1) + (n-2)*NUD(n-2)
# Maybe W(n) = (n + a)*W(n-1) + (n + b)*W(n-2) for some a, b?
print("=== Direct check: W(n) = (n+a)*W(n-1) + (n+b)*W(n-2) ===")
# Two unknowns, use n=3,4
# W(3) = (3+a)*W(2) + (3+b)*W(1) => 8 = 2(3+a) + (3+b) => 8 = 6+2a+3+b => 2a+b = -1
# W(4) = (4+a)*W(3) + (4+b)*W(2) => 32 = 8(4+a) + 2(4+b) => 32 = 32+8a+8+2b => 8a+2b = -8 => 4a+b=-4
# From these: 2a = -3, a = -3/2, b = -1+3 = 2
a_coeff = F(-3, 2)
b_coeff = F(2)
print(f"  a = {a_coeff}, b = {b_coeff}")
print("  Checking:")
for n in range(3, 28):
    pred = (n + a_coeff) * wn[n-1] + (n + b_coeff) * wn[n-2]
    match = pred == wn[n]
    if not match:
        print(f"  n={n}: FAIL (pred={pred}, actual={wn[n]}, diff={pred - wn[n]})")
        break
    else:
        print(f"  n={n}: OK")

# What about W(n) = (n-1)*W(n-1) + (n-2)*W(n-2) + correction?
print("\n=== Correction: W(n) - (n-1)*W(n-1) - (n-2)*W(n-2) ===")
for n in range(3, 28):
    corr = wn[n] - (n-1)*wn[n-1] - (n-2)*wn[n-2]
    print(f"  n={n}: {corr}")
