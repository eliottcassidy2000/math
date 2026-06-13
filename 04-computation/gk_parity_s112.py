#!/usr/bin/env python3
"""
gk_parity_s112.py — The parity structure of TM g_k
kind-pasteur-2026-03-15-S112

DISCOVERY: D_k * g_k(m) is a polynomial with ONLY even/odd powers matching parity of k.
D_k appears to be related to double factorials.

Verify and extend to k=8,9,10.
"""

from fractions import Fraction
from math import factorial, gcd

def transfer_gk_values(k, m_max):
    results = {}
    for m in range(0, m_max + 1):
        n = m + 2*k
        num_edges = n - 2
        k_max = (n - 1) // 2
        if k > k_max:
            continue
        state = [[Fraction(1)] + [Fraction(0)] * k_max,
                 [Fraction(0)] * (k_max + 1),
                 [Fraction(0)] * (k_max + 1)]
        for step in range(num_edges):
            A, B, C = state
            nA = [A[i] + C[i] for i in range(k_max + 1)]
            nB = [Fraction(0)] + [2*A[i] + C[i] for i in range(k_max)]
            nC = list(B)
            state = [nA, nB, nC]
        total = [state[0][i] + state[1][i] + state[2][i] for i in range(k_max + 1)]
        if k < len(total) and total[k] != 0:
            results[m] = total[k] / 2
        else:
            results[m] = Fraction(0)
    return results

# Solve for polynomial of degree k for each k
for k in range(1, 11):
    tm = transfer_gk_values(k, k + 6)

    n_pts = k + 1
    ms = list(range(1, n_pts + 1))
    vs = [tm[m] for m in ms]

    # Vandermonde solve
    mat = [[Fraction(m**j) for j in range(n_pts)] for m in ms]
    rhs = list(vs)

    for col in range(n_pts):
        for row in range(col, n_pts):
            if mat[row][col] != 0:
                mat[col], mat[row] = mat[row], mat[col]
                rhs[col], rhs[row] = rhs[row], rhs[col]
                break
        pivot = mat[col][col]
        for j in range(n_pts):
            mat[col][j] /= pivot
        rhs[col] /= pivot
        for row in range(n_pts):
            if row != col and mat[row][col] != 0:
                factor = mat[row][col]
                for j in range(n_pts):
                    mat[row][j] -= factor * mat[col][j]
                rhs[row] -= factor * rhs[col]

    # Find common denominator
    denoms = [c.denominator for c in rhs if c != 0]
    lcm = 1
    for d in denoms:
        lcm = lcm * d // gcd(lcm, d)

    scaled = [int(c * lcm) for c in rhs]

    # Verify
    ok = True
    for m in range(1, min(k + 6, 15)):
        if m in tm:
            pred = sum(rhs[j] * Fraction(m**j) for j in range(n_pts))
            if pred != tm[m]:
                ok = False

    # Check parity
    even_powers = [scaled[j] for j in range(0, n_pts, 2)]
    odd_powers = [scaled[j] for j in range(1, n_pts, 2)]
    k_parity = k % 2  # 0 for even k, 1 for odd k

    if k_parity == 0:  # even k: only even powers should be nonzero
        wrong_parity = [scaled[j] for j in range(1, n_pts, 2) if scaled[j] != 0]
    else:  # odd k: only odd powers
        wrong_parity = [scaled[j] for j in range(0, n_pts, 2) if scaled[j] != 0]

    parity_ok = len(wrong_parity) == 0

    terms = []
    for j in range(n_pts-1, -1, -1):
        if scaled[j] != 0:
            terms.append(f"{scaled[j]}*m^{j}" if j > 0 else f"{scaled[j]}")

    double_fact = 1
    for i in range(1, 2*k, 2):
        double_fact *= i

    print(f"k={k:2d}: {lcm}*g_{k}(m) = {' + '.join(terms)}")
    print(f"       parity={'ODD' if k_parity else 'EVEN'} powers only: {parity_ok}")
    print(f"       D_k={lcm}, (2k-1)!!={double_fact}, D_k/(2k-1)!!={Fraction(lcm, double_fact)}")
    print(f"       verified={ok}")

# Check the sequence of denominators
print("\nDenominator sequence:")
denoms = []
for k in range(1, 11):
    tm = transfer_gk_values(k, k + 6)
    n_pts = k + 1
    ms = list(range(1, n_pts + 1))
    vs = [tm[m] for m in ms]
    mat = [[Fraction(m**j) for j in range(n_pts)] for m in ms]
    rhs = list(vs)
    for col in range(n_pts):
        for row in range(col, n_pts):
            if mat[row][col] != 0:
                mat[col], mat[row] = mat[row], mat[col]
                rhs[col], rhs[row] = rhs[row], rhs[col]
                break
        pivot = mat[col][col]
        for j in range(n_pts):
            mat[col][j] /= pivot
        rhs[col] /= pivot
        for row in range(n_pts):
            if row != col and mat[row][col] != 0:
                factor = mat[row][col]
                for j in range(n_pts):
                    mat[row][j] -= factor * mat[col][j]
                rhs[row] -= factor * rhs[col]
    ds = [c.denominator for c in rhs if c != 0]
    lcm = 1
    for d in ds:
        lcm = lcm * d // gcd(lcm, d)
    denoms.append(lcm)

print(f"D_k = {denoms}")
print(f"(2k-1)!! = {[1] + [int(factorial(2*k)/(factorial(k)*2**k)) for k in range(1, 10)]}")

# Check: are g_k(m) related to m * some function of m^2?
# For odd k: g_k(m) = m * h_k(m^2) for some polynomial h_k
# For even k: g_k(m) = h_k(m^2)
print("\nSubstitution u = m^2:")
for k in range(1, 8):
    tm = transfer_gk_values(k, k + 4)
    if k % 2 == 1:
        # g_k(m) = m * h(m^2). Compute h(u) = g_k(sqrt(u))/sqrt(u) at u=1,4,9,16,...
        u_vals = [(m, tm[m] / m) for m in range(1, k + 3)]
        print(f"k={k} (odd): g_{k}(m)/m at m=1..{k+2}: {[str(v) for _, v in u_vals]}")
    else:
        # g_k(m) = h(m^2). Compute h(u) at u=1,4,9,...
        u_vals = [(m**2, tm[m]) for m in range(1, k + 3)]
        print(f"k={k} (even): g_{k}(m) at m^2=1,4,9,...: {[str(v) for _, v in u_vals]}")

print("\nDone!")
