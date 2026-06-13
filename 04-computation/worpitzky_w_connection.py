#!/usr/bin/env python3
"""
worpitzky_w_connection.py — Connecting Worpitzky expansion to W-polynomial.

THM-K gives: W(T,r) = (r-1/2)^{n-1} * F(T, (2r+1)/(2r-1))

The Worpitzky expansion gives: F(T,x)/(1-x)^n = sum a_m x^m
where a_m is polynomial in m with universal top coefficients.

QUESTION: What is the Worpitzky expansion in the W-world?

In the W-world, we expand: W(T,r) = sum_{k=0}^{n-1} w_k * r^k
where w_k = w_{n-1-2j} for appropriate j (odd coefficients vanish at odd n).

The W-coefficients form a hierarchy (INV-082):
- w_{n-1} = n! (universal)
- w_{n-3} = (n-2)! * [2*t3 - C(n,3)/2]
etc.

The Worpitzky coefficients also have this structure:
- c_{n-1} = n (universal)
- c_{n-2} = C(n,2)
- c_{n-3} = C(n,3) + 2(n-2)*t3

Is there a direct algebraic relationship between c_j and w_k?

Author: opus-2026-03-07-S46
"""
from itertools import permutations, combinations
from math import comb, factorial
import numpy as np

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

def compute_F(adj, n):
    F = [0]*n
    for P in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if adj[P[i]][P[i+1]])
        F[fwd] += 1
    return F

def count_3cycles(adj, n):
    t3 = 0
    for triple in combinations(range(n), 3):
        i, j, k = triple
        if (adj[i][j] and adj[j][k] and adj[k][i]) or \
           (adj[i][k] and adj[k][j] and adj[j][i]):
            t3 += 1
    return t3

def compute_W(F, n, r):
    """W(T,r) = (r-1/2)^{n-1} * F(T, (2r+1)/(2r-1))"""
    if abs(2*r - 1) < 1e-10:
        return None
    x = (2*r + 1) / (2*r - 1)
    Fx = sum(F[k] * x**k for k in range(n))
    return (r - 0.5)**(n-1) * Fx

def worpitzky_a(F, n, m):
    return sum(F[k] * comb(m + n - 1 - k, n - 1) for k in range(n))

# ============================================================
# COMPARE WORPITZKY AND W-COEFFICIENTS
# ============================================================
print("=" * 60)
print("WORPITZKY vs W-POLYNOMIAL COEFFICIENTS")
print("=" * 60)

for n in [5]:
    m_vals = n*(n-1)//2
    seen = set()

    for bits in range(1 << m_vals):
        adj = tournament_from_bits(n, bits)
        F = compute_F(adj, n)
        key = tuple(F)
        if key in seen:
            continue
        seen.add(key)

        t3 = count_3cycles(adj, n)

        # W-coefficients: W(T,r) as polynomial in r
        r_pts = [0.1 * i + 1 for i in range(n + 2)]
        w_pts = [compute_W(F, n, r) for r in r_pts]
        w_coeffs = np.polyfit(r_pts, w_pts, n - 1)  # high to low

        # Worpitzky coefficients
        m_pts = list(range(n + 2))
        a_pts = [worpitzky_a(F, n, m) for m in m_pts]
        c_coeffs = np.polyfit(m_pts, a_pts, n - 1)  # high to low

        print(f"\nF={F}, t3={t3}")
        print(f"  W-coeffs (r^4..r^0): {[f'{c:.1f}' for c in w_coeffs]}")
        print(f"  Worpitzky (m^4..m^0): {[f'{c:.1f}' for c in c_coeffs]}")

        # Ratios
        ratios = []
        for i in range(n):
            if abs(c_coeffs[i]) > 0.1 and abs(w_coeffs[i]) > 0.1:
                ratios.append(f"{w_coeffs[i]/c_coeffs[i]:.4f}")
            else:
                ratios.append("N/A")
        print(f"  W/Worpitzky ratios: {ratios}")

# ============================================================
# THE KEY TRANSFORM: F(x) ↔ W(r) ↔ a_m
# ============================================================
print("\n" + "=" * 60)
print("TRIANGLE OF TRANSFORMS")
print("=" * 60)

# We have THREE representations of the same tournament polynomial:
# 1. F(T,x) = sum F_k x^k (forward-edge polynomial)
# 2. W(T,r) = sum w_j r^j (W-polynomial)
# 3. a_m(T) = sum c_j m^j (Worpitzky polynomial)
#
# Relations:
# F ↔ W via THM-K: W(r) = (r-1/2)^{n-1} F((2r+1)/(2r-1))
# F ↔ a via Worpitzky: a_m = sum_k F_k C(m+n-1-k, n-1)
# W ↔ a: indirect via F
#
# But is there a DIRECT relation between W and a?
# Let's find out.

# Consider: what value of r corresponds to m in the Worpitzky sense?
# a_m = [coefficient of x^m in F(x)/(1-x)^n]
# W(r) = (r-1/2)^{n-1} * F((2r+1)/(2r-1))
# If x = (2r+1)/(2r-1), then r = (x+1)/(2(x-1)).
# At x → 1: r → ∞ (pole)
# At x = 0: r = -1/2
# At x = -1: r = 0
# At x = 2: r = 3/2

# The generating function F(x)/(1-x)^n diverges at x=1, which is where W is easiest.
# So the Worpitzky expansion works near x=0, while W works near x=1 (r=∞).
# They're complementary viewpoints!

print("\nWorpitzky works near x=0 (r=-1/2)")
print("W-polynomial works near r=∞ (x=1)")
print("These are DUAL expansions of the same function.")

# Verify: W at large r should give the leading Worpitzky structure
# At large r: x = (2r+1)/(2r-1) ≈ 1 + 1/r
# F(T, 1+1/r) ≈ F(T,1) + F'(T,1)/r + ...
# W(r) = (r-1/2)^{n-1} F(1+1/r) ≈ r^{n-1} * n! * [1 + corrections/r + ...]

# So the LARGE-r expansion of W gives the SAME hierarchy as Worpitzky!
# w_{n-1} = n! is analogous to c_{n-1} = n (both are universal top coefficients)
# The scaling differs by (n-1)! since W uses r^j while Worpitzky uses m^j/(j!)

print("\n\nScaling relationship: w_j ≈ (n-1)! * c_j * something?")
n = 5
for bits in [0]:
    adj = tournament_from_bits(n, bits)
    F = compute_F(adj, n)

    r_pts = np.linspace(1, 10, 20)
    w_pts = [compute_W(F, n, r) for r in r_pts]
    w_coeffs = np.polyfit(r_pts, w_pts, n-1)

    m_pts = list(range(n+2))
    a_pts = [worpitzky_a(F, n, m) for m in m_pts]
    c_coeffs = np.polyfit(m_pts, a_pts, n-1)

    print(f"\nTransitive T5: F={F}")
    print(f"  W-coeffs: {[f'{c:.1f}' for c in w_coeffs]}")
    print(f"  C-coeffs: {[f'{c:.1f}' for c in c_coeffs]}")
    for j in range(n):
        if abs(c_coeffs[j]) > 0.01:
            ratio = w_coeffs[j] / c_coeffs[j]
            print(f"  w/c at position {j} (m^{n-1-j}): ratio = {ratio:.4f}")

# ============================================================
# CRITICAL: Worpitzky + Variance = t3 formula
# ============================================================
print("\n" + "=" * 60)
print("VARIANCE OF fwd(P) AND t3")
print("=" * 60)

# Var[fwd] = E[fwd^2] - (E[fwd])^2
# E[fwd] = (n-1)/2 always
# E[fwd^2] = sum_k k^2 F_k / n!
# Var = E[fwd^2] - ((n-1)/2)^2
#
# The Worpitzky third coefficient c_{n-3} depends on E[fwd^2].
# c_{n-3} = C(n,3) + 2(n-2)*t3
# So 2(n-2)*t3 = c_{n-3} - C(n,3) is determined by E[fwd^2].
# This means: Var[fwd] = f(t3, n) for some function f.

for n in [4, 5, 6]:
    m_vals = n*(n-1)//2
    seen = set()

    import random
    random.seed(42)
    num = min(1 << m_vals, 50000)

    print(f"\nn={n}:")
    for trial in range(num):
        if n <= 5:
            bits = trial
        else:
            bits = random.getrandbits(m_vals)

        adj = tournament_from_bits(n, bits)
        F = compute_F(adj, n)
        key = tuple(F)
        if key in seen:
            continue
        seen.add(key)

        t3 = count_3cycles(adj, n)
        total = sum(F)
        mean = sum(k * F[k] for k in range(n)) / total
        efwd2 = sum(k**2 * F[k] for k in range(n)) / total
        var = efwd2 - mean**2

        if len(seen) <= 10:
            # Check formula: Var = A + B*t3 for some A, B
            print(f"  F={F}, t3={t3}, E[fwd]={(n-1)/2:.1f}, Var={var:.4f}")

    # Compute A, B from all data
    all_data = []
    seen2 = set()
    for trial in range(num):
        if n <= 5:
            bits = trial
        else:
            bits = random.getrandbits(m_vals)
        adj = tournament_from_bits(n, bits)
        F = compute_F(adj, n)
        key = tuple(F)
        if key in seen2:
            continue
        seen2.add(key)
        t3 = count_3cycles(adj, n)
        total = sum(F)
        var = sum(k**2 * F[k] for k in range(n)) / total - ((n-1)/2)**2
        all_data.append((t3, var))

    # Linear fit
    X = np.array([[1, d[0]] for d in all_data])
    y = np.array([d[1] for d in all_data])
    coeffs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    max_err = max(abs(X @ coeffs - y))
    print(f"  Var[fwd] = {coeffs[0]:.6f} + {coeffs[1]:.6f} * t3, max_err = {max_err:.8f}")

    # Theoretical prediction:
    # E[fwd] = (n-1)/2
    # E[X_i] = 1/2 for each of (n-1) steps
    # Var[sum X_i] = sum Var(X_i) + 2 sum_{i<j} Cov(X_i, X_j)
    # Var(X_i) = 1/4, sum = (n-1)/4
    # Cov(X_i, X_{i+1}) depends on 3-paths through shared vertex
    # Other Cov(X_i, X_j) for |i-j|>=2: all are the same? No, they depend on the specific pair.
    print(f"  Predicted base variance (independent): {(n-1)/4:.4f}")
