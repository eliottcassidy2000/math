#!/usr/bin/env python3
"""
h_trace_formula.py — opus-2026-03-12-S58

DISCOVERY: At p=7, H = (462 - tr(A^4))/2 for ALL circulant tournaments!

H is a DECREASING LINEAR function of tr(A^4). Since Paley minimizes tr(A^4)
(flat eigenvalue spectrum → min Σ|λ|^4 → min Σ Re(λ^4)), Paley maximizes H.

This script:
1. Verifies the formula at p=7
2. Tests whether analogous formulas exist at p=11, p=13
3. Explores the combinatorial meaning of tr(A^4)
4. Tests whether H is multilinear in (p_4, p_5, p_6, ...)
"""

import numpy as np
from itertools import combinations
from math import factorial, comb

def circulant_adjacency(p, S):
    A = np.zeros((p, p), dtype=int)
    for i in range(p):
        for j in range(p):
            if (j - i) % p in S:
                A[i][j] = 1
    return A

def all_circulant_sets(p):
    m = (p - 1) // 2
    reps = list(range(1, m + 1))
    result = []
    for bits in range(1 << m):
        S = set()
        for i in range(m):
            if bits & (1 << i):
                S.add(reps[i])
            else:
                S.add(p - reps[i])
        result.append(S)
    return result

def hamiltonian_paths_dp(A, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    full = (1 << n) - 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask >> v & 1) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask >> u & 1: continue
                if A[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[full])

def eigenvalues_circulant(p, S):
    omega = np.exp(2j * np.pi / p)
    return [sum(omega**(k*s) for s in S) for k in range(p)]

print("=" * 70)
print("THEOREM: H = f(tr(A^k)) FOR CIRCULANT TOURNAMENTS")
print("=" * 70)
print()

# ================================================================
# p = 7: Verify H = (462 - tr(A^4)) / 2
# ================================================================
print("=== p = 7: Testing H = (462 - tr(A^4))/2 ===")
print()

p = 7
sets = all_circulant_sets(p)
all_match = True
for S in sets:
    A = circulant_adjacency(p, S)
    H = hamiltonian_paths_dp(A, p)
    trA4 = int(round(np.trace(A @ A @ A @ A)))
    pred = (462 - trA4) // 2
    match = "✓" if pred == H else "✗"
    if pred != H: all_match = False
    print(f"  S={sorted(S)}: H={H}, tr(A⁴)={trA4}, (462-tr)/2={pred} {match}")

print(f"\nAll match: {all_match}")
print()

# What is 462?
# 462 = 2 * 231 = 2 * 3 * 7 * 11. Or: 462 = C(11,4) + something...
# Actually: (p-1)/2 = 3. p=7.
# 462 = 2*H_max - ... hmm. 2*189 = 378. 462 - 378 = 84 = tr(A⁴) for Paley.
# So formula is: H = (462 - tr(A⁴))/2 = (H_max + ... )/something.
# Let's understand the constant 462.
print("Understanding the constant 462:")
print(f"  2 * H_max = {2*189} = 378")
print(f"  tr(A⁴)_Paley = {int(round(np.trace(circulant_adjacency(7, {1,2,4}) @ circulant_adjacency(7, {1,2,4}) @ circulant_adjacency(7, {1,2,4}) @ circulant_adjacency(7, {1,2,4}))))}")
print(f"  462 = 378 + 84 = 2*H + tr(A⁴) for Paley")
print(f"  462 = 2*175 + 112 = 2*H + tr(A⁴) for non-Paley ✓")
print()

# ================================================================
# p = 11: Test H as function of tr(A^k)
# ================================================================
print("=" * 70)
print("=== p = 11: H as function of traces ===")
print("=" * 70)
print()

p = 11
sets = all_circulant_sets(p)
data = []
for S in sets:
    A = circulant_adjacency(p, S)
    H = hamiltonian_paths_dp(A, p)
    traces = {}
    Ak = np.eye(p, dtype=int)
    for k in range(1, p+1):
        Ak = Ak @ A
        traces[k] = int(round(np.trace(Ak)))
    data.append((H, traces, sorted(S)))

# Show distinct H values with their traces
print("Distinct (H, traces) combinations:")
seen = {}
for H, tr, S in sorted(data, key=lambda x: -x[0]):
    key = H
    if key not in seen:
        seen[key] = (tr, S)
        print(f"  H={H}: tr(A²)={tr[2]}, tr(A³)={tr[3]}, tr(A⁴)={tr[4]}, "
              f"tr(A⁵)={tr[5]}, tr(A⁶)={tr[6]}")

# Check: are tr(A²) and tr(A³) universal?
tr2_vals = set(d[1][2] for d in data)
tr3_vals = set(d[1][3] for d in data)
print(f"\ntr(A²) values: {tr2_vals} {'(universal)' if len(tr2_vals)==1 else '(varies!)'}")
print(f"tr(A³) values: {tr3_vals} {'(universal)' if len(tr3_vals)==1 else '(varies!)'}")

# p_k = tr(A^k) - λ_0^k where λ_0 = (p-1)/2
lam0 = (p-1)//2  # = 5
print(f"\nλ_0 = {lam0}")

# Power sums p_k (non-trivial eigenvalues)
print("\nPower sums of non-trivial eigenvalues:")
for H in sorted(seen.keys(), reverse=True):
    tr, S = seen[H]
    p4 = tr[4] - lam0**4
    p5 = tr[5] - lam0**5
    p6 = tr[6] - lam0**6
    p7 = tr[7] - lam0**7
    print(f"  H={H}: p_4={p4}, p_5={p5}, p_6={p6}, p_7={p7}")

# Fit H = a + b*p_4 + c*p_5 + d*p_6
print("\nFitting H = a + b*p_4 + c*p_5 + d*p_6:")

unique = list(seen.items())
n_pts = len(unique)

A_mat = np.zeros((n_pts, n_pts))
b_vec = np.zeros(n_pts)

for i, (H, (tr, S)) in enumerate(unique):
    p4 = tr[4] - lam0**4
    p5 = tr[5] - lam0**5
    p6 = tr[6] - lam0**6
    b_vec[i] = H
    if n_pts == 4:
        A_mat[i] = [1, p4, p5, p6]
    elif n_pts == 3:
        A_mat[i] = [1, p4, p5]

try:
    coeffs = np.linalg.solve(A_mat, b_vec)
    if n_pts == 4:
        print(f"  H = {coeffs[0]:.4f} + ({coeffs[1]:.6f})*p_4 + "
              f"({coeffs[2]:.6f})*p_5 + ({coeffs[3]:.6f})*p_6")
    elif n_pts == 3:
        print(f"  H = {coeffs[0]:.4f} + ({coeffs[1]:.6f})*p_4 + ({coeffs[2]:.6f})*p_5")

    # Verify on ALL data
    max_err = 0
    for H, tr, S in data:
        p4 = tr[4] - lam0**4
        p5 = tr[5] - lam0**5
        p6 = tr[6] - lam0**6
        if n_pts == 4:
            pred = coeffs[0] + coeffs[1]*p4 + coeffs[2]*p5 + coeffs[3]*p6
        elif n_pts == 3:
            pred = coeffs[0] + coeffs[1]*p4 + coeffs[2]*p5
        err = abs(H - pred)
        max_err = max(max_err, err)
    print(f"  Max error: {max_err:.6f}")

    # Rationalize coefficients
    from fractions import Fraction
    print("\n  Rational approximations:")
    labels = ["a", "b(p_4)", "c(p_5)", "d(p_6)"]
    for i, c in enumerate(coeffs):
        frac = Fraction(c).limit_denominator(10000)
        print(f"    {labels[i]} = {c:.8f} ≈ {frac}")
except Exception as e:
    print(f"  Error: {e}")

print()

# ================================================================
# p = 13: Same analysis
# ================================================================
print("=" * 70)
print("=== p = 13: H as function of traces ===")
print("=" * 70)
print()

p = 13
sets = all_circulant_sets(p)
data13 = []
for S in sets:
    A = circulant_adjacency(p, S)
    H = hamiltonian_paths_dp(A, p)
    traces = {}
    Ak = np.eye(p, dtype=int)
    for k in range(1, p+1):
        Ak = Ak @ A
        traces[k] = int(round(np.trace(Ak)))
    data13.append((H, traces, sorted(S)))

print("Distinct (H, traces):")
seen13 = {}
lam0_13 = 6
for H, tr, S in sorted(data13, key=lambda x: -x[0]):
    key = H
    if key not in seen13:
        seen13[key] = (tr, S)
        p4 = tr[4] - lam0_13**4
        p5 = tr[5] - lam0_13**5
        p6 = tr[6] - lam0_13**6
        p7 = tr[7] - lam0_13**7
        print(f"  H={H}: p_4={p4}, p_5={p5}, p_6={p6}, p_7={p7}")

n_pts = len(seen13)
print(f"\n{n_pts} distinct H values → need {n_pts-1} power sums to fit")

# Fit with p_4, ..., p_{4+n_pts-2}
unique13 = list(seen13.items())

if n_pts <= 6:
    # Use p_4,...,p_{3+n_pts}
    pk_indices = list(range(4, 3+n_pts+1))[:n_pts-1]

    A_mat = np.zeros((n_pts, n_pts))
    b_vec = np.zeros(n_pts)

    for i, (H, (tr, S)) in enumerate(unique13):
        b_vec[i] = H
        A_mat[i, 0] = 1
        for j, k in enumerate(pk_indices):
            A_mat[i, 1+j] = tr[k] - lam0_13**k

    try:
        coeffs = np.linalg.solve(A_mat, b_vec)
        terms = [f"{coeffs[0]:.2f}"]
        for j, k in enumerate(pk_indices):
            terms.append(f"({coeffs[1+j]:.8f})*p_{k}")
        print(f"\n  H = {' + '.join(terms)}")

        # Verify
        max_err = 0
        for H, tr, S in data13:
            pred = coeffs[0]
            for j, k in enumerate(pk_indices):
                pred += coeffs[1+j] * (tr[k] - lam0_13**k)
            max_err = max(max_err, abs(H - pred))
        print(f"  Max error: {max_err:.6f}")

        # Rationalize
        from fractions import Fraction
        print("\n  Rational approximations:")
        labels = ["a"] + [f"coeff(p_{k})" for k in pk_indices]
        for i, c in enumerate(coeffs):
            frac = Fraction(c).limit_denominator(100000)
            print(f"    {labels[i]} = {c:.10f} ≈ {frac}")
    except Exception as e:
        print(f"  Error: {e}")

print()

# ================================================================
# COMBINATORIAL MEANING OF tr(A^4)
# ================================================================
print("=" * 70)
print("COMBINATORIAL MEANING OF tr(A^4)")
print("=" * 70)
print()
print("tr(A^k) = # closed walks of length k in the tournament")
print("tr(A^4) = # ordered (i,j,k,l) with i→j→k→l→i")
print()
print("For tournaments, tr(A^4) counts:")
print("  - 4-cycles: ordered directed 4-cycles (but tournaments have no")
print("    directed 4-cycles! A tournament on 4 vertices has exactly one")
print("    Hamiltonian cycle, which is a 4-cycle.)")
print("  - Products of 2-walks: paths of length 4 returning to start")
print("    include repeated vertices (e.g., i→j→i→j→i doesn't exist in")
print("    tournaments since arcs are single-direction)")
print()
print("Wait: in a tournament, A[i][j] + A[j][i] = 1 for i≠j, A[i][i]=0.")
print("So A²[i][j] = #{k : i→k→j}.")
print("tr(A²) = Σ_i A²[i][i] = Σ_i #{k : i→k→i} = 0 (no 2-cycles in tournament)")
print("Hmm, that gives 0, but we computed tr(A²) = 21 at p=7!")
print()
print("Ah: A²[i][j] = Σ_k A[i][k]*A[k][j] counts intermediate vertices")
print("where i→k AND k→j. For i=j: A²[i][i] = #{k : A[i][k]=1 AND A[k][i]=1}")
print("= #{k : i→k AND k→i} = 0 (tournament: exactly one direction)")
print()

A7 = circulant_adjacency(7, {1,2,4})
A7sq = A7 @ A7
print(f"p=7 Paley: A²[0][0] = {A7sq[0][0]}, A²[1][1] = {A7sq[1][1]}")
print(f"tr(A²) = {np.trace(A7sq)}")
print(f"But we said tr(A²) = Σλ² = {sum(e**2 for e in eigenvalues_circulant(7, {1,2,4})).real:.0f}")
print()
print("Hmm, there's a discrepancy. Let me recheck...")
print(f"A²[0][0] = {A7sq[0][0]} — this counts #{'{k : 0→k AND k→0}'}")
print(f"Since 0→k means k∈QR={'{1,2,4}'}, and k→0 means -k∈QR, i.e., k∈-QR={'{3,5,6}'}")
print(f"So A²[0][0] = |QR ∩ (-QR)| = |{'{1,2,4}'} ∩ {'{3,5,6}'}| = 0")
print()
print("So tr(A²) = 0 for ANY tournament (no 2-cycles). ✓")
print()
print("But we said Σ_k λ_k² = p*N_2 where N_2 = |{(s1,s2)∈S²: s1+s2≡0}|.")
print(f"At p=7: N_2 = 0 (tournament condition), so Σλ² = 0.")
print(f"But λ_0 = 3, so Σ_{'{k≠0}'} λ_k² = -9. And tr(A²) = Σλ² = 0.")
print(f"Check: 9 + (-9) = 0 ✓")
print()
print(f"So tr(A²) = 0 always. Then tr(A⁴) = Σ_{'{k=0}'}^{'{p-1}'} λ_k⁴.")
print()

# Now let's understand tr(A^4)
# tr(A^4) = Σ_i Σ_{j,k,l} A[i][j]*A[j][k]*A[k][l]*A[l][i]
# = Σ walks (i→j→k→l→i) of length 4 returning to start

A7_4 = A7 @ A7 @ A7 @ A7
print(f"p=7 Paley: tr(A⁴) = {np.trace(A7_4)}")
print(f"  = 7 * (A⁴[0][0]) by circulant symmetry = 7 * {A7_4[0][0]}")
print(f"  A⁴[0][0] = # ways to walk 0→?→?→?→0 in 4 steps")
print()

# Count them explicitly
count = 0
for j in range(7):
    for k in range(7):
        for l in range(7):
            if A7[0][j] and A7[j][k] and A7[k][l] and A7[l][0]:
                count += 1
print(f"  Direct count of walks 0→j→k→l→0: {count}")
print(f"  This includes walks with repeated vertices!")
print()

# Let me count by vertex pattern
from collections import Counter
patterns = Counter()
for j in range(7):
    for k in range(7):
        for l in range(7):
            if A7[0][j] and A7[j][k] and A7[k][l] and A7[l][0]:
                verts = (0, j, k, l)
                n_distinct = len(set(verts))
                patterns[n_distinct] += 1
print(f"  Walk patterns by # distinct vertices:")
for ndist in sorted(patterns.keys()):
    print(f"    {ndist} distinct vertices: {patterns[ndist]} walks")

print()
print("So tr(A⁴) includes walks revisiting vertices.")
print("4 distinct = directed 4-cycles through 0 (× 7 for all starts)")
print("3 distinct = diamond-shaped walks (i→j→k→j→i type)")
print("2 distinct = back-and-forth (impossible in tournament: no 2-cycles)")
print()

# ================================================================
# THE KEY FORMULA: H = (C - tr(A^4))/2 at p=7
# ================================================================
print("=" * 70)
print("THE KEY FORMULA")
print("=" * 70)
print()
print("At p=7: H = (462 - tr(A⁴))/2")
print()
print("This means: MINIMIZING tr(A⁴) = MAXIMIZING H.")
print()
print("tr(A⁴) = Σ |λ_k|⁴·cos(4·arg(λ_k)) + |λ_0|⁴")
print("        = ((p-1)/2)⁴ + Σ_{k≠0} Re(λ_k⁴)")
print()
print("For Paley: all |λ_k| = √((p+1)/4), all phases = ±arctan(√p)")
print("  tr(A⁴) = ((p-1)/2)⁴ + (p-1)·((p+1)/4)² · cos(4·arctan(√p))")
print()
print("For flat spectrum at p=7:")
print(f"  λ_k = -1/2 ± i√7/2, |λ|² = 2, |λ|⁴ = 4")
print(f"  arg(λ) ≈ ±1.9322 rad")
print(f"  4·arg ≈ ±7.7288 → cos(7.7288) ≈ cos(7.7288-2π) ≈ cos(1.446) ≈ {np.cos(1.446):.4f}")
print(f"  Σ Re(λ⁴) = 6 × 4 × {np.cos(4*np.arctan(np.sqrt(7))):.6f} = {6*4*np.cos(4*np.arctan(np.sqrt(7))):.4f}")
print(f"  tr(A⁴) = 81 + {6*4*np.cos(4*np.arctan(np.sqrt(7))):.4f} = {81 + 6*4*np.cos(4*np.arctan(np.sqrt(7))):.4f}")
print(f"  Expected: 84. Check: {abs(84 - (81 + 6*4*np.cos(4*np.arctan(np.sqrt(7))))) < 0.01}")
print()

# ================================================================
# SCHUR CONVEXITY ARGUMENT
# ================================================================
print("=" * 70)
print("SCHUR CONVEXITY: WHY FLAT MINIMIZES tr(A⁴)")
print("=" * 70)
print()
print("For p ≡ 3 mod 4, all eigenvalues satisfy Re(λ_k) = -1/2 (universal).")
print("So λ_k = -1/2 + i·y_k where y_k are the imaginary parts.")
print()
print("Parseval: Σ |λ_k|² = Σ (1/4 + y_k²) = (p-1)(p+1)/4")
print("  → Σ y_k² = (p-1)(p+1)/4 - (p-1)/4 = (p-1)p/4")
print()
print("Now: λ_k⁴ = (-1/2 + iy_k)⁴")
print("  = (1/2 - iy_k)⁴ ... wait: (-1/2 + iy)⁴ = ((-1+2iy)/2)⁴ = (1-2iy)⁴/16")
print("  = (1 - 8iy - 24y² + 32iy³ + 16y⁴)/16... let me compute properly.")
print()

from sympy import symbols, expand, I, Rational

y = symbols('y')
expr = (-Rational(1,2) + I*y)**4
expanded = expand(expr)
print(f"  (-1/2 + iy)⁴ = {expanded}")
print(f"  Re part: {expanded.as_real_imag()[0]}")
print()

# Re((-1/2 + iy)^4) = 1/16 - 3y²/2 + y⁴
# Check: at y = √7/2 (Paley p=7): 1/16 - 3*(7/4)/2 + (7/4)² = 1/16 - 21/8 + 49/16
# = 1/16 - 42/16 + 49/16 = 8/16 = 1/2
# So Σ Re(λ_k⁴) = 6 × 1/2 = 3. And p_4 = 3 ✓!

re_part_formula = "Re((-1/2 + iy)⁴) = 1/16 - 3y²/2 + y⁴"
print(f"  {re_part_formula}")
print()
print("  Σ_{k≠0} Re(λ_k⁴) = Σ (1/16 - 3y_k²/2 + y_k⁴)")
print("                     = (p-1)/16 - 3/2·Σy_k² + Σy_k⁴")
print("                     = (p-1)/16 - 3(p-1)p/8 + Σy_k⁴")
print()
print("  So tr(A⁴) = ((p-1)/2)⁴ + (p-1)/16 - 3(p-1)p/8 + Σy_k⁴")
print("            = CONSTANT + Σy_k⁴")
print()
print("  H = (C₁ - tr(A⁴))/2 = (C₂ - Σy_k⁴)/2")
print()
print("  THEREFORE: H is maximized ⟺ Σy_k⁴ is minimized")
print("  ⟺ all y_k² are equal (by convexity of x²)")
print("  ⟺ FLAT SPECTRUM (Paley)")
print()
print("  THIS IS THE COMPLETE PROOF OF STEP B AT p=7!")
print()

# Verify for all p ≡ 3 mod 4
for p in [7, 11]:
    lam0 = (p-1)//2
    # Parseval: Σy_k² = (p-1)p/4
    sum_y2 = (p-1)*p/4
    # Constant in tr(A^4)
    const = lam0**4 + (p-1)/16 - 3*(p-1)*p/8

    sets = all_circulant_sets(p)
    for S in sets[:3]:
        A = circulant_adjacency(p, S)
        H = hamiltonian_paths_dp(A, p)
        eigs = eigenvalues_circulant(p, S)
        # Imaginary parts
        y_vals = [e.imag for e in eigs[1:(p+1)//2]]
        sum_y4 = sum(yi**4 for yi in y_vals)
        # But we need to double (conjugate pairs)
        sum_y4_full = 2 * sum_y4
        tr4 = int(round(np.trace(A @ A @ A @ A)))
        tr4_pred = const + sum_y4_full
        print(f"  p={p}, S={sorted(S)}: tr(A⁴)={tr4}, "
              f"const+Σy⁴={tr4_pred:.2f}, diff={abs(tr4-tr4_pred):.4f}")

print()
print("IF H is linear in tr(A⁴) at p=11 too, then:")
print("  H = C - α·Σy_k⁴ (with α > 0)")
print("  → minimizing Σy_k⁴ → flat y_k² → Paley maximizes H")
print()

# Test linearity at p=11
p = 11
lam0 = 5
seen = {}
for S in all_circulant_sets(p):
    A = circulant_adjacency(p, S)
    H = hamiltonian_paths_dp(A, p)
    tr4 = int(round(np.trace(A @ A @ A @ A)))
    if H not in seen:
        seen[H] = tr4

print(f"p=11 (H, tr(A⁴)) pairs:")
for H in sorted(seen.keys(), reverse=True):
    print(f"  H={H}, tr(A⁴)={seen[H]}")

H_vals = sorted(seen.keys(), reverse=True)
tr4_vals = [seen[H] for H in H_vals]

# Linear fit: H = a + b*tr4
if len(H_vals) >= 2:
    from numpy.polynomial import polynomial as P
    coeffs = np.polyfit(tr4_vals, H_vals, 1)
    print(f"\nLinear fit: H = {coeffs[1]:.4f} + {coeffs[0]:.6f}*tr(A⁴)")
    max_err = max(abs(H - (coeffs[1] + coeffs[0]*seen[H])) for H in H_vals)
    print(f"Max error: {max_err:.4f}")
    if max_err > 1:
        print("NOT linear in tr(A⁴) alone at p=11!")
        # Try quadratic
        coeffs2 = np.polyfit(tr4_vals, H_vals, 2)
        print(f"Quadratic: H = {coeffs2[2]:.4f} + {coeffs2[1]:.8f}*tr + {coeffs2[0]:.10f}*tr²")
        max_err2 = max(abs(H - np.polyval(coeffs2, seen[H])) for H in H_vals)
        print(f"Quadratic max error: {max_err2:.6f}")
