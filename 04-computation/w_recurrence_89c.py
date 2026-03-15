#!/usr/bin/env python3
"""
w_recurrence_89c.py — Find recurrence for W(n) = Σ_{σ∈NUD(n)} 2^{adj1(σ)}
opus-2026-03-14-S89c

Strategy: compute W(n) exactly for small n, then find a recurrence
using the "guess and verify" approach.

W(n) values from nud_weight.c:
  W(3)=8, W(4)=32, W(5)=158, W(6)=928, W(7)=6350, W(8)=49752
  W(9)=439670, W(10)=4327904, W(11)=46963358

Also try to prove CV² = 2/n analytically by finding W(n)/n!.
"""

from fractions import Fraction
from itertools import permutations
from math import factorial

def compute_W(n):
    """Exact W(n) by brute force enumeration."""
    W = 0
    nud_count = 0
    for perm in permutations(range(n)):
        # Check NUD: no position with perm[i+1] = perm[i] - 1
        is_nud = True
        adj1 = 0
        for i in range(n-1):
            if perm[i+1] == perm[i] - 1:
                is_nud = False
                break
            if perm[i+1] == perm[i] + 1:
                adj1 += 1
        if is_nud:
            nud_count += 1
            W += 2**adj1
    return W, nud_count

# Compute W(n) for small n
print("Computing W(n) exactly...")
W_vals = {}
for n in range(2, 11):
    W, nud = compute_W(n)
    nf = factorial(n)
    cv2 = Fraction(W, nf) - 1
    print(f"  W({n}) = {W}, NUD={nud}, n!={nf}, W/n!={float(Fraction(W,nf)):.10f}, CV²={float(cv2):.10f}")
    W_vals[n] = W

print("\n" + "="*60)
print("RECURRENCE SEARCH")
print("="*60)

# Try W(n) = a(n)*W(n-1) + b(n)*W(n-2)
# where a(n), b(n) are polynomial in n
print("\nTrying W(n) = a·W(n-1) + b·W(n-2):")
for n in range(4, 11):
    # W(n) = a*W(n-1) + b*W(n-2)
    # Two equations from consecutive n would overdetermine a,b if they're constants
    pass

# For each n, compute r(n) = W(n)/W(n-1)
print("\nRatios W(n)/W(n-1):")
for n in range(3, 11):
    r = Fraction(W_vals[n], W_vals[n-1])
    print(f"  W({n})/W({n-1}) = {float(r):.10f} = {r}")

# Try to express as W(n) = (an+b)*W(n-1) + (cn+d)*W(n-2)
# System: for each n >= 4
# W(n) - (an+b)*W(n-1) - (cn+d)*W(n-2) = 0
# 4 unknowns, need 4 equations (n=4,5,6,7)
print("\nTrying W(n) = (an+b)·W(n-1) + (cn+d)·W(n-2):")

from fractions import Fraction as F

# Set up linear system
# W(n) = (an+b)*W(n-1) + (cn+d)*W(n-2)
# Rearrange: a*n*W(n-1) + b*W(n-1) + c*n*W(n-2) + d*W(n-2) = W(n)
# Variables: [a, b, c, d]
# Row for each n: [n*W(n-1), W(n-1), n*W(n-2), W(n-2)] = W(n)

import numpy as np

ns = [4, 5, 6, 7]
A = []
rhs = []
for n in ns:
    row = [n*W_vals[n-1], W_vals[n-1], n*W_vals[n-2], W_vals[n-2]]
    A.append(row)
    rhs.append(W_vals[n])

A = np.array(A, dtype=float)
rhs = np.array(rhs, dtype=float)
sol = np.linalg.solve(A, rhs)
print(f"  a={sol[0]:.10f}, b={sol[1]:.10f}, c={sol[2]:.10f}, d={sol[3]:.10f}")

# Verify with exact fractions
from sympy import symbols, Eq, solve as sym_solve, Rational, Matrix

a, b, c, d = symbols('a b c d')
eqs = []
for n in ns:
    eq = Eq(a*n*W_vals[n-1] + b*W_vals[n-1] + c*n*W_vals[n-2] + d*W_vals[n-2], W_vals[n])
    eqs.append(eq)

exact_sol = sym_solve(eqs, [a, b, c, d])
print(f"  Exact: a={exact_sol[a]}, b={exact_sol[b]}, c={exact_sol[c]}, d={exact_sol[d]}")

# Verify on n=8,9,10
print("\n  Verification:")
aa, bb, cc, dd = exact_sol[a], exact_sol[b], exact_sol[c], exact_sol[d]
for n in range(4, 11):
    predicted = (aa*n + bb)*W_vals[n-1] + (cc*n + dd)*W_vals[n-2]
    actual = W_vals[n]
    match = "✓" if predicted == actual else "✗"
    print(f"    n={n}: predicted={predicted}, actual={actual} {match}")

# If that doesn't work, try 3-term recurrence
# W(n) = (an+b)*W(n-1) + (cn+d)*W(n-2) + (en+f)*W(n-3)
if exact_sol[a] != aa:  # just continue anyway
    pass

print("\n\nTrying W(n) = (an+b)·W(n-1) + (cn+d)·W(n-2) + (en+f)·W(n-3):")
ns = [5, 6, 7, 8, 9, 10]
a2, b2, c2, d2, e2, f2 = symbols('a2 b2 c2 d2 e2 f2')
eqs2 = []
for n in ns:
    eq = Eq(a2*n*W_vals[n-1] + b2*W_vals[n-1] + c2*n*W_vals[n-2] + d2*W_vals[n-2]
            + e2*n*W_vals[n-3] + f2*W_vals[n-3], W_vals[n])
    eqs2.append(eq)

exact_sol2 = sym_solve(eqs2, [a2, b2, c2, d2, e2, f2])
if exact_sol2:
    print(f"  a={exact_sol2[a2]}, b={exact_sol2[b2]}, c={exact_sol2[c2]}, d={exact_sol2[d2]}, e={exact_sol2[e2]}, f={exact_sol2[f2]}")

# Also try: define R(n) = W(n)/n! and look for recurrence
print("\n" + "="*60)
print("RATIO R(n) = W(n)/n! ANALYSIS")
print("="*60)

R = {}
for n in range(2, 11):
    R[n] = Fraction(W_vals[n], factorial(n))
    print(f"  R({n}) = W({n})/n! = {R[n]} = {float(R[n]):.12f}")

# R(n) - 1 = CV²
print("\n  R(n) - 1 (= CV²):")
for n in range(2, 11):
    cv2 = R[n] - 1
    print(f"  CV²({n}) = {cv2} = {float(cv2):.12f}, n×CV² = {float(n*cv2):.12f}")

# Try: R(n) = 1 + 2/n + c/n² + ...
# So n²(R(n) - 1 - 2/n) should converge
print("\n  n²(R(n) - 1 - 2/n):")
for n in range(3, 11):
    val = R[n] - 1 - Fraction(2, n)
    scaled = n*n * val
    print(f"    n={n}: {float(scaled):.10f} (exact: {scaled})")

# Maybe try n²(R(n) - 1 - 2/n) → some constant
# And n³(R(n) - 1 - 2/n - c/n²) → some other constant
print("\n  Successive asymptotic extraction:")
print("  n²(CV² - 2/n):")
coeffs = []
for n in range(3, 11):
    cv2 = R[n] - 1
    val = n*n * (cv2 - Fraction(2, n))
    coeffs.append((n, float(val)))
    print(f"    n={n}: {float(val):.10f}")

# The second coefficient seems to converge
# Try extracting: CV² = 2/n + a₂/n² + a₃/n³ + ...
# where a₂ = lim n²(CV² - 2/n)
print("\n  Guessing a₂ from trend...")
# Values should converge to a₂

# Try to fit CV² = 2/n + a/n² + b/n³ using n=8,9,10
from sympy import Rational as Rat

n8, n9, n10 = 8, 9, 10
cv2_8 = R[8] - 1
cv2_9 = R[9] - 1
cv2_10 = R[10] - 1

a_sym, b_sym = symbols('a_coeff b_coeff')
eq_8 = Eq(cv2_8, Fraction(2, 8) + a_sym/64 + b_sym/512)
eq_9 = Eq(cv2_9, Fraction(2, 9) + a_sym/81 + b_sym/729)
eq_10 = Eq(cv2_10, Fraction(2, 10) + a_sym/100 + b_sym/1000)

sol_ab = sym_solve([eq_8, eq_9], [a_sym, b_sym])
if sol_ab:
    print(f"  From n=8,9: a₂ = {sol_ab[a_sym]} = {float(sol_ab[a_sym]):.10f}")
    print(f"              a₃ = {sol_ab[b_sym]} = {float(sol_ab[b_sym]):.10f}")
    # Check on n=10
    pred_10 = Fraction(2, 10) + sol_ab[a_sym]/100 + sol_ab[b_sym]/1000
    print(f"  Predicted CV²(10) = {float(pred_10):.12f}")
    print(f"  Actual    CV²(10) = {float(cv2_10):.12f}")

# Try from n=9,10
sol_ab2 = sym_solve([eq_9, eq_10], [a_sym, b_sym])
if sol_ab2:
    print(f"  From n=9,10: a₂ = {sol_ab2[a_sym]} = {float(sol_ab2[a_sym]):.10f}")
    print(f"               a₃ = {sol_ab2[b_sym]} = {float(sol_ab2[b_sym]):.10f}")

print("\n" + "="*60)
print("ANALYTIC APPROACH: EGF CONNECTION")
print("="*60)

# Key insight: NUD permutations avoid the consecutive pattern "descent by 1"
# Their EGF is exp(-x)/(1-x)² (= EGF of A000255)
#
# The bivariate EGF tracking unit ascents (successions) should be:
# F(x,z) = Σ_n (Σ_{σ∈NUD(n)} z^{adj1(σ)}) x^n/n!
#
# At z=1: F(x,1) = exp(-x)/(1-x)² (just counting NUD perms)
# At z=2: F(x,2) = Σ W(n) x^n/n! (what we want)
#
# For UNRESTRICTED permutations, the EGF tracking successions is:
# G(x,z) = exp((z-1)x)/(1-x)   [classical result]
# At z=1: G(x,1) = 1/(1-x) (all perms)
# At z=2: G(x,2) = exp(x)/(1-x) (EGF of Σ 2^{succ(σ)})

# Verify unrestricted succession EGF
print("\nVerifying unrestricted succession EGF G(x,z)=exp((z-1)x)/(1-x):")
for n in range(2, 9):
    # Σ_{σ∈S_n} z^{succ(σ)} at z=2 should be n! × [x^n] exp(x)/(1-x)
    # exp(x)/(1-x) = Σ x^n/n! × Σ x^n = Σ_n (Σ_{k=0}^n n!/k!) x^n/n!
    # So coefficient of x^n/n! in exp(x)/(1-x) = Σ_{k=0}^n 1/k!  NO...
    # Actually exp(x)/(1-x) = (Σ x^k/k!) × (Σ x^j)
    # [x^n] = Σ_{k=0}^n 1/k!
    # So n! × [x^n] = n! × Σ_{k=0}^n 1/k! = Σ_{k=0}^n n!/k!
    predicted = sum(factorial(n)//factorial(k) for k in range(n+1))

    # Brute force
    actual = 0
    for perm in permutations(range(n)):
        succ = sum(1 for i in range(n-1) if perm[i+1] == perm[i] + 1)
        actual += 2**succ

    print(f"  n={n}: predicted={predicted}, actual={actual}, match={'✓' if predicted==actual else '✗'}")

# Now for NUD perms with successions tracked:
# The NUD constraint forbids the "descent by 1" pattern.
# What is F(x,z) for NUD perms with z^{succession}?
#
# Think of it as a transfer matrix on a graph:
# States: the value at current position (0..n-1)
# Transitions: from v to u ≠ v-1, with weight z if u = v+1
#
# For a fixed alphabet [n], this is equivalent to paths in a weighted graph.
# But the EGF approach considers all n simultaneously.
#
# Alternative: use the "barred pattern" or "marked succession" technique.
#
# Consider the composition of a permutation into maximal ascending runs
# of consecutive values. In a NUD permutation, descents by 1 are forbidden,
# but longer descents are allowed.

# Let me try a different approach: verify if F(x,2) has a nice closed form
print("\nF(x,2) = Σ W(n) x^n/n! coefficients:")
from sympy import Rational, symbols, series, exp, simplify

x = symbols('x')
# Build the series from data
F2_series = sum(Fraction(W_vals[n], factorial(n)) * x**n for n in range(2, 11))
print(f"  F(x,2) = {F2_series}")

# Compare with some candidate closed forms:
# Guess 1: exp(-x)/(1-x)² × exp(x) = 1/(1-x)²
# At z=2 with the z-1 shift: exp((2-1)x)·exp(-x)/(1-x)² = 1/(1-x)²
# 1/(1-x)² has [x^n] = n+1, so n!×[x^n/n!] = n+1... no that's [x^n],
# the EGF coeff would be (n+1)/n!... not matching

# Guess 2: exp(x)/(1-x)³
# [x^n] exp(x)/(1-x)³ = Σ_{k=0}^n C(n-k+2,2)/k!
# n! × this = Σ_{k=0}^n n!/k! × C(n-k+2,2)

# Let me just compute several candidate EGF coefficients
candidates = {
    "1/(1-x)^2": lambda n: 1,  # [x^n/n!] = (n+1)/n!, wait no
}

# Actually let me be more careful
# [x^n] f(x) means the coefficient of x^n
# For EGF: a_n = n! × [x^n] f(x)
#
# 1/(1-x)^2: [x^n] = n+1, so a_n = (n+1)·n!  ... too big
# exp(x)/(1-x)^2: [x^n] = Σ_{k=0}^n (n-k+1)/k! = Σ_{j=0}^n (j+1)/((n-j)!)
#   a_n = n! × Σ_{j=0}^n (j+1)/(n-j)! = Σ_{j=0}^n (j+1)·n!/(n-j)!

print("\nComparing with exp(x)/(1-x)²:")
for n in range(2, 11):
    val = sum((j+1) * factorial(n) // factorial(n-j) for j in range(n+1))
    print(f"  n={n}: exp(x)/(1-x)² gives {val}, W(n) = {W_vals[n]}, ratio = {val/W_vals[n]:.6f}")

print("\nComparing with exp(ax)/(1-x)² for various a:")
for a_val in [Fraction(1,2), Fraction(1,3), Fraction(2,3), Fraction(-1,2)]:
    print(f"\n  a = {a_val}:")
    for n in range(2, 7):
        # [x^n] exp(ax)/(1-x)² = Σ_{k=0}^n (n-k+1) a^k / k!
        val = sum((n-k+1) * a_val**k / factorial(k) for k in range(n+1))
        egf_coeff = val * factorial(n)
        print(f"    n={n}: {float(egf_coeff):.4f} vs W={W_vals[n]}")

# Try to find R(n) = W(n)/n! recurrence
print("\n" + "="*60)
print("RECURRENCE FOR R(n) = W(n)/n!")
print("="*60)

# If W(n) = (an+b)W(n-1) + (cn+d)W(n-2), then
# R(n)·n! = (an+b)·R(n-1)·(n-1)! + (cn+d)·R(n-2)·(n-2)!
# R(n) = (an+b)/n · R(n-1) + (cn+d)/(n(n-1)) · R(n-2)

# Use the exact recurrence we found above
print(f"\nUsing exact solution: a={exact_sol[a]}, b={exact_sol[b]}, c={exact_sol[c]}, d={exact_sol[d]}")
print(f"\nR(n) recurrence:")
print(f"  R(n) = ({exact_sol[a]}n + {exact_sol[b]})/n × R(n-1) + ({exact_sol[c]}n + {exact_sol[d]})/(n(n-1)) × R(n-2)")

# Simplify the coefficients
aa_val = exact_sol[a]
bb_val = exact_sol[b]
cc_val = exact_sol[c]
dd_val = exact_sol[d]

print(f"\n  Leading coefficient of R(n-1): ({aa_val}n+{bb_val})/n = {aa_val} + {bb_val}/n")
print(f"  Leading coefficient of R(n-2): ({cc_val}n+{dd_val})/(n(n-1))")

# If a=1, then R(n) ≈ R(n-1) + correction → R converges
# R(∞) = 1 (since CV² → 0)

print("\nDone!")
