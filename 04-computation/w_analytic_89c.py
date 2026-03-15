#!/usr/bin/env python3
"""
w_analytic_89c.py — Analytical approach to proving CV² = 2/n
opus-2026-03-14-S89c

Key identity: CV² = W(n)/n! - 1
where W(n) = Σ_{σ∈NUD(n)} 2^{adj1(σ)}.

Define B(n) = (n+2)(n-1)! - W(n). Then n·CV² = 2 - B(n)/(n-1)!
We need B(n)/(n-1)! → 0.

B(n) sequence: 2, 4, 10, 32, 130, 648, 3850, 26656
"""

from fractions import Fraction
from itertools import permutations
from math import factorial
import sys

# Compute W(n) exactly
def compute_W(n):
    W = 0
    nud_count = 0
    for perm in permutations(range(n)):
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

# Also compute decomposed statistics
def compute_detailed(n):
    """Count NUD perms by (adj1, adj_desc) where adj_desc would be unit descents if allowed."""
    from collections import Counter
    succ_hist = Counter()  # histogram of succession count in NUD perms
    total_succ = 0  # total successions across all NUD perms
    total_succ_sq = 0
    nud_count = 0
    W = 0

    for perm in permutations(range(n)):
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
            succ_hist[adj1] += 1
            total_succ += adj1
            total_succ_sq += adj1**2

    return W, nud_count, succ_hist, total_succ, total_succ_sq

# Also compute for ALL perms (not just NUD)
def compute_all_perms(n):
    """Succession statistics for all perms of [n]."""
    from collections import Counter
    succ_hist = Counter()
    total = 0
    total_succ = 0

    for perm in permutations(range(n)):
        adj1 = sum(1 for i in range(n-1) if perm[i+1] == perm[i] + 1)
        desc1 = sum(1 for i in range(n-1) if perm[i+1] == perm[i] - 1)
        succ_hist[(adj1, desc1)] += 1
        total += 1
        total_succ += adj1

    return succ_hist, total

print("="*70)
print("DETAILED SUCCESSION ANALYSIS IN NUD PERMUTATIONS")
print("="*70)

W_vals = {}
for n in range(2, 10):
    W, nud, succ_hist, ts, ts2 = compute_detailed(n)
    W_vals[n] = W
    nf = factorial(n)

    print(f"\nn={n}: |NUD|={nud}, W={W}, n!={nf}")
    print(f"  E[succ | NUD] = {Fraction(ts, nud)} = {ts/nud:.6f}")
    print(f"  E[succ²| NUD] = {Fraction(ts2, nud)} = {ts2/nud:.6f}")
    print(f"  Var[succ|NUD] = {ts2/nud - (ts/nud)**2:.6f}")
    print(f"  Succession distribution in NUD perms:")
    for k in sorted(succ_hist):
        cnt = succ_hist[k]
        print(f"    succ={k}: {cnt} ({cnt/nud*100:.1f}%), contribution to W: {cnt*2**k} ({cnt*2**k/W*100:.1f}%)")

# B(n) = (n+2)(n-1)! - W(n)
print("\n" + "="*70)
print("B(n) = (n+2)(n-1)! - W(n) SEQUENCE ANALYSIS")
print("="*70)

B_vals = {}
for n in range(3, 10):
    B = (n+2)*factorial(n-1) - W_vals[n]
    B_vals[n] = B
    print(f"  B({n}) = {B}")

# Search for recurrence of B(n) with n-dependent coefficients
print("\nTrying B(n) = (an+b)B(n-1) + (cn+d)B(n-2):")

from sympy import symbols, Eq, solve, Rational

a, b, c, d = symbols('a b c d')
eqs = []
for n in [5, 6, 7, 8]:
    eqs.append(Eq(B_vals[n], (a*n+b)*B_vals[n-1] + (c*n+d)*B_vals[n-2]))

sol = solve(eqs, [a, b, c, d])
print(f"  Solution: a={sol[a]}, b={sol[b]}, c={sol[c]}, d={sol[d]}")

# Verify
print("  Verification:")
for n in range(5, 10):
    pred = (sol[a]*n + sol[b]) * B_vals[n-1] + (sol[c]*n + sol[d]) * B_vals[n-2]
    match = "✓" if pred == B_vals[n] else "✗"
    print(f"    B({n}): predicted={pred}, actual={B_vals[n]} {match}")

# Check W(n) directly for recurrence with n-dependent coefficients
print("\n" + "="*70)
print("W(n) RECURRENCE WITH n-DEPENDENT COEFFICIENTS")
print("="*70)

# Try W(n) = (an² + bn + c)W(n-1) + (dn² + en + f)W(n-2)
a1, b1, c1, d1, e1, f1 = symbols('a1 b1 c1 d1 e1 f1')
eqs2 = []
for n in [4, 5, 6, 7, 8, 9]:
    eqs2.append(Eq(W_vals[n],
                   (a1*n**2 + b1*n + c1)*W_vals[n-1] +
                   (d1*n**2 + e1*n + f1)*W_vals[n-2]))

sol2 = solve(eqs2, [a1, b1, c1, d1, e1, f1])
if sol2:
    print(f"  Solution: a={sol2[a1]}, b={sol2[b1]}, c={sol2[c1]}")
    print(f"            d={sol2[d1]}, e={sol2[e1]}, f={sol2[f1]}")

# Simpler: try W(n) = (n+α)W(n-1) + β·(n-1)·W(n-2)
# (inspired by A000255 recurrence a(n) = n·a(n-1) + (n-1)·a(n-2))
print("\nTrying W(n) = (n+α)W(n-1) + β(n-1)W(n-2):")
alpha, beta = symbols('alpha beta')
eqs3 = []
for n in [4, 5]:
    eqs3.append(Eq(W_vals[n], (n+alpha)*W_vals[n-1] + beta*(n-1)*W_vals[n-2]))

sol3 = solve(eqs3, [alpha, beta])
if sol3:
    print(f"  α={sol3[alpha]}, β={sol3[beta]}")
    for n in range(4, 10):
        pred = (n+sol3[alpha])*W_vals[n-1] + sol3[beta]*(n-1)*W_vals[n-2]
        match = "✓" if pred == W_vals[n] else "✗"
        print(f"    W({n}): predicted={pred}, actual={W_vals[n]} {match}")

# More general: W(n) = (n+α)W(n-1) + (β·n+γ)W(n-2)
alpha2, beta2, gamma2 = symbols('alpha2 beta2 gamma2')
eqs4 = []
for n in [4, 5, 6]:
    eqs4.append(Eq(W_vals[n], (n+alpha2)*W_vals[n-1] + (beta2*n+gamma2)*W_vals[n-2]))

sol4 = solve(eqs4, [alpha2, beta2, gamma2])
if sol4:
    print(f"\n  α={sol4[alpha2]}, β={sol4[beta2]}, γ={sol4[gamma2]}")
    for n in range(4, 10):
        pred = (n+sol4[alpha2])*W_vals[n-1] + (sol4[beta2]*n+sol4[gamma2])*W_vals[n-2]
        match = "✓" if pred == W_vals[n] else "✗"
        print(f"    W({n}): predicted={pred}, actual={W_vals[n]} {match}")

# Try: W(n) = n·W(n-1) + (n-1)·W(n-2) - correction
# This would match A000255 recurrence if correction=0
print("\nResiduals from A000255-type recurrence W(n) - n·W(n-1) - (n-1)·W(n-2):")
for n in range(4, 10):
    res = W_vals[n] - n*W_vals[n-1] - (n-1)*W_vals[n-2]
    print(f"  n={n}: residual = {res}")

# Maybe: W(n) = n·W(n-1) + (n-3)·W(n-2)?
print("\nResiduals from W(n) - n·W(n-1) - (n-3)·W(n-2):")
for n in range(4, 10):
    res = W_vals[n] - n*W_vals[n-1] - (n-3)*W_vals[n-2]
    print(f"  n={n}: residual = {res}")

# Try various simple (an+b)W(n-1) + (cn+d)W(n-2) until one sticks
# Brute force search over small integer coefficients
print("\nBrute force search for W(n) = (an+b)W(n-1) + (cn+d)W(n-2):")
best = None
best_error = float('inf')
for a_try in range(-3, 4):
    for b_try in range(-5, 6):
        for c_try in range(-3, 4):
            for d_try in range(-5, 6):
                error = 0
                for n in range(4, 10):
                    pred = (a_try*n + b_try)*W_vals[n-1] + (c_try*n + d_try)*W_vals[n-2]
                    error += abs(pred - W_vals[n])
                if error < best_error:
                    best_error = error
                    best = (a_try, b_try, c_try, d_try)

if best_error > 0:
    print(f"  Best: a={best[0]}, b={best[1]}, c={best[2]}, d={best[3]}, error={best_error}")
else:
    print(f"  EXACT: a={best[0]}, b={best[1]}, c={best[2]}, d={best[3]}")
    for n in range(4, 10):
        pred = (best[0]*n + best[1])*W_vals[n-1] + (best[2]*n + best[3])*W_vals[n-2]
        print(f"    W({n}): predicted={pred}, actual={W_vals[n]}")

# Maybe we need a 3-term recurrence with n-linear coefficients
print("\n3-term recurrence search W(n) = (an+b)W(n-1) + (cn+d)W(n-2) + (en+f)W(n-3):")
a3, b3, c3, d3, e3, f3 = symbols('a3 b3 c3 d3 e3 f3')
eqs5 = []
for n in [5, 6, 7, 8, 9]:
    eqs5.append(Eq(W_vals[n], (a3*n+b3)*W_vals[n-1] + (c3*n+d3)*W_vals[n-2] + (e3*n+f3)*W_vals[n-3]))

sol5 = solve(eqs5[:5], [a3, b3, c3, d3, e3, f3])
# 5 equations, 6 unknowns — underdetermined. Need n=10 too or fix one param.
# Let me add n=5..9 = 5 equations... but 6 unknowns. Fix a3=1:
eqs5b = [eq.subs(a3, 1) for eq in eqs5]
sol5b = solve(eqs5b, [b3, c3, d3, e3, f3])
if sol5b:
    print(f"  With a=1: b={sol5b[b3]}, c={sol5b[c3]}, d={sol5b[d3]}, e={sol5b[e3]}, f={sol5b[f3]}")
    # Verify on all n
    for n in range(5, 10):
        pred = (n+sol5b[b3])*W_vals[n-1] + (sol5b[c3]*n+sol5b[d3])*W_vals[n-2] + (sol5b[e3]*n+sol5b[f3])*W_vals[n-3]
        match = "✓" if pred == W_vals[n] else "✗"
        print(f"    W({n}): predicted={pred}, actual={W_vals[n]} {match}")

print("\n" + "="*70)
print("ANALYTIC PROOF STRATEGY")
print("="*70)

print("""
Key insight from Poisson approximation:

For a random permutation of [n], let A = # unit ascents, D = # unit descents.
Both A and D are approximately Poisson(1) for large n, with correlation.

W(n)/n! = E[2^A · 1_{D=0}]

where expectation is over uniform random permutations.

By independence (approximately):
E[2^A · 1_{D=0}] ≈ E[2^A] · P(D=0)

For large n:
- E[2^A] → E_Poisson[2^X] = Σ e^{-1}/k! · 2^k = e^{-1} · e^2 = e
- P(D=0) → e^{-1}

So W(n)/n! → e · e^{-1} = 1. ✓

For the correction term:
- E[A] = (n-1)/n ≈ 1 - 1/n  (prob 1/(n-1) at each of n-1 positions... wait)

Actually, let me be precise about the first moment:
""")

# Exact E[A] and P(D=0) for permutations of [n]
print("Exact statistics for random permutation of [n]:")
for n in range(3, 10):
    total_A = 0
    total_A2 = 0
    total_D = 0
    count_D0 = 0
    count_D0_A = 0  # sum of A when D=0
    count_D0_2A = 0  # sum of 2^A when D=0
    N = factorial(n)

    for perm in permutations(range(n)):
        A = sum(1 for i in range(n-1) if perm[i+1] == perm[i] + 1)
        D = sum(1 for i in range(n-1) if perm[i+1] == perm[i] - 1)
        total_A += A
        total_A2 += A**2
        total_D += D
        if D == 0:
            count_D0 += 1
            count_D0_A += A
            count_D0_2A += 2**A

    EA = Fraction(total_A, N)
    ED = Fraction(total_D, N)
    PD0 = Fraction(count_D0, N)
    EA_given_D0 = Fraction(count_D0_A, count_D0) if count_D0 > 0 else 0
    E2A = Fraction(sum(2**sum(1 for i in range(n-1) if perm[i+1]==perm[i]+1) for perm in permutations(range(n))), N)
    E2A_D0 = Fraction(count_D0_2A, N)  # This is W(n)/n!

    print(f"\n  n={n}:")
    print(f"    E[A] = {EA} = {float(EA):.6f}")
    print(f"    E[D] = {ED} = {float(ED):.6f}")
    print(f"    P(D=0) = {PD0} = {float(PD0):.6f} (1/e = {1/2.718281828:.6f})")
    print(f"    E[2^A] = {E2A} = {float(E2A):.6f} (e = 2.718282)")
    print(f"    E[2^A|D=0] = {Fraction(count_D0_2A, count_D0) if count_D0 else 'N/A'} = {count_D0_2A/count_D0 if count_D0 else 0:.6f}")
    print(f"    E[2^A·1(D=0)] = {E2A_D0} = {float(E2A_D0):.6f} = W/n!")
    print(f"    E[2^A]·P(D=0) = {float(E2A*PD0):.6f}")
    print(f"    Ratio W(n)/n! / (E[2^A]·P(D=0)) = {float(E2A_D0/(E2A*PD0)):.8f}")
    print(f"    Covariance factor = {float(E2A_D0 - E2A*PD0):.8f}")

print("\n" + "="*70)
print("CUMULANT EXPANSION")
print("="*70)

# The key quantity is the covariance between 2^A and 1(D=0).
# Write E[2^A · 1(D=0)] = E[2^A] · P(D=0) + Cov(2^A, 1(D=0))
# So W(n)/n! = e·e^{-1} + Cov = 1 + Cov
# If Cov = 2/n + o(1/n), we're done.

# More precisely, for each position j:
# - P(unit ascent at j) = 1/n  (need to choose σ(j) and σ(j+1) = σ(j)+1 among n values)
# Wait, for a random position j in a random perm:
# P(σ(j+1) = σ(j)+1) ≈ 1/n (not exactly: depends on what's available)
# Actually for perm of [0..n-1], P(σ(j+1) = σ(j)+1) = (n-1)/(n(n-1)) = 1/n for each j
# No, it's: there are n! perms, and the number with σ(j+1)=σ(j)+1 for a fixed j
# is (n-1) × (n-2)! = (n-1)! So P = (n-1)!/n! = 1/n.
# E[A] = (n-1) × 1/n = (n-1)/n = 1-1/n

# Similarly P(unit descent at j) = 1/n (by symmetry σ → n-1-σ)
# E[D] = (n-1)/n

for n in range(3, 10):
    print(f"  n={n}: E[A]=(n-1)/n = {(n-1)/n:.6f}")

print("\nDone!")
