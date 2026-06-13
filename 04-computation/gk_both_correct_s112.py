#!/usr/bin/env python3
"""
gk_both_correct_s112.py — Show that TM and opus g_k differ but sum identically
kind-pasteur-2026-03-15-S112

KEY FINDING: The decomposition CV^2 = sum 2*g_k(n-2k)/(n)_{2k} is NOT unique.
- Transfer matrix g_k: degree k, natural combinatorial interpretation (matching count * 2^c)
- Opus g_k: degree 3 for k>=3, more compact but no direct combinatorial meaning for k>=4
Both are valid and give the same CV^2 for all n.
"""

from fractions import Fraction
from math import factorial

def falling_factorial(n, k):
    r = Fraction(1)
    for i in range(k):
        r *= (n - i)
    return r

def transfer_gk(n):
    num_edges = n - 2
    k_max = (n - 1) // 2
    if k_max <= 0:
        return {}
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
    return {k: total[k] / 2 for k in range(1, k_max + 1) if total[k] != 0}

def opus_gk(k, m):
    if k == 1: return Fraction(m)
    if k == 2: return Fraction(m * m)
    C = {3:(2,0,1,0), 4:(10,-33,50,-24), 5:(388,-2040,3431,-1776),
         6:(69660,-380445,653748,-342960), 7:(19826270,-109486152,189674605,-100014720)}
    if k not in C: return None
    a,b,c,d = C[k]
    return Fraction(a*m**3 + b*m**2 + c*m + d, 3)

print("n=13: Individual g_k values and their contributions")
print("="*60)
n = 13
for k in range(1, 7):
    m = n - 2*k
    if m < 1: continue
    tm_g = transfer_gk(n).get(k, Fraction(0))
    op_g = opus_gk(k, m)
    denom = falling_factorial(n, 2*k)
    tm_term = 2 * tm_g / denom
    op_term = 2 * op_g / denom if op_g is not None else None

    print(f"k={k}, m={m}:")
    print(f"  TM:   g_{k}({m}) = {tm_g:>10}, term = {float(tm_term):.12e}")
    print(f"  Opus: g_{k}({m}) = {str(op_g):>10}, term = {float(op_term):.12e}" if op_g else "  Opus: N/A")

# Sum
tm_total = sum(2*g/falling_factorial(n,2*k) for k,g in transfer_gk(n).items() if n-2*k >= 1)
op_total = Fraction(0)
for k in range(1, 7):
    m = n - 2*k
    if m < 1: continue
    g = opus_gk(k, m)
    if g: op_total += 2*g/falling_factorial(n, 2*k)

print(f"\nTotal CV^2 at n=13:")
print(f"  TM:   {float(tm_total):.15f}")
print(f"  Opus: {float(op_total):.15f}")
print(f"  Match: {tm_total == op_total}")

# The key: show the redistribution
print(f"\n{'='*60}")
print("REDISTRIBUTION between k=4 and k=5:")
d8 = falling_factorial(13, 8)
d10 = falling_factorial(13, 10)

tm_4 = 2 * Fraction(225) / d8
op_4 = 2 * Fraction(217) / d8
tm_5 = 2 * Fraction(51) / d10
op_5 = 2 * Fraction(211) / d10

print(f"  k=4: TM={float(tm_4):.12e}, Opus={float(op_4):.12e}, diff={float(tm_4-op_4):.12e}")
print(f"  k=5: TM={float(tm_5):.12e}, Opus={float(op_5):.12e}, diff={float(tm_5-op_5):.12e}")
print(f"  Sum diff: {float((tm_4-op_4)+(tm_5-op_5)):.15e}")

print(f"\n{'='*60}")
print("CONCLUSION:")
print("Transfer matrix g_k has degree k (natural combinatorial meaning).")
print("Opus g_k has degree 3 for k>=3 (valid reparametrization).")
print("Both give the same CV^2. Neither is 'wrong'.")
print(f"The degree-3 universality is a property of the DECOMPOSITION,")
print(f"not of the individual matching weights.")

# Key insight: why does degree-3 reparametrization exist?
print(f"\n{'='*60}")
print("WHY DEGREE-3 REPARAMETRIZATION EXISTS:")
print("The transfer matrix M(x) is 3x3. Its characteristic polynomial")
print("p(lambda) = lambda^3 - lambda^2 - x*lambda - x constrains")
print("lambda_1(x) = 1 + 2x + O(x^2). The 3 eigenvalues generate")
print("at most 3 exponential terms in M^N, which gives at most")
print("3 independent functions of N (hence m). Any higher-degree")
print("polynomial in m can be re-expressed using these 3 functions,")
print("hence the degree-3 reparametrization exists.")
print("This is fundamentally because the Z_j process has TRIDIAGONAL")
print("covariance (nearest-neighbor interaction), requiring only 3 states.")
