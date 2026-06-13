#!/usr/bin/env python3
"""
dk_generating_function.py — D_k as coefficients of a generating function.

KEY INSIGHT: Define F(T, x) = sum_P x^{fwd(P)} where fwd(P) = #{forward edges}.
Then D_k = [x^k contributes via binomial] -- actually:

sum_P x^{fwd(P)} = sum_P x^{fwd(P)} = F(T, x)

And sum_k D_k * y^k = sum_P sum_k C(fwd(P),k) y^k = sum_P (1+y)^{fwd(P)} = F(T, 1+y)

So the generating function for D_k is F(T, 1+y), i.e., D_k = [y^k] F(T, 1+y).

Note: F(T, 1) = n! (all permutations), F(T, x) evaluated at x=1 gives n!.

Also: H(T) = F(T, 1) restricted to fwd = n-1 ... no.
H(T) = D_{n-1} = [y^{n-1}] F(T, 1+y) = coefficient of y^{n-1} in F(T, 1+y).

And S(T) = sum_P prod (2A[i][j]-1) = sum_P (-1)^{n-1-fwd(P)} * 1^{fwd(P)} * 2^...
Wait: prod (2A-1) = prod_{edge} (2*1-1)^{fwd} * (2*0-1)^{bwd}
    = 1^{fwd} * (-1)^{bwd} = (-1)^{n-1-fwd}
Hmm no: B[i][j] = 2A[i][j]-1, so:
  if A[i][j]=1: B=1
  if A[i][j]=0: B=-1

prod B[P_i][P_{i+1}] = 1^{fwd} * (-1)^{n-1-fwd} = (-1)^{n-1-fwd}

Wait that can't be right — S would just be sum_P (-1)^{n-1-fwd(P)}.

Let me verify: if ALL edges forward, prod = 1^{n-1} = 1. fwd = n-1, (-1)^0 = 1. ✓
If one backward: prod = (-1)^1 = -1. fwd = n-2, (-1)^1 = -1. ✓

So S(T) = sum_P (-1)^{n-1-fwd(P)} = F(T, -1) * (-1)^{n-1} ... no.

S(T) = sum_P (-1)^{n-1-fwd(P)} = (-1)^{n-1} sum_P (-1)^{-fwd(P)}
     = (-1)^{n-1} sum_P (-1)^{-fwd(P)}

Hmm, (-1)^{-fwd} = (-1)^{fwd} when fwd is integer (since (-1)^{-k} = (-1)^k).

So S(T) = (-1)^{n-1} sum_P (-1)^{fwd(P)} = (-1)^{n-1} F(T, -1).

AND: H(T) = #{P : fwd(P) = n-1} = [x^{n-1}] F(T, x) ... wait.

F(T, x) = sum_P x^{fwd(P)} is a polynomial in x of degree n-1.
H(T) = coefficient of x^{n-1} in F(T, x). ✓

S(T) = (-1)^{n-1} F(T, -1). ✓ (Verified below)

Also: D_k = [y^k] F(T, 1+y).

So the FULL information is in F(T, x), a degree-(n-1) polynomial!

New connection to OCF:
H(T) = leading coefficient of F(T, x), times x^{n-1}.
H(T) = I(Omega(T), 2) by OCF.

Is there a relation between F(T, x) and I(Omega(T), something)?

Let me compute F(T, x) for several tournaments and see.

Author: opus-2026-03-07-S43b
"""
from itertools import permutations, combinations
from collections import Counter
import math

def tournament_from_bits(bits, n):
    A = [[0]*n for _ in range(n)]
    pos = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> pos) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            pos += 1
    return A

def compute_F(A, n):
    """Compute F(T, x) = sum_P x^fwd(P), returning coefficient list."""
    coeffs = [0] * n  # F(x) = sum_{k=0}^{n-1} coeffs[k] * x^k
    for P in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if A[P[i]][P[i+1]])
        coeffs[fwd] += 1
    return coeffs

def poly_eval(coeffs, x):
    return sum(c * x**k for k, c in enumerate(coeffs))

def find_odd_cycles(A, n):
    cycles = []
    for length in range(3, n+1, 2):
        for verts in combinations(range(n), length):
            for perm in permutations(verts[1:]):
                path = (verts[0],) + perm
                is_cycle = all(A[path[i]][path[(i+1) % length]] for i in range(length))
                if is_cycle:
                    normalized = min(path[i:] + path[:i] for i in range(length))
                    if normalized not in cycles:
                        cycles.append(normalized)
    return cycles

def independence_polynomial(cycles, x):
    m = len(cycles)
    if m > 20:
        return None
    cycle_sets = [set(c) for c in cycles]
    result = 0
    for mask in range(1 << m):
        bits_list = [i for i in range(m) if mask & (1 << i)]
        independent = True
        for i in range(len(bits_list)):
            for j in range(i+1, len(bits_list)):
                if cycle_sets[bits_list[i]] & cycle_sets[bits_list[j]]:
                    independent = False
                    break
            if not independent:
                break
        if independent:
            result += x ** len(bits_list)
    return result

# Compute for n=5
print("=== n=5: F(T,x) polynomials ===")
n = 5
m = n*(n-1)//2

seen = {}
for bits in range(1 << m):
    A = tournament_from_bits(bits, n)
    F = compute_F(A, n)
    key = tuple(F)
    if key not in seen:
        cycles = find_odd_cycles(A, n)
        H = F[n-1]
        S_computed = (-1)**(n-1) * poly_eval(F, -1)
        I2 = independence_polynomial(cycles, 2)
        t3 = sum(1 for c in cycles if len(c) == 3)
        seen[key] = (F, H, S_computed, I2, t3)

print(f"{'F(x) coefficients':>40} {'H':>4} {'S':>5} {'I(2)':>5} {'t3':>3} {'F(-1)':>6} {'F(0)':>5} {'F(1)':>5}")
for key in sorted(seen.keys()):
    F, H, S, I2, t3 = seen[key]
    Fm1 = poly_eval(F, -1)
    F0 = F[0]
    F1 = poly_eval(F, 1)
    print(f"{str(F):>40} {H:4d} {S:5d} {I2:5d} {t3:3d} {Fm1:6d} {F0:5d} {F1:5d}")

# Verify S = (-1)^{n-1} F(-1)
print(f"\nVerification: S = (-1)^{n-1} * F(-1)")
for bits in range(1 << m):
    A = tournament_from_bits(bits, n)
    F = compute_F(A, n)
    S_from_F = (-1)**(n-1) * poly_eval(F, -1)
    # Direct S computation
    S_direct = 0
    for P in permutations(range(n)):
        prod_b = 1
        for i in range(n-1):
            prod_b *= (2*A[P[i]][P[i+1]] - 1)
        S_direct += prod_b
    if S_from_F != S_direct:
        print(f"  MISMATCH at bits={bits}")
        break
else:
    print("  CONFIRMED for all n=5 tournaments!")

# Key observation: F(T, x) at x=2 gives...?
print(f"\n=== F(T, 2) vs I(Omega, ?) ===")
for key in sorted(seen.keys()):
    F, H, S, I2, t3 = seen[key]
    F2 = poly_eval(F, 2)
    print(f"  F(2) = {F2:6d}, H = {H:4d}, I(Omega,2) = {I2:4d}, F(2)/H = {F2/H:.4f}")

# F(T, 2) = sum_P 2^{fwd(P)}. Each perm weighted by 2^{#forward edges}.
# Compare: I(Omega, 2) = H. So F(2) != I(Omega, 2).

# What about F(T, x) / [something]?
# The ratio F(2)/H varies. Not a clean relationship.

# Let me look at the transform: if G(y) = F(1+y), then G_k = D_k.
# G(y) = sum_k D_k y^k
# G(0) = D_0 = n!
# G(1) = F(2) = sum_P 2^{fwd(P)}
# G(-1) = F(0) = #{P with fwd(P)=0} = #{P where ALL edges backward}
#        = H(T^comp) where T^comp is complement (reverse all arcs)

# Actually F(0) = #{perms with 0 forward edges} = H(T^op) where T^op is converse.
# And H(T^op) = H(T) (since reversing all arcs preserves path count by path reversal).
# Wait: that gives F[0] = H(T). But F[n-1] = H(T) too!

# F[0] = #{perms with 0 forward} = #{perms with ALL backward}
# = #{perms P where A[P_i][P_{i+1}]=0 for all i}
# = #{perms P where P_{i+1} → P_i for all i}
# = #{anti-Hamiltonian paths of T}
# = #{Hamiltonian paths of T^op}
# = H(T^op) = H(T) (by reversal)

# So F[0] = F[n-1] = H(T)! A palindrome-like property!
print(f"\n=== F[0] vs F[n-1] (both should = H) ===")
for key in sorted(seen.keys()):
    F, H, S, I2, t3 = seen[key]
    print(f"  F[0]={F[0]}, F[n-1]={F[n-1]}, H={H}, match={F[0]==F[n-1]==H}")
