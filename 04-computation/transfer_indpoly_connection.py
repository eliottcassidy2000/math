#!/usr/bin/env python3
"""
Can the transfer matrix M encode the FULL independence polynomial I(Omega,x)?

We know H = I(Omega,2). The transfer matrix M (at {0,1} tournament values)
gives H via trace (odd n) or off-diagonal sum (even n).

But I(Omega,x) = 1 + alpha_1*x + alpha_2*x^2 + ...
Can we recover alpha_1, alpha_2, ... from M?

Approach: The c-tournament parameter r = c/2 is a FREE parameter.
At r = 1/2 (c=1, standard tournament): H = I(Omega(T), 2).
But what about other values of r?

For a c-tournament with t_{ij} + t_{ji} = c = 2r:
  T(W; r) = total Hamiltonian path weight (polynomial in r and s_{ij})
  At specific s_{ij} values (determined by the tournament): T(r) is a polynomial in r.

KEY QUESTION: Does T(r) evaluated at different r values recover I(Omega, x)?
For example, does T(r)|_{s=s_T} = I(Omega(T), 2r) for all r?

If this were true, then T(r) as a function of r would encode the full
independence polynomial! Let's test.
"""

from itertools import permutations, combinations
from collections import defaultdict
import numpy as np

def count_paths_subset(A, verts, start=None, end=None):
    count = 0
    for p in permutations(verts):
        if start is not None and p[0] != start:
            continue
        if end is not None and p[-1] != end:
            continue
        valid = True
        for i in range(len(p)-1):
            if A[p[i]][p[i+1]] != 1:
                valid = False
                break
        if valid:
            count += 1
    return count

def T_weighted(s_vals, r_val, n):
    """Compute T(W) = total Hamiltonian path weight for c-tournament with
    given s_{ij} values and r value."""
    total = 0.0
    for p in permutations(range(n)):
        w = 1.0
        for i in range(n-1):
            a, b = p[i], p[i+1]
            t_ab = r_val + s_vals.get((a,b), 0)
            w *= t_ab
        total += w
    return total

def tournament_to_s(A):
    """Convert {0,1} tournament to s values (r = 1/2)."""
    n = len(A)
    s = {}
    for i in range(n):
        for j in range(n):
            if i != j:
                # t(i,j) = A[i][j], and t(i,j) = r + s_{ij} = 1/2 + s_{ij}
                s[(i,j)] = A[i][j] - 0.5
    return s

def odd_cycles(A):
    n = len(A)
    cycles = []
    for length in range(3, n+1, 2):
        for verts in combinations(range(n), length):
            for perm in permutations(verts):
                valid = True
                for i in range(length):
                    if A[perm[i]][perm[(i+1) % length]] != 1:
                        valid = False
                        break
                if valid:
                    min_idx = perm.index(min(perm))
                    canonical = perm[min_idx:] + perm[:min_idx]
                    if canonical not in cycles:
                        cycles.append(canonical)
    return cycles

def independence_polynomial_coeffs(A):
    n = len(A)
    oc = odd_cycles(A)
    m = len(oc)
    adj = [[False]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1, m):
            if set(oc[i]) & set(oc[j]):
                adj[i][j] = adj[j][i] = True
    alpha = defaultdict(int)
    for mask in range(2**m):
        verts = [i for i in range(m) if (mask >> i) & 1]
        indep = True
        for i in range(len(verts)):
            for j in range(i+1, len(verts)):
                if adj[verts[i]][verts[j]]:
                    indep = False
                    break
            if not indep:
                break
        if indep:
            alpha[len(verts)] += 1
    return dict(alpha)

def I_poly_eval(alphas, x):
    return sum(alphas.get(k, 0) * x**k for k in range(max(alphas.keys())+1 if alphas else 1))

def all_tournaments(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    for bits in range(2**len(edges)):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        yield [list(row) for row in A]

print("=" * 70)
print("TEST: Does T(r) = I(Omega(T), 2r) ?")
print("=" * 70)
print()
print("At r=1/2: T(1/2) = H(T) = I(Omega(T), 2*1/2) = I(Omega(T), 1)... no,")
print("H = I(Omega(T), 2), not I(Omega(T), 1).")
print("So if T(r) = I(Omega(T), f(r)) for some function f, then f(1/2) = 2.")
print()
print("Hypothesis 1: T(r) = I(Omega(T), 2r)")
print("  At r=1/2: I(Omega, 1) -- but we need H = I(Omega, 2). WRONG.")
print()
print("Hypothesis 2: T(r) relates to I(Omega, 2) differently for each r.")
print("  The s-values are FIXED (determined by tournament). Only r varies.")
print()

# Test: for a specific tournament, compute T(r) for various r values
# and compare with I(Omega(T), x) for various x.

# n=4 tournament with 2 odd cycles (alpha_1 = 2, H = 5)
A4 = [[0,1,1,1],[0,0,0,1],[0,1,0,1],[0,0,0,0]]  # transitive
s4 = tournament_to_s(A4)
n = 4
H = count_paths_subset(A4, list(range(n)))
alphas = independence_polynomial_coeffs(A4)
print(f"Transitive n=4: H={H}, alphas={alphas}")

# T(r) for various r values
r_vals = [0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0]
for rv in r_vals:
    Tr = T_weighted(s4, rv, n)
    print(f"  T({rv}) = {Tr:.4f}")

print()

# Now try a 4-tournament with more cycles
A4b = [[0,1,1,0],[0,0,1,1],[0,0,0,1],[1,0,0,0]]  # has cycles
s4b = tournament_to_s(A4b)
H4b = count_paths_subset(A4b, list(range(n)))
alphas4b = independence_polynomial_coeffs(A4b)
print(f"Tournament with cycles n=4: H={H4b}, alphas={alphas4b}")
I_poly_str = " + ".join(f"{alphas4b.get(k,0)}*x^{k}" for k in range(max(alphas4b.keys())+1))
print(f"  I(Omega, x) = {I_poly_str}")

for rv in r_vals:
    Tr = T_weighted(s4b, rv, n)
    # Try to match I(Omega, f(r))
    print(f"  T({rv}) = {Tr:.6f}")

print()
print("--- Checking T(r) structure ---")
print("T(r) with fixed s values is a polynomial in r of degree n-1.")
print("I(Omega,x) is a polynomial in x of degree <= alpha_max(Omega).")
print("These have DIFFERENT degrees in general, so T(r) != I(Omega, g(r)).")
print()

# But wait: for the TRANSFER MATRIX at general r, what about
# tr(M(r)) (odd n) or Sigma(r) (even n)?
# These are polynomials in r of degree n-2.
# At r=1/2: tr(M) = H = I(Omega, 2) for odd n.
# At r=1: tr(M) = ?? for odd n.

print("=" * 70)
print("TEST: tr(M(r)) at various r (odd n)")
print("=" * 70)

# For a c-tournament, M[a,b] is computed using t(i,j) = r + s_{ij}
# We can compute it numerically for any r

def transfer_matrix_r(s_vals, r_val, n):
    """Compute transfer matrix for c-tournament with given s and r."""
    M = np.zeros((n, n))
    vertices = list(range(n))
    for a in range(n):
        for b in range(n):
            U = [v for v in vertices if v != a and v != b]
            total = 0.0
            for k in range(len(U)+1):
                for S in combinations(U, k):
                    S_set = set(S)
                    R = [v for v in U if v not in S_set]
                    S_verts = sorted(list(S) + [a])
                    R_verts = sorted(R + [b])
                    # E_a through S_verts
                    ea = 0.0
                    for p in permutations(S_verts):
                        if p[-1] != a: continue
                        w = 1.0
                        for i in range(len(p)-1):
                            w *= r_val + s_vals.get((p[i],p[i+1]), 0)
                        ea += w
                    # B_b through R_verts
                    bb = 0.0
                    for p in permutations(R_verts):
                        if p[0] != b: continue
                        w = 1.0
                        for i in range(len(p)-1):
                            w *= r_val + s_vals.get((p[i],p[i+1]), 0)
                        bb += w
                    total += ((-1)**k) * ea * bb
            M[a][b] = total
    return M

# n=5 tournament with alpha_1 = 4 (H=9)
# Use the cyclic tournament: 0->1->2->3->4->0, plus 0->2, 1->3, 2->4
A5 = [[0]*5 for _ in range(5)]
# 0->1, 1->2, 2->3, 3->4, 4->0, 0->3, 1->4, 2->0... no, let me use a specific one
# Let me try the almost-regular tournament
# 0->1, 0->2, 0->3, 1->2, 1->4, 2->4, 3->1, 3->2, 4->3 (score 3,2,2,2,1)
A5 = [[0,1,1,1,0],[0,0,1,0,1],[0,0,0,0,1],[0,1,1,0,0],[1,0,0,1,0]]
H5 = count_paths_subset(A5, list(range(5)))
alphas5 = independence_polynomial_coeffs(A5)
s5 = tournament_to_s(A5)
print(f"\nn=5 tournament: H={H5}, alphas={alphas5}")
I5_str = " + ".join(f"{alphas5.get(k,0)}*x^{k}" for k in range(max(alphas5.keys())+1))
print(f"  I(Omega, x) = {I5_str}")
print(f"  I(Omega, 2) = {I_poly_eval(alphas5, 2)}")

print(f"\n  tr(M(r)) and T(r) for various r:")
for rv in [0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0]:
    M_r = transfer_matrix_r(s5, rv, 5)
    tr_r = np.trace(M_r)
    T_r = T_weighted(s5, rv, 5)
    sig_r = M_r.sum() - tr_r
    print(f"  r={rv}: tr(M)={tr_r:.4f}, T={T_r:.4f}, Sigma={sig_r:.4f}")

print(f"\n  I(Omega, x) for various x:")
for x in [0, 0.5, 1, 2, 3, 4]:
    print(f"  I(Omega, {x}) = {I_poly_eval(alphas5, x)}")

# KEY CHECK: does tr(M(r)) = I(Omega, 2r) at r=1/2?
# tr(M(0.5)) should = H = I(Omega, 2).
# tr(M(1)) should = I(Omega, ???).
print(f"\n  tr(M(0.5)) = {np.trace(transfer_matrix_r(s5, 0.5, 5)):.4f}, I(Omega,2) = {I_poly_eval(alphas5, 2)}")

# Check for EVEN n too
print(f"\n--- n=4 ---")
A4c = [[0,1,0,1],[0,0,1,0],[1,0,0,1],[0,1,0,0]]  # some 4-tournament
H4c = count_paths_subset(A4c, list(range(4)))
alphas4c = independence_polynomial_coeffs(A4c)
s4c = tournament_to_s(A4c)
print(f"n=4 tournament: H={H4c}, alphas={alphas4c}")

print(f"\n  Sigma(r)/2 and T(r)/r for various r:")
for rv in [0.25, 0.5, 0.75, 1.0, 1.5, 2.0]:
    M_r = transfer_matrix_r(s4c, rv, 4)
    tr_r = np.trace(M_r)
    sig_r = M_r.sum() - tr_r
    T_r = T_weighted(s4c, rv, 4)
    print(f"  r={rv}: tr(M)={tr_r:.4f}, Sigma={sig_r:.4f}, Sigma/2={sig_r/2:.4f}, T={T_r:.4f}, T/r={T_r/rv:.4f}")

print(f"\n  I(Omega, x) for various x:")
for x in [0, 0.5, 1, 2, 3, 4]:
    print(f"  I(Omega, {x}) = {I_poly_eval(alphas4c, x)}")

print("\n" + "=" * 70)
print("CONCLUSION")
print("=" * 70)
print("""
The transfer matrix M(r, s_T) is a polynomial in r with FIXED s-values.
The independence polynomial I(Omega(T), x) is a polynomial in x.

At r = 1/2:  tr(M) = H = I(Omega, 2)  [odd n]
             Sigma = 2H = 2*I(Omega,2) [even n]

But T(r) and I(Omega,x) have DIFFERENT degrees:
  - T(r) has degree n-1 in r
  - I(Omega,x) has degree alpha(Omega) in x

So T(r) != I(Omega, f(r)) for any function f in general.
The transfer matrix captures I(Omega,2) but NOT the full polynomial.
""")
